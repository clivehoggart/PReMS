#!/usr/bin/env Rscript

# Command-line wrapper for parallel PReMS cross-validation and checkpointed
# model-search continuation.
#
# Commands:
#   fit        Fit a full-data model or one CV training fold. If --out already
#              contains a PReMS fit, continue it to --target-k using
#              ModelSearchIncrease().
#   summarise  Combine fold-specific fits into a cv.prems()-compatible summary.
#
# The analysis RDS supplied with --data should be a list containing:
#   y        outcome vector or survival::Surv object
#   x        numeric candidate-predictor matrix with column names
#   x.fixed  optional numeric matrix of fixed covariates (or NULL)
#   family   "gaussian", "binomial" or "cox"
#   foldid   required only for parallel CV; integer/vector fold assignments

`%||%` <- function(x, y) if( is.null(x) ) y else x

usage <- function(status=0) {
    cat(
"PReMS parallel wrapper\n\n",
"Usage:\n",
"  Rscript prems_parallel.R fit [options]\n",
"  Rscript prems_parallel.R summarise [options]\n\n",
"Common options:\n",
"  --data FILE            Analysis-data RDS (required)\n",
"  --verbose TRUE|FALSE   Print progress [TRUE]\n\n",
"fit options:\n",
"  --out FILE             Output PReMS fit RDS (required)\n",
"  --fold INTEGER         CV fold to leave out; omit for full-data fit\n",
"  --initial-k INTEGER    Initial search maximum for a new fit [target-k]\n",
"  --target-k INTEGER     Target maximum model size (required)\n",
"  --max-s INTEGER        Models retained per size [50 for new fits;\n",
"                         current retained breadth when continuing]\n",
"  --max2way all|INTEGER  Initial two-predictor search [all]\n",
"  --cores INTEGER        Cores used within each PReMS search [1]\n",
"  --tau NUMERIC          Fixed tau; omit to estimate within training data\n",
"  --tau-folds INTEGER    glmnet folds used to estimate tau [10]\n",
"  --tau-attempts INTEGER Maximum attempts to estimate tau [20]\n",
"  --lasso-factor NUMERIC Multiplier applied to estimated tau [1]\n",
"  --standardize TRUE|FALSE  Standardise predictors [TRUE]\n",
"  --seed INTEGER         Reproducibility seed [1]\n\n",
"summarise options:\n",
"  --fit-template STRING  Fold fit path containing {fold}, e.g.\n",
"                         cv/fold{fold}.rds (required)\n",
"  --out FILE             Output CV summary RDS (required)\n",
"  --k-min INTEGER        Smallest model size to evaluate [1]\n",
"  --k-max INTEGER        Largest size; default is largest size common\n",
"                         to all fold fits\n",
"  --criteria ML|aic      Within-size ranking criterion [ML]\n\n",
sep="")
    quit(save="no", status=status)
}

parse_options <- function(x) {
    out <- list()
    i <- 1L
    while( i <= length(x) ) {
        key <- x[[i]]
        if( !startsWith(key, "--") ) {
            stop("Unexpected argument: ", key)
        }
        key <- sub("^--", "", key)
        if( i == length(x) || startsWith(x[[i+1L]], "--") ) {
            out[[key]] <- TRUE
            i <- i + 1L
        }else{
            out[[key]] <- x[[i+1L]]
            i <- i + 2L
        }
    }
    out
}

as_bool <- function(x, default=FALSE) {
    if( is.null(x) ) return(default)
    if( is.logical(x) ) return(x)
    x <- tolower(as.character(x))
    if( x %in% c("true", "t", "1", "yes", "y") ) return(TRUE)
    if( x %in% c("false", "f", "0", "no", "n") ) return(FALSE)
    stop("Cannot interpret logical value: ", x)
}

as_int <- function(x, default=NULL) {
    if( is.null(x) ) return(default)
    ans <- suppressWarnings(as.integer(x))
    if( is.na(ans) ) stop("Expected integer, got: ", x)
    ans
}

as_num <- function(x, default=NULL) {
    if( is.null(x) ) return(default)
    ans <- suppressWarnings(as.numeric(x))
    if( is.na(ans) ) stop("Expected numeric value, got: ", x)
    ans
}

subset_y <- function(y, idx, family) {
    if( family == "cox" ) y[idx,,drop=FALSE] else y[idx]
}

validate_data <- function(dat, need_foldid=FALSE) {
    required <- c("y", "x", "family")
    missing <- required[!required %in% names(dat)]
    if( length(missing) ) {
        stop("Analysis RDS is missing: ", paste(missing, collapse=", "))
    }
    dat$x <- as.matrix(dat$x)
    if( is.null(colnames(dat$x)) || anyDuplicated(colnames(dat$x)) ) {
        stop("x must have unique column names.")
    }
    dat$x.fixed <- dat$x.fixed %||% NULL
    if( !is.null(dat$x.fixed) ) dat$x.fixed <- as.matrix(dat$x.fixed)
    dat$family <- match.arg(dat$family, c("gaussian", "binomial", "cox"))

    n <- NROW(dat$y)
    if( nrow(dat$x) != n ) stop("NROW(y) must equal nrow(x).")
    if( !is.null(dat$x.fixed) && nrow(dat$x.fixed) != n ) {
        stop("NROW(y) must equal nrow(x.fixed).")
    }
    if( need_foldid ) {
        if( is.null(dat$foldid) ) stop("Analysis RDS must contain foldid for CV.")
        if( length(dat$foldid) != n ) stop("length(foldid) must equal NROW(y).")
        if( anyNA(dat$foldid) ) stop("foldid contains missing values.")
    }
    dat
}

order_for_existing_fit <- function(fit, x, x.fixed) {
    p <- match(fit$cnames, colnames(x))
    if( anyNA(p) ) {
        stop("Existing fit contains predictor names not found in x: ",
             paste(fit$cnames[is.na(p)], collapse=", "))
    }
    x <- x[,p,drop=FALSE]

    if( length(fit$cnames.fixed) == 0 ) {
        x.fixed <- NULL
    }else{
        if( is.null(x.fixed) ) stop("Existing fit contains fixed covariates but x.fixed is NULL.")
        if( !is.null(colnames(x.fixed)) && !is.null(fit$cnames.fixed) ) {
            q <- match(fit$cnames.fixed, colnames(x.fixed))
            if( anyNA(q) ) {
                stop("Existing fit contains fixed covariates not found in x.fixed: ",
                     paste(fit$cnames.fixed[is.na(q)], collapse=", "))
            }
            x.fixed <- x.fixed[,q,drop=FALSE]
        }else if( ncol(x.fixed) != length(fit$cnames.fixed) ) {
            stop("x.fixed does not have the same number of columns as the existing fit.")
        }
    }
    list(x=x, x.fixed=x.fixed)
}

estimate_tau <- function(y, x, x.fixed, family, nfolds, attempts,
                         lasso.factor, seed, verbose) {
    nfolds <- min(nfolds, NROW(y))
    if( nfolds < 2 ) stop("At least two observations are required to estimate tau.")

    for( a in seq_len(attempts) ) {
        set.seed(seed + a - 1L)
        ans <- try(
            TauEst(
                y=y, x=x, x.fixed=x.fixed, family=family,
                nfolds=nfolds, parallel=FALSE
            ),
            silent=TRUE
        )
        ok <- !inherits(ans, "try-error") &&
              !is.null(ans$tau.opt) && length(ans$tau.opt) == 1L &&
              is.numeric(ans$tau.opt) && is.finite(ans$tau.opt) && ans$tau.opt > 0
        if( ok ) {
            tau <- ans$tau.opt * lasso.factor
            if( verbose ) message("tau = ", signif(tau, 6))
            return(tau)
        }
    }
    stop("Unable to estimate tau after ", attempts,
         " attempts. Supply --tau explicitly or inspect the glmnet fits.")
}

run_fit <- function(opt) {
    data.file <- opt$data %||% stop("--data is required.")
    out.file <- opt$out %||% stop("--out is required.")
    target.k <- as_int(opt[["target-k"]])
    if( is.null(target.k) ) stop("--target-k is required.")
    if( target.k < 2 ) stop("--target-k must be at least 2.")

    cores <- as_int(opt$cores, 1L)
    if( cores < 1 ) stop("--cores must be >= 1.")
    max.s.opt <- as_int(opt[["max-s"]], NULL)
    max2way.raw <- opt$max2way %||% "all"
    max2way <- if( tolower(max2way.raw) == "all" ) "all" else as_int(max2way.raw)
    standardize <- as_bool(opt$standardize, TRUE)
    verbose <- as_bool(opt$verbose, TRUE)
    seed <- as_int(opt$seed, 1L)
    tau.user <- as_num(opt$tau, NULL)
    tau.folds <- as_int(opt[["tau-folds"]], 10L)
    tau.attempts <- as_int(opt[["tau-attempts"]], 20L)
    lasso.factor <- as_num(opt[["lasso-factor"]], 1)

    dat <- validate_data(readRDS(data.file), need_foldid=!is.null(opt$fold))
    n <- NROW(dat$y)

    if( is.null(opt$fold) ) {
        idx <- seq_len(n)
        label <- "full data"
    }else{
        fold <- as_int(opt$fold)
        if( !fold %in% unique(dat$foldid) ) stop("Fold not present in foldid: ", fold)
        idx <- which(dat$foldid != fold)
        label <- paste0("CV fold ", fold, " training data")
    }

    y <- subset_y(dat$y, idx, dat$family)
    x <- dat$x[idx,,drop=FALSE]
    x.fixed <- if( is.null(dat$x.fixed) ) NULL else dat$x.fixed[idx,,drop=FALSE]

    if( file.exists(out.file) ) {
        fit <- readRDS(out.file)
        if( !identical(fit$family, dat$family) ) {
            stop("Existing fit family is ", fit$family,
                 " but analysis data family is ", dat$family, ".")
        }
        ordered <- order_for_existing_fit(fit, x, x.fixed)
        x <- ordered$x
        x.fixed <- ordered$x.fixed
        current.k <- length(fit$fitted.models)
        if( verbose ) message("Loaded ", out.file, " at model size ", current.k, ".")
    }else{
        initial.k <- as_int(opt[["initial-k"]], target.k)
        initial.k <- min(initial.k, target.k)
        if( initial.k < 2 ) stop("--initial-k must be at least 2.")
        max.s.initial <- max.s.opt %||% 50L

        tau <- tau.user
        if( is.null(tau) ) {
            tau <- estimate_tau(
                y=y, x=x, x.fixed=x.fixed, family=dat$family,
                nfolds=tau.folds, attempts=tau.attempts,
                lasso.factor=lasso.factor, seed=seed,
                verbose=verbose
            )
        }

        if( verbose ) {
            message("Starting PReMS on ", label,
                    " through model size ", initial.k, ".")
        }
        set.seed(seed)
        fit <- prems(
            y=y, x=x, x.fixed=x.fixed, family=dat$family,
            tau=tau, k.max=initial.k, max.s=max.s.initial,
            max2way=max2way, no.cores=cores,
            standardize=standardize, verbose=verbose
        )
        dir.create(dirname(out.file), recursive=TRUE, showWarnings=FALSE)
        saveRDS(fit, out.file)
        current.k <- length(fit$fitted.models)
        if( verbose ) message("Saved model sizes 1-", current.k, " to ", out.file, ".")
    }

    usable <- sum(fit$sd != 0)
    if( target.k > usable ) {
        stop("--target-k = ", target.k,
             " exceeds the number of non-constant candidate predictors (", usable, ").")
    }

    while( current.k < target.k ) {
        if( verbose ) message("Extending search to model size ", current.k + 1L, ".")
        fit <- ModelSearchIncrease(
            fitted.models=fit,
            y=y,
            x=x,
            x.fixed=x.fixed,
            no.cores=cores,
            max.s=max.s.opt
        )
        current.k <- length(fit$fitted.models)
        saveRDS(fit, out.file)
        if( verbose ) message("Model size ", current.k, " complete and checkpointed.")
    }

    if( verbose && current.k >= target.k ) {
        message("Fit complete through model size ", current.k, ": ", out.file)
    }
    invisible(fit)
}

resolve_fit_file <- function(template, fold) {
    if( !grepl("{fold}", template, fixed=TRUE) ) {
        stop("--fit-template must contain {fold}.")
    }
    gsub("{fold}", as.character(fold), template, fixed=TRUE)
}

run_summarise <- function(opt) {
    data.file <- opt$data %||% stop("--data is required.")
    out.file <- opt$out %||% stop("--out is required.")
    template <- opt[["fit-template"]] %||% stop("--fit-template is required.")
    criteria <- opt$criteria %||% "ML"
    if( !criteria %in% c("ML", "aic") ) stop("--criteria must be ML or aic.")
    k.min <- as_int(opt[["k-min"]], 1L)
    k.max.user <- as_int(opt[["k-max"]], NULL)
    verbose <- as_bool(opt$verbose, TRUE)

    dat <- validate_data(readRDS(data.file), need_foldid=TRUE)
    folds <- sort(unique(dat$foldid))
    fit.files <- vapply(folds, function(f) resolve_fit_file(template, f), character(1))
    missing <- fit.files[!file.exists(fit.files)]
    if( length(missing) ) stop("Missing fold fit(s): ", paste(missing, collapse=", "))

    fits <- lapply(fit.files, readRDS)
    common.k <- min(vapply(fits, function(z) length(z$fitted.models), integer(1)))
    k.max <- k.max.user %||% common.k
    if( k.max > common.k ) {
        stop("Requested k-max = ", k.max,
             " but the largest model size common to all folds is ", common.k, ".")
    }
    if( k.min < 1 || k.min > k.max ) stop("Require 1 <= k-min <= k-max.")

    sizes <- k.min:k.max
    pwll <- matrix(
        NA_real_, nrow=length(folds), ncol=length(sizes),
        dimnames=list(paste0("fold", folds), as.character(sizes))
    )
    selected.coef <- vector("list", k.max)
    for( k in sizes ) {
        selected.coef[[k]] <- matrix(
            NA_character_, nrow=length(folds), ncol=k,
            dimnames=list(paste0("fold", folds), NULL)
        )
    }

    for( ii in seq_along(folds) ) {
        fold <- folds[[ii]]
        fit <- fits[[ii]]
        if( !identical(fit$family, dat$family) ) {
            stop("Family mismatch in ", fit.files[[ii]], ".")
        }
        test <- which(dat$foldid == fold)
        y.test <- subset_y(dat$y, test, dat$family)
        x.test <- dat$x[test,,drop=FALSE]
        x.fixed.test <- if( is.null(dat$x.fixed) ) NULL else dat$x.fixed[test,,drop=FALSE]

        for( jj in seq_along(sizes) ) {
            k <- sizes[[jj]]
            pred <- predict.prems(
                fit, newx=x.test, newx.fixed=x.fixed.test,
                size=k, criteria=criteria
            )
            pred <- as.numeric(pred)

            if( dat$family == "binomial" ) {
                p <- pmin(pmax(pred, .Machine$double.eps), 1-.Machine$double.eps)
                pwll[ii,jj] <- sum(y.test*log(p) + (1-y.test)*log(1-p))
            }else if( dat$family == "gaussian" ) {
                pwll[ii,jj] <- sum(-(y.test-pred)^2)
            }else{
                pwll[ii,jj] <- cox_pl_left_trunc(y.test, pred)
            }

            coef.names <- names(getModelFit(fit, size=k, rank=1, criteria=criteria)$beta)
            selected <- setdiff(coef.names, c("I", fit$cnames.fixed))
            if( length(selected) != k ) {
                stop("Expected ", k, " selected predictors in fold ", fold,
                     ", model size ", k, "; found ", length(selected), ".")
            }
            selected.coef[[k]][ii,] <- selected
        }
        if( verbose ) message("Scored fold ", fold, ".")
    }

    cvm <- colMeans(pwll)
    cvsd <- apply(pwll, 2, stats::sd) / sqrt(nrow(pwll))
    names(cvm) <- names(cvsd) <- as.character(sizes)
    opt.idx <- prems.optim(cvm, cvsd)
    selected.sizes <- sizes[opt.idx]

    ret <- list(
        best=selected.sizes[[1]],
        one.se=selected.sizes[[2]],
        cvm=cvm,
        cvsd=cvsd,
        selected.coef=selected.coef,
        pwll=pwll,
        foldid=dat$foldid,
        fit.files=fit.files,
        family=dat$family,
        criteria=criteria
    )
    class(ret) <- c("cv.prems.parallel", "list")

    dir.create(dirname(out.file), recursive=TRUE, showWarnings=FALSE)
    saveRDS(ret, out.file)
    if( verbose ) {
        message("Saved CV summary to ", out.file, ".")
        message("Best model size: ", ret$best,
                "; one-SE model size: ", ret$one.se, ".")
    }
    invisible(ret)
}

args <- commandArgs(trailingOnly=TRUE)
if( length(args) == 0 || args[[1]] %in% c("-h", "--help", "help") ) usage(0)
command <- args[[1]]
opt <- parse_options(args[-1])

# Locate prems.R relative to this wrapper rather than relative to the
# user's working directory. This keeps prems_parallel.R and prems.R as a
# self-contained pair that can be called from anywhere.
full.args <- commandArgs(trailingOnly=FALSE)
file.arg <- grep("^--file=", full.args, value=TRUE)
if( length(file.arg) != 1L ) {
    stop("Unable to determine the location of prems_parallel.R. Run this script with Rscript.")
}
script.file <- sub("^--file=", "", file.arg[[1]])
script.dir <- dirname(normalizePath(script.file, mustWork=TRUE))
prems.file <- file.path(script.dir, "prems.R")
if( !file.exists(prems.file) ) {
    stop("Cannot find prems.R in the same directory as prems_parallel.R: ", script.dir)
}
source(prems.file)

if( command == "fit" ) {
    run_fit(opt)
}else if( command == "summarise" ) {
    run_summarise(opt)
}else{
    stop("Unknown command: ", command, ". Use fit or summarise.")
}
