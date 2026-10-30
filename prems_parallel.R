#!/usr/bin/env Rscript

# R and command-line interfaces for fold-level PReMS fitting, checkpointed
# model-size continuation and aggregation of independently fitted CV folds.
#
# Primary R interface:
#   fit <- prems_parallel_fit(...)
#   cv.fit <- prems_parallel_summarise(...)
#
# Command-line wrappers:
#   Rscript prems_parallel.R fit [options]
#   Rscript prems_parallel.R summarise [options]

`%||%` <- function(x, y) if( is.null(x) ) y else x

.wrapper_dir <- local({
    ofile <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
    if( !is.null(ofile) && nzchar(ofile) ){
        dirname(normalizePath(ofile, mustWork=FALSE))
    }else{
        full.args <- commandArgs(trailingOnly=FALSE)
        file.arg <- grep("^--file=", full.args, value=TRUE)
        if( length(file.arg)==1L ){
            dirname(normalizePath(sub("^--file=", "", file.arg[[1]]), mustWork=FALSE))
        }else{
            getwd()
        }
    }
})

.ensure_prems <- function() {
    if( exists("prems", mode="function", inherits=TRUE) &&
        exists("ModelSearchIncrease", mode="function", inherits=TRUE) )
        return(invisible(TRUE))
    prems.file <- file.path(.wrapper_dir, "prems.R")
    if( !file.exists(prems.file) )
        stop("Source prems.R first, or place it in the same directory as prems_parallel.R.")
    source(prems.file, local=.GlobalEnv)
    invisible(TRUE)
}

subset_y <- function(y, idx, family) {
    if( family=="cox" ) y[idx,,drop=FALSE] else y[idx]
}

validate_parallel_data <- function(dat, need.foldid=FALSE) {
    if( is.character(dat) && length(dat)==1L ) dat <- readRDS(dat)
    if( !is.list(dat) ) stop("data must be an analysis list or path to an RDS file.")

    required <- c("y", "x", "family")
    missing <- required[!required %in% names(dat)]
    if( length(missing) )
        stop("Analysis data are missing: ", paste(missing, collapse=", "))

    dat$x <- as.matrix(dat$x)
    if( is.null(colnames(dat$x)) || anyDuplicated(colnames(dat$x)) )
        stop("x must have unique column names.")

    dat$x.fixed <- dat$x.fixed %||% NULL
    if( !is.null(dat$x.fixed) ) dat$x.fixed <- as.matrix(dat$x.fixed)
    dat$family <- match.arg(dat$family, c("gaussian", "binomial", "cox"))

    n <- NROW(dat$y)
    if( nrow(dat$x)!=n ) stop("NROW(y) must equal nrow(x).")
    if( !is.null(dat$x.fixed) && nrow(dat$x.fixed)!=n )
        stop("NROW(y) must equal nrow(x.fixed).")

    if( need.foldid ){
        if( is.null(dat$foldid) ) stop("foldid is required for fold-level CV fitting.")
        if( length(dat$foldid)!=n ) stop("length(foldid) must equal NROW(y).")
        if( anyNA(dat$foldid) ) stop("foldid contains missing values.")
    }
    dat
}

order_for_existing_fit <- function(fit, x, x.fixed) {
    p <- match(fit$cnames, colnames(x))
    if( anyNA(p) )
        stop("Existing fit contains predictor names not found in x: ",
             paste(fit$cnames[is.na(p)], collapse=", "))
    x <- x[,p,drop=FALSE]

    if( length(fit$cnames.fixed)==0 ){
        x.fixed <- NULL
    }else{
        if( is.null(x.fixed) )
            stop("Existing fit contains fixed covariates but x.fixed is NULL.")
        if( !is.null(colnames(x.fixed)) && !is.null(fit$cnames.fixed) ){
            q <- match(fit$cnames.fixed, colnames(x.fixed))
            if( anyNA(q) )
                stop("Existing fit contains fixed covariates not found in x.fixed: ",
                     paste(fit$cnames.fixed[is.na(q)], collapse=", "))
            x.fixed <- x.fixed[,q,drop=FALSE]
        }else if( ncol(x.fixed)!=length(fit$cnames.fixed) ){
            stop("x.fixed does not have the same number of columns as the existing fit.")
        }
    }
    list(x=x, x.fixed=x.fixed)
}

estimate_parallel_tau <- function(y, x, x.fixed, family, nfolds=10,
                                  attempts=20, lasso.factor=1, seed=1,
                                  verbose=TRUE) {
    nfolds <- min(nfolds, NROW(y))
    if( nfolds<2 ) stop("At least two observations are required to estimate tau.")

    for( a in seq_len(attempts) ){
        set.seed(seed+a-1L)
        ans <- try(
            TauEst(y=y, x=x, x.fixed=x.fixed, family=family,
                   nfolds=nfolds, parallel=FALSE),
            silent=TRUE
        )
        ok <- !inherits(ans, "try-error") && !is.null(ans$tau.opt) &&
              length(ans$tau.opt)==1L && is.numeric(ans$tau.opt) &&
              is.finite(ans$tau.opt) && ans$tau.opt>0
        if( ok ){
            tau <- ans$tau.opt*lasso.factor
            if( verbose ) message("tau = ", signif(tau,6))
            return(tau)
        }
    }
    stop("Unable to estimate tau after ", attempts,
         " attempts. Supply tau explicitly or inspect the glmnet fits.")
}

# Fit either one CV training fold or the complete data. Existing out.file
# objects are continued rather than refitted from model size 1.
prems_parallel_fit <- function(data, out.file, target.k, fold=NULL,
                               initial.k=target.k, max.s=NULL,
                               max2way="all", no.cores=1, tau=NULL,
                               tau.folds=10, tau.attempts=20,
                               lasso.factor=1, standardize=TRUE,
                               seed=1, verbose=TRUE) {
    .ensure_prems()
    if( missing(out.file) || !nzchar(out.file) ) stop("out.file is required for checkpointing.")
    if( length(target.k)!=1L || !is.finite(target.k) || target.k<2 )
        stop("target.k must be an integer >= 2.")
    target.k <- as.integer(target.k)
    initial.k <- as.integer(min(initial.k, target.k))
    if( initial.k<2 ) stop("initial.k must be >= 2.")
    if( no.cores<1 ) stop("no.cores must be >= 1.")
    if( !identical(max2way, "all") ) max2way <- as.integer(max2way)

    dat <- validate_parallel_data(data, need.foldid=!is.null(fold))
    n <- NROW(dat$y)

    if( is.null(fold) ){
        idx <- seq_len(n)
        label <- "full data"
    }else{
        if( !fold %in% unique(dat$foldid) ) stop("Fold not present in foldid: ", fold)
        idx <- which(dat$foldid!=fold)
        label <- paste0("CV fold ", fold, " training data")
    }

    y <- subset_y(dat$y, idx, dat$family)
    x <- dat$x[idx,,drop=FALSE]
    x.fixed <- if( is.null(dat$x.fixed) ) NULL else dat$x.fixed[idx,,drop=FALSE]

    if( file.exists(out.file) ){
        fit <- readRDS(out.file)
        if( !identical(fit$family, dat$family) )
            stop("Existing fit family does not match analysis data.")
        ordered <- order_for_existing_fit(fit, x, x.fixed)
        x <- ordered$x
        x.fixed <- ordered$x.fixed
        current.k <- length(fit$fitted.models)
        if( verbose ) message("Loaded ", out.file, " through model size ", current.k, ".")
    }else{
        max.s.initial <- max.s %||% 50L
        tau.use <- tau
        if( is.null(tau.use) ){
            tau.use <- estimate_parallel_tau(
                y=y, x=x, x.fixed=x.fixed, family=dat$family,
                nfolds=tau.folds, attempts=tau.attempts,
                lasso.factor=lasso.factor, seed=seed, verbose=verbose
            )
        }
        if( verbose ) message("Starting PReMS on ", label,
                              " through model size ", initial.k, ".")
        set.seed(seed)
        fit <- prems(
            y=y, x=x, x.fixed=x.fixed, family=dat$family,
            tau=tau.use, k.max=initial.k, max.s=max.s.initial,
            max2way=max2way, no.cores=no.cores,
            standardize=standardize, verbose=verbose
        )
        dir.create(dirname(out.file), recursive=TRUE, showWarnings=FALSE)
        saveRDS(fit, out.file)
        current.k <- length(fit$fitted.models)
    }

    usable <- sum(fit$sd!=0)
    if( target.k>usable )
        stop("target.k exceeds the number of non-constant candidate predictors (", usable, ").")

    if( current.k<target.k ){
        if( verbose ) message("Extending search from model size ", current.k,
                              " through model size ", target.k, ".")
        checkpoint <- function(fitted.models, k) {
            saveRDS(fitted.models, out.file)
            if( verbose ) message("Model size ", k, " checkpointed.")
            invisible(NULL)
        }
        fit <- ModelSearchIncrease(
            fitted.models=fit, y=y, x=x, x.fixed=x.fixed,
            k.max=target.k, no.cores=no.cores, max.s=max.s,
            verbose=FALSE, checkpoint=checkpoint
        )
        current.k <- length(fit$fitted.models)
    }
    invisible(fit)
}

resolve_fit_file <- function(template, fold) {
    if( !grepl("{fold}", template, fixed=TRUE) )
        stop("fit.template must contain {fold}.")
    gsub("{fold}", as.character(fold), template, fixed=TRUE)
}

# Aggregate independently fitted folds into the same summary structure used by
# cv.prems(), including the fold-by-model-size score matrix required by
# prems.cv.stop().
prems_parallel_summarise <- function(data, fit.template, out.file=NULL,
                                     k.min=1, k.max=NULL,
                                     criteria="ML", verbose=TRUE) {
    .ensure_prems()
    criteria <- match.arg(criteria, c("ML", "aic"))
    dat <- validate_parallel_data(data, need.foldid=TRUE)
    folds <- sort(unique(dat$foldid))
    fit.files <- vapply(folds, function(f) resolve_fit_file(fit.template,f), character(1))
    missing <- fit.files[!file.exists(fit.files)]
    if( length(missing) ) stop("Missing fold fit(s): ", paste(missing, collapse=", "))

    fits <- lapply(fit.files, readRDS)
    common.k <- min(vapply(fits, function(z) length(z$fitted.models), integer(1)))
    k.max <- k.max %||% common.k
    if( k.max>common.k )
        stop("k.max exceeds the largest model size common to all fold fits (", common.k, ").")
    if( k.min<1 || k.min>k.max ) stop("Require 1 <= k.min <= k.max.")

    sizes <- k.min:k.max
    pwll <- matrix(NA_real_, nrow=length(folds), ncol=length(sizes),
                   dimnames=list(paste0("fold",folds), as.character(sizes)))
    selected.coef <- vector("list", k.max)
    for( k in sizes )
        selected.coef[[k]] <- matrix(NA_character_, nrow=length(folds), ncol=k,
                                     dimnames=list(paste0("fold",folds),NULL))

    for( ii in seq_along(folds) ){
        fold <- folds[[ii]]
        fit <- fits[[ii]]
        if( !identical(fit$family, dat$family) )
            stop("Family mismatch in ", fit.files[[ii]], ".")

        test <- which(dat$foldid==fold)
        y.test <- subset_y(dat$y, test, dat$family)
        x.test <- dat$x[test,,drop=FALSE]
        x.fixed.test <- if( is.null(dat$x.fixed) ) NULL else dat$x.fixed[test,,drop=FALSE]

        for( jj in seq_along(sizes) ){
            k <- sizes[[jj]]
            pred <- as.numeric(predict.prems(
                fit, newx=x.test, newx.fixed=x.fixed.test,
                size=k, criteria=criteria
            ))

            if( dat$family=="binomial" ){
                p <- pmin(pmax(pred,.Machine$double.eps),1-.Machine$double.eps)
                pwll[ii,jj] <- sum(y.test*log(p)+(1-y.test)*log(1-p))
            }else if( dat$family=="gaussian" ){
                pwll[ii,jj] <- sum(-(y.test-pred)^2)
            }else{
                pwll[ii,jj] <- cox_pl_left_trunc(y.test,pred)
            }

            coef.names <- names(getModelFit(fit,size=k,rank=1,criteria=criteria)$beta)
            selected <- setdiff(coef.names,c("I",fit$cnames.fixed))
            if( length(selected)!=k )
                stop("Expected ", k, " selected predictors in fold ", fold,
                     ", model size ", k, "; found ", length(selected), ".")
            selected.coef[[k]][ii,] <- selected
        }
        if( verbose ) message("Scored fold ", fold, ".")
    }

    cvm <- colMeans(pwll)
    cvsd <- apply(pwll,2,stats::sd)/sqrt(nrow(pwll))
    names(cvm) <- names(cvsd) <- as.character(sizes)
    opt.idx <- prems.optim(cvm,cvsd)
    selected.sizes <- sizes[opt.idx]

    ret <- list(
        best=selected.sizes[[1]], one.se=selected.sizes[[2]],
        cvm=cvm, cvsd=cvsd, selected.coef=selected.coef,
        pwll=pwll, foldid=dat$foldid, fit.files=fit.files,
        family=dat$family, criteria=criteria
    )
    class(ret) <- c("cv.prems.parallel","list")

    if( !is.null(out.file) ){
        dir.create(dirname(out.file), recursive=TRUE, showWarnings=FALSE)
        saveRDS(ret,out.file)
        if( verbose ) message("Saved CV summary to ",out.file,".")
    }
    ret
}

# --------------------------- CLI wrappers ---------------------------

parse_options <- function(x) {
    out <- list(); i <- 1L
    while( i<=length(x) ){
        key <- x[[i]]
        if( !startsWith(key,"--") ) stop("Unexpected argument: ",key)
        key <- sub("^--","",key)
        if( i==length(x) || startsWith(x[[i+1L]],"--") ){
            out[[key]] <- TRUE; i <- i+1L
        }else{
            out[[key]] <- x[[i+1L]]; i <- i+2L
        }
    }
    out
}

as_bool <- function(x, default=FALSE) {
    if( is.null(x) ) return(default)
    if( is.logical(x) ) return(x)
    z <- tolower(as.character(x))
    if( z %in% c("true","t","1","yes","y") ) return(TRUE)
    if( z %in% c("false","f","0","no","n") ) return(FALSE)
    stop("Cannot interpret logical value: ",x)
}
as_int <- function(x, default=NULL) { if(is.null(x)) return(default); z<-as.integer(x); if(is.na(z)) stop("Expected integer: ",x); z }
as_num <- function(x, default=NULL) { if(is.null(x)) return(default); z<-as.numeric(x); if(is.na(z)) stop("Expected numeric value: ",x); z }

usage <- function(status=0L) {
    cat(
"PReMS parallel/checkpoint wrapper\n\n",
"Primary R interface:\n",
"  prems_parallel_fit(...)\n",
"  prems_parallel_summarise(...)\n\n",
"Command-line use:\n",
"  Rscript prems_parallel.R fit [options]\n",
"  Rscript prems_parallel.R summarise [options]\n\n",
"fit: --data FILE --out FILE --target-k INTEGER [--fold INTEGER]\n",
"     [--initial-k INTEGER] [--max-s INTEGER] [--max2way all|INTEGER]\n",
"     [--cores INTEGER] [--tau NUMERIC] [--tau-folds INTEGER]\n",
"     [--tau-attempts INTEGER] [--lasso-factor NUMERIC]\n",
"     [--standardize TRUE|FALSE] [--seed INTEGER]\n\n",
"summarise: --data FILE --fit-template STRING --out FILE\n",
"           [--k-min INTEGER] [--k-max INTEGER] [--criteria ML|aic]\n",
sep="")
    quit(save="no",status=status)
}

if( sys.nframe()==0L ){
    args <- commandArgs(trailingOnly=TRUE)
    if( !length(args) || args[[1]] %in% c("-h","--help","help") ) usage(0L)
    command <- args[[1]]
    opt <- parse_options(args[-1])

    if( command=="fit" ){
        data.file <- opt$data %||% stop("--data is required.")
        out.file <- opt$out %||% stop("--out is required.")
        target.k <- as_int(opt[["target-k"]]) %||% stop("--target-k is required.")
        max2way.raw <- opt$max2way %||% "all"
        max2way <- if(tolower(max2way.raw)=="all") "all" else as_int(max2way.raw)
        prems_parallel_fit(
            data=data.file, out.file=out.file, target.k=target.k,
            fold=as_int(opt$fold,NULL),
            initial.k=as_int(opt[["initial-k"]],target.k),
            max.s=as_int(opt[["max-s"]],NULL), max2way=max2way,
            no.cores=as_int(opt$cores,1L), tau=as_num(opt$tau,NULL),
            tau.folds=as_int(opt[["tau-folds"]],10L),
            tau.attempts=as_int(opt[["tau-attempts"]],20L),
            lasso.factor=as_num(opt[["lasso-factor"]],1),
            standardize=as_bool(opt$standardize,TRUE),
            seed=as_int(opt$seed,1L), verbose=as_bool(opt$verbose,TRUE)
        )
    }else if( command=="summarise" ){
        prems_parallel_summarise(
            data=opt$data %||% stop("--data is required."),
            fit.template=opt[["fit-template"]] %||% stop("--fit-template is required."),
            out.file=opt$out %||% stop("--out is required."),
            k.min=as_int(opt[["k-min"]],1L),
            k.max=as_int(opt[["k-max"]],NULL),
            criteria=opt$criteria %||% "ML",
            verbose=as_bool(opt$verbose,TRUE)
        )
    }else{
        stop("Unknown command: ",command,". Use fit or summarise.")
    }
}
