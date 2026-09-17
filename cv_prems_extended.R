cv.prems <- function( y, x, x.fixed=NULL, no.cores=10, k.min=1, k.max, tau.i=NULL,
                      max.s=50, max2way='all', standardize=TRUE, nfolds=NULL, foldid=NULL,
                      n.waic=100, lasso.factor=1,
                      criteria='ML', fit='mode', family='binomial', verbose=TRUE,
                      fit.dir=NULL, fit.action=c('restart','continue') ){

    fit.action <- match.arg(fit.action)
    use.saved.fits <- !is.null(fit.dir)

    ## Set up folds. When continuing, reuse the original fold allocation.
    if( use.saved.fits && fit.action=='continue' ){

        if( !dir.exists(fit.dir) )
            stop('fit.dir does not exist: ', fit.dir)

        fold.file <- file.path(fit.dir, 'foldid.rds')

        if( file.exists(fold.file) ){
            saved.foldid <- readRDS(fold.file)

            if( !is.null(foldid) &&
                !identical(as.integer(foldid), as.integer(saved.foldid)) )
                stop('Supplied foldid differs from foldid saved in fit.dir.')

            foldid <- saved.foldid

        }else if( is.null(foldid) ){
            stop('No foldid.rds found in fit.dir. Supply the original foldid or restart the cross-validation search.')
        }else{
            ## Allows existing saved fits to be adopted if the original
            ## fold allocation is supplied explicitly.
            saveRDS(foldid, fold.file)
        }

        saved.nfolds <- length(unique(foldid))

        if( !is.null(nfolds) && nfolds != saved.nfolds )
            stop('nfolds differs from the saved fold allocation.')

        nfolds <- saved.nfolds

    }else{

        if( is.null(foldid) & is.null(nfolds) ){
            nfolds <- NROW(y)
            foldid <- 1:nfolds
        }

        if( !is.null(foldid) & is.null(nfolds) )
            nfolds <- length(unique(foldid))

        if( is.null(foldid) ){
            yy <- y

            if( family=='cox' )
                yy <- y[,3]

            if( family=='cox' | family=='binomial' ){
                foldid <- make.folds2( yy, folds=nfolds )
                if( verbose )
                    print( table( yy, foldid ) )
            }

            if( family=='gaussian' ){
                foldid <- make.folds.continuous( NROW(y), nfolds )
                if( verbose )
                    print( table( foldid ) )
            }
        }

        if( use.saved.fits ){
            if( !dir.exists(fit.dir) )
                dir.create(fit.dir, recursive=TRUE)

            ## Restart replaces fold fits created by this function.
            old.files <- list.files(
                fit.dir,
                pattern='^prems_fold[0-9]+\\.rds$',
                full.names=TRUE
            )
            if( length(old.files)>0 )
                unlink(old.files)

            saveRDS(foldid, file.path(fit.dir, 'foldid.rds'))
        }
    }

    if( length(foldid) != NROW(y) )
        stop('foldid must contain one value per observation.')

    if( anyNA(foldid) )
        stop('foldid contains missing values.')

    if( !identical(sort(unique(as.integer(foldid))), seq_len(nfolds)) )
        stop('foldid values must be the integers 1,...,nfolds.')

    pwll <- matrix( nrow=nfolds, ncol=(k.max-k.min+1) )

    if( verbose )
        print(paste(nfolds, 'fold cross-validation'))

    selected.coef <- vector('list', k.max)
    for( k in k.min:k.max )
        selected.coef[[k]] <- matrix(nrow=nfolds, ncol=k)

    for( i in 1:nfolds ){

        train <- which(foldid != i)
        test  <- which(foldid == i)

        x.fixed.train <- if( is.null(x.fixed) )
            NULL else x.fixed[train,,drop=FALSE]

        x.fixed.test <- if( is.null(x.fixed) )
            NULL else x.fixed[test,,drop=FALSE]

        outfile <- if( use.saved.fits )
            file.path(fit.dir, paste0('prems_fold', i, '.rds')) else NULL

        if( use.saved.fits && fit.action=='continue' ){

            if( !file.exists(outfile) )
                stop('No saved PReMS fit found for fold ', i, ': ', outfile)

            my.fit <- readRDS(outfile)

            if( !identical(my.fit$family, family) )
                stop('Saved PReMS fit for fold ', i,
                     " was fitted using family='", my.fit$family,
                     "', not family='", family, "'.")

            current.k <- length(my.fit$fitted.models)

            if( verbose )
                print(paste('Fold', i, ': continuing from model size',
                            current.k, 'to', k.max))

            while( current.k < k.max ){

                my.fit <- ModelSearchIncrease(
                    fitted.models=my.fit,
                    y=y[train],
                    x=x[train,,drop=FALSE],
                    x.fixed=x.fixed.train,
                    no.cores=no.cores,
                    max.s=max.s
                )

                current.k <- length(my.fit$fitted.models)

                ## Checkpoint after every added model size.
                saveRDS(my.fit, outfile)

                if( verbose )
                    print(paste('Fold', i, ': model size',
                                current.k, 'complete.'))
            }

        }else{

            if( is.null(tau.i) ){

                tauest <- NULL
                attempt <- 1
                max.attempts <- 20
                tau.nfolds <- 10

                while( attempt <= max.attempts ){

                    tauest.try <- try(
                        TauEst(
                            y=y[train],
                            x=x[train,,drop=FALSE],
                            x.fixed=x.fixed.train,
                            family=family,
                            nfolds=tau.nfolds,
                            parallel=TRUE
                        ),
                        silent=TRUE
                    )

                    ok <- !inherits(tauest.try, 'try-error') &&
                          !is.null(tauest.try$tau.opt) &&
                          length(tauest.try$tau.opt)==1 &&
                          is.numeric(tauest.try$tau.opt) &&
                          is.finite(tauest.try$tau.opt) &&
                          tauest.try$tau.opt > 0

                    if( ok ){
                        tauest <- tauest.try
                        break
                    }

                    attempt <- attempt + 1
                }

                if( is.null(tauest) )
                    stop('TauEst failed for fold ', i, ' after ',
                         max.attempts, ' attempts.')

                tau <- tauest$tau.opt * lasso.factor

                if( verbose )
                    print(paste0('tau=', tau))

            }else{
                tau <- tau.i
            }

            my.fit <- prems(
                y=y[train],
                x=x[train,,drop=FALSE],
                x.fixed=x.fixed.train,
                family=family,
                tau=tau,
                k.max=k.max,
                max.s=max.s,
                standardize=standardize,
                max2way=max2way,
                no.cores=no.cores,
                verbose=FALSE
            )

            if( use.saved.fits )
                saveRDS(my.fit, outfile)
        }

        ## Re-evaluate all requested model sizes. On continuation the
        ## existing model searches are not repeated.
        for( k in k.min:k.max ){

            kk <- k - k.min + 1

            pred <- predict.prems(
                my.fit,
                newx=x[test,,drop=FALSE],
                newx.fixed=x.fixed.test,
                size=k,
                criteria=criteria,
                fit=fit
            )

            if( family=='binomial' ){

                lp1 <- log(pred)
                lp1 <- ifelse(is.finite(lp1), lp1, -1000)

                lp0 <- log(1-pred)
                lp0 <- ifelse(is.finite(lp0), lp0, -1000)

                pwll[i,kk] <- sum(y[test]*lp1 + (1-y[test])*lp0)

            }else if( family=='gaussian' ){

                pwll[i,kk] <- sum(-(y[test] - pred)^2)

            }else if( family=='cox' ){

                pwll[i,kk] <- cox_pl_left_trunc(y[test], pred)
            }

            ## The last k coefficients are the selected non-fixed predictors:
            ## binomial/gaussian: intercept, fixed covariates, selected predictors
            ## cox:               fixed covariates, selected predictors
            selected.coef[[k]][i,] <- tail(
                names(getModelFit(my.fit, size=k, rank=1)$beta), k
            )
        }

        if( verbose )
            print(paste('Fold', i, 'complete.'))
    }

    cvm <- apply(pwll, 2, mean)
    cvsd <- apply(pwll, 2, sd)/sqrt(nfolds)

    names(cvm) <- k.min:k.max
    names(cvsd) <- k.min:k.max

    sizes <- (k.min:k.max)[prems.optim(cvm, cvsd)]

    ret <- list(
        best=sizes[1],
        one.se=sizes[2],
        cvm=cvm,
        cvsd=cvsd,
        selected.coef=selected.coef
    )

    return(ret)
}
