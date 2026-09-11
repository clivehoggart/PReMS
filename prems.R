# PReMS dependencies
library(parallel)
library(glmnet)
library(pROC)


cox_pl_left_trunc <- function( Surv_obj, eta ) {
    if( length(eta) != nrow(Surv_obj) ) stop("eta and surv must be same length")

    if( ncol(Surv_obj) == 2 ){
        start <- rep(-Inf, nrow(Surv_obj))
        stop <- Surv_obj[,1]
        status <- Surv_obj[,2]
    }else if( ncol(Surv_obj) == 3 ){
        start <- Surv_obj[,1]
        stop <- Surv_obj[,2]
        status <- Surv_obj[,3]
    }else{
        stop("Survival response must be Surv(time, status) or Surv(start, stop, status).")
    }

    ord <- order(stop)
    start <- start[ord]
    stop <- stop[ord]
    status <- status[ord]
    eta <- eta[ord]

    loglik <- 0
    for( i in which(status == 1) ){
        t_i <- stop[i]
        riskset <- which(start <= t_i & stop >= t_i)
        loglik <- loglik + eta[i] - log(sum(exp(eta[riskset])))
    }
    return(loglik)
}

plot.cv.prems <- function( cv.fit, ylim=NULL, cex=1 ){
    if( is.null(ylim) ){
        ylim <- range(c(cv.fit$cvm-cv.fit$cvsd, cv.fit$cvm+cv.fit$cvsd))
    }
    k <- as.numeric(names(cv.fit$cvm))
    plot(k, cv.fit$cvm, xlab='Model size', ylab='Cross-validated score',
         ylim=ylim, pch=19, cex.lab=cex, cex.axis=cex)
    arrows(k, cv.fit$cvm-cv.fit$cvsd, k, cv.fit$cvm+cv.fit$cvsd,
           angle=90, code=3, length=0.05)
    abline(v=k[cv.fit$one.se], lty=2)
    abline(v=k[cv.fit$best], lty=2)
}

getCoefGlmnet <- function( fit, s="lambda.min" ){
    beta <- as.matrix(stats::coef(fit, s=s))
    if( nrow(beta) > 0 && rownames(beta)[1] == "(Intercept)" ){
        beta <- beta[-1,,drop=FALSE]
    }
    beta <- beta[,1]
    beta <- beta[beta!=0]
    return(beta)
}

cv.auc <- function( y, pred, folds ){
    r <- vector()
    for( i in 1:max(folds) ){
        r[i] <- pROC::roc(y[folds==i], pred[folds==i], quiet=TRUE)$auc
    }
    return(mean(r))
}

getEta <- function( y, x, beta ){
  eta <- y * x %*% beta
  return(eta)
}

stepUP <- function( old.modelS, P, old.ml, max.s=10 ){
    s <- order( old.ml, decreasing=FALSE )
    new.modelS <- NULL
#  new.modelS <- expand.model( old.modelS[s[1],], P )
    max.s = min( max.s, length(s) )
    for( k in 1:max.s ){
        new.modelS <- rbind( new.modelS, expand.model( old.modelS[s[k],], P ) )
    }
    new.modelS <- unique(new.modelS)
    return(new.modelS)
}

getLogPost <- function( y, x, beta, tau ){
  eta <- getEta( y, x, beta )
  ptr <- eta > -709.5
  loglike <- -sum(log(1 + exp(-eta))[ptr]) + sum(eta[!ptr])
  logprior <- 0.5*sum(log(tau)) - 0.5*sum(tau*beta^2)
  # return minus log posterior for use by optim function which by default minimises
  return( -(loglike + logprior) )
}

getDLogPost <- function( y, x, beta, tau ){
  eta <- getEta( y, x, beta )
  dloglike <- t(as.matrix( y / (1 + exp(eta)) )) %*% x
  dlogprior <-  -tau*beta
  return( -(dloglike + dlogprior) )
}

getAll2Way <- function(n){
  models22 <- matrix(ncol=2,nrow=choose(n,2))
  ii <- 0
  for( I in 1:(n-1) ){
    for( J in (I+1):n ){
      ii <- ii+1
      models22[ii,] <- c(I,J)
    }
  }
  return(models22)
}

getHessian_gpt <- function( x, beta, tau ){# In matrix notation
    eta <- as.numeric(x %*% beta)
    mu <- 1 / ( 1 + exp(-eta) )

    omega <- pmax(mu * (1 - mu), 1e-12)
    H <- crossprod(x, x * omega) + diag(tau,nrow=length(tau))
    
    cholH <- chol(H)
    logdetH <- 2 * sum(log(diag(cholH)))

    return(logdetH)
}

expand.model <- function( old.model, P ){
  k <- length(old.model) + 1
  gamma <- setdiff( 1:P, old.model )
  new.model <- matrix( ncol=k, nrow=length(gamma) )
  for( i in 1:length(gamma) ){
    new.model[i,] <- sort( c(old.model,gamma[i]) )
  }
  return(new.model)
}

getMargLikelihood2 <- function( x.select=NULL, x.fixed=NULL, y, tau=1, family='gaussian', m1, sd1, m.fixed, sd.fixed, m.y, s.y ){
    n <- NROW(y)

    k1 <- ifelse( is.null(x.select), 0, ncol(x.select) ) # no. of selected covs
    k2 <- ifelse( is.null(x.fixed), 0, ncol(x.fixed) ) # no. fixed covs
    if( family=='gaussian' | family=='binomial' ){
        x1 <- cbind( rep(1,n), x.fixed, x.select )
        k <- ncol(x1) # total covs + intercept
        tau1 <- c( rep(1e-12,k2+1), rep(tau,k1) )
    }

    if( family=='cox' ){
        if( k1==0 & k2==0 ){
            fit <- survival::coxph( y ~ 1, model = FALSE, x = FALSE, y = FALSE )
        }else if( k1==0 & k2!=0 ){
            fit <- survival::coxph( y ~ x.fixed, model = FALSE, x = FALSE, y = FALSE )
        }else if( k1!=0 & k2==0 ){
            fit <- survival::coxph( y ~ survival::ridge( x.select, theta = tau, scale=FALSE ),
                         model = FALSE, x = FALSE, y = FALSE )
        }else{
            fit <- survival::coxph( y ~ x.fixed + survival::ridge( x.select, theta = tau, scale=FALSE ),
                         model = FALSE, x = FALSE, y = FALSE )
        }
        beta.tilde <- fit$coef
        l.gamma1 <- tail(fit$loglik, 1)
        aic <- -2*l.gamma1 + 2*length(beta.tilde)
    }
    if( family=='gaussian' ){
        penalty <- diag(c(rep(0, k2+1), rep(tau, k1)), k)
        A <- penalty + crossprod( x1, x1 )
        b = crossprod( x1, y )

        cholA <- chol(A)
        logdetA <- 2 * sum(log(diag(cholA)))
        Ainv_b <- backsolve(cholA, forwardsolve(t(cholA), b))
        S.tau <- as.numeric((n - 1) - crossprod(b, Ainv_b))
        q <- k2 + 1
        a.n <- 1 + (n - q) / 2
        l.gamma1 <- -a.n * log(1 + S.tau / 2) - 0.5 * logdetA
        
        beta.tilde <- as.numeric(Ainv_b)
        rss <- sum((y - as.numeric(x1 %*% beta.tilde))^2)
        df <- length(beta.tilde) - 1
        aic <- n * log(rss / n) + 2 * df
        beta.tilde <- beta.tilde * s.y
        beta.tilde[1] = beta.tilde[1] + m.y
    }
    if( family=='binomial' ){
        yy <- 2*y-1

        tmp <- optim( rep(0,k), fn=getLogPost, gr=getDLogPost, y=yy, x=x1, tau=tau1, method="L-BFGS" )
        beta.tilde <- tmp$par
        # Log-posterior is NEGATIVE of value which is returned by optim -- by default optim minimises
        logPost <- -tmp$value
#            hess <- getHessian2( x1, beta.tilde, tau1 )
#            l.gamma1 <- logPost - 0.5*log(det(hess))
        logdetH <- getHessian_gpt( x1, beta.tilde, tau1 )
        l.gamma1 <- logPost - 0.5 * logdetH
        eta <- getEta( yy, x1, beta.tilde )
        aic <- 2 * (sum(log(1 + exp(-eta))) + length(beta.tilde) - 1)
    }
    if( family=='binomial' | family=='gaussian' ){
        if( k2>0 ){
            for( ii in 1:k2 ){
                i <- ii + 1
                beta.tilde[i] <- beta.tilde[i] / sd.fixed[ii]
                beta.tilde[1] <- beta.tilde[1] - m.fixed[ii]*beta.tilde[i]
            }
        }
        if( k1>0 ){
            for( ii in 1:k1 ){
                i <- ii + k2 + 1
                beta.tilde[i] <- beta.tilde[i] / sd1[ii]
                beta.tilde[1] <- beta.tilde[1] - m1[ii]*beta.tilde[i]
            }
        }
    }else if( family=='cox' ){
        if( k2>0 ){
            for( i in 1:k2 ){
                beta.tilde[i] <- beta.tilde[i] / sd.fixed[i]
            }
        }
        if( k1>0 ){
            for( ii in 1:k1 ){
                i <- ii + k2
                beta.tilde[i] <- beta.tilde[i] / sd1[ii]
            }
        }
    }
    ret <- list( -l.gamma1, aic, beta.tilde )
######################################################################
# Returning MINUS log-posterior to be consistent with other measures #
# of model fit, ie ICs, which are minimised for best fit             #
######################################################################
    names(ret) <- c( 'ML', 'aic', 'beta' )
    return( ret )
}

############# Public functions below #############

prems <- function( y, x, x.fixed=NULL, max2way="all", k.max=5,
                  family='gaussian', tau=1, max.s=10, no.cores=10,
                  standardize=TRUE, verbose=TRUE ){
    model.indicator <- list()
    fitted.models <- list()

    m1 <- apply( x, 2, mean )
    s1 <- apply( x, 2, sd )
    if( !standardize ){
        s1 <- rep(1,ncol(x))
    }
    x <- t(t(x)-m1)
    x <- t(t(x)/s1)

    if( !is.null(x.fixed) ){
        m.fixed <- apply( x.fixed, 2, mean )
        s.fixed <- apply( x.fixed, 2, sd )
        if( !standardize ){
            s.fixed <- rep(1,ncol(x.fixed))
        }
        x.fixed <- t(t(x.fixed)-m.fixed)
        x.fixed <- t(t(x.fixed)/s.fixed)
    }else{
        m.fixed <- vector(length=0)
        s.fixed <- vector(length=0)
    }
    if( family=='gaussian' ){
        m.y <- mean(y)
        s.y <- sd(y)
        y <- ( y - m.y ) / s.y
    }else{
        m.y <- NULL
        s.y <- NULL
    }

    ptr.covs.use <- which( s1!=0 )
    if( length(ptr.covs.use) != ncol(x) & verbose ){
        print( paste('WARNING: Variables', setdiff( 1:ncol(x), ptr.covs.use ), 'are monomorphic.') )
    }
    Ncov <- length(ptr.covs.use)

    null <- getMargLikelihood2( y=y, x.fixed=x.fixed, family=family, tau=tau,
                               m1=vector(length=0), sd1=vector(length=0),
                               m.fixed=m.fixed, sd.fixed=s.fixed,
                               m.y=m.y, s.y=s.y )

    model.indicator[[1]] <- cbind(1:Ncov)
    if( verbose ){
        print( paste('Searching',Ncov,'1D models (all possible)') )
    }

    fitted.models[[1]] <- parallel::mclapply(1:Ncov, function(ptr)
    {getMargLikelihood2( x.select=x[,ptr.covs.use[ptr],drop=FALSE], x.fixed=x.fixed, y=y,
                        family=family, tau=tau,
                        m1=m1[ptr.covs.use[ptr]], sd1=s1[ptr.covs.use[ptr]],
                        m.fixed=m.fixed, sd.fixed=s.fixed,
                        m.y=m.y, s.y=s.y )},
    mc.cores=no.cores)
    if( verbose ){
        print("Finished 1D models")
    }

    if( max2way=='all' ){
        model.indicator[[2]] <- getAll2Way(Ncov)
        if( verbose ){
            print( paste('Searching',nrow(model.indicator[[2]]),'2D models (all possible)') )
        }
    }
    if( max2way!='all' ){
        ML <- unlist(parallel::mclapply( fitted.models[[1]], getElement, 'ML', mc.cores=no.cores ))
        model.indicator[[2]] <- stepUP( model.indicator[[1]], Ncov, ML, max.s=max2way )
        if( verbose ){
            print( paste('Searching',nrow(model.indicator[[2]]),'2D models') )
        }
    }
    k <- 2
    fitted.models[[k]] <- parallel::mclapply( 1:nrow(model.indicator[[k]]), function(i)
    {getMargLikelihood2( x.select=x[,ptr.covs.use[model.indicator[[k]][i,]]], x.fixed=x.fixed, y=y,
                        family=family, tau=tau,
                        m1=m1[ptr.covs.use[model.indicator[[k]][i,]]],
                        sd1=s1[ptr.covs.use[model.indicator[[k]][i,]]],
                        m.fixed=m.fixed, sd.fixed=s.fixed,
                        m.y=m.y, s.y=s.y )}, mc.cores=no.cores)
    if( verbose ){
        print("Finished 2D models")
    }
    ML <- unlist(parallel::mclapply( fitted.models[[2]], getElement, 'ML', mc.cores=no.cores ))

    s <- order( ML, decreasing=FALSE )
    n.keep <- min(max.s, length(s))
    tmp.fits <- vector("list", n.keep)
    tmp.indicator <- matrix(ncol=k, nrow=n.keep)
    iML <- numeric(n.keep)
    for( j in seq_len(n.keep) ){
        iML[j] <- ML[s[j]]
        tmp.fits[[j]] <- fitted.models[[2]][[s[j]]]
        tmp.indicator[j,] <- model.indicator[[2]][s[j],]
    }
    ML <- iML
    fitted.models[[2]] <- tmp.fits
    model.indicator[[2]] <- tmp.indicator

    if( k.max>2 ){
        for( k in 3:k.max ){
            model.indicator[[k]] <- stepUP( model.indicator[[(k-1)]], Ncov, ML, max.s=max.s )
            if( verbose ){
                print( paste('Searching ',nrow(model.indicator[[k]]),' ',k,'D models',sep='') )
            }
            fitted.models[[k]] <- parallel::mclapply(1:nrow(model.indicator[[k]]), function(i)
            {getMargLikelihood2( x.select=x[,ptr.covs.use[model.indicator[[k]][i,]]],
                                x.fixed=x.fixed, y=y,
                                family=family, tau=tau,
                                m1=m1[ptr.covs.use[model.indicator[[k]][i,]]],
                                sd1=s1[ptr.covs.use[model.indicator[[k]][i,]]],
                                m.fixed=m.fixed, sd.fixed=s.fixed,
                                m.y=m.y, s.y=s.y )},
            mc.cores=no.cores)
            if( verbose ){
                print( paste('Finished ',k,'D models',sep='') )
            }
            ML <- unlist(parallel::mclapply( fitted.models[[k]], getElement, 'ML', mc.cores=no.cores ))

            s <- order( ML, decreasing=FALSE )
            n.keep <- min(max.s, length(s))
            tmp.fits <- vector("list", n.keep)
            tmp.indicator <- matrix(ncol=k, nrow=n.keep)
            iML <- numeric(n.keep)
            for( j in seq_len(n.keep) ){
                iML[j] <- ML[s[j]]
                tmp.fits[[j]] <- fitted.models[[k]][[s[j]]]
                tmp.indicator[j,] <- model.indicator[[k]][s[j],]
            }
            ML <- iML
            fitted.models[[k]] <- tmp.fits
            model.indicator[[k]] <- tmp.indicator
        }
    }

    ret <- list( null, fitted.models, model.indicator, colnames(x.fixed), colnames(x),
                m1, s1, m.fixed, s.fixed, tau, standardize, family )
    names(ret) <- c('null','fitted.models','model.indicator', 'cnames.fixed', 'cnames',
                    'm', 'sd', 'm.fixed', 'sd.fixed', 'tau', 'standardize', 'family' )

    return( ret )
}

ModelSearchIncrease <- function( fitted.models, y, x, x.fixed=NULL, no.cores=10, max.s=NULL ){
    x <- t(t(x)-fitted.models$m)
    x <- t(t(x)/fitted.models$sd)

    if( !is.null(x.fixed) ){
        x.fixed <- t(t(x.fixed)-fitted.models$m.fixed)
        x.fixed <- t(t(x.fixed)/fitted.models$sd.fixed)
    }

    if( fitted.models$family=='gaussian' ){
        m.y <- mean(y)
        s.y <- sd(y)
        y.fit <- ( y - m.y ) / s.y
    }else{
        m.y <- NULL
        s.y <- NULL
        y.fit <- y
    }

    k <- length(fitted.models$fitted.models)
    ptr.covs.use <- which( fitted.models$sd!=0 )
    Ncov <- length(ptr.covs.use)

    # By default preserve the search breadth retained at the current largest
    # model size. max.s can be supplied explicitly to widen or narrow expansion.
    if( is.null(max.s) ){
        max.s <- length(fitted.models$fitted.models[[k]])
    }

    ML <- unlist(parallel::mclapply(
        fitted.models$fitted.models[[k]], getElement, 'ML', mc.cores=no.cores
    ))

    fitted.models$model.indicator[[k+1]] <- stepUP(
        fitted.models$model.indicator[[k]], Ncov, ML, max.s=max.s
    )

    fitted.models$fitted.models[[k+1]] <- parallel::mclapply(
        seq_len(nrow(fitted.models$model.indicator[[k+1]])),
        function(i){
            getMargLikelihood2(
                x.select=x[,ptr.covs.use[fitted.models$model.indicator[[k+1]][i,]],drop=FALSE],
                x.fixed=x.fixed,
                y=y.fit, family=fitted.models$family,
                tau=fitted.models$tau,
                m1=fitted.models$m[ptr.covs.use[fitted.models$model.indicator[[k+1]][i,]]],
                sd1=fitted.models$sd[ptr.covs.use[fitted.models$model.indicator[[k+1]][i,]]],
                m.fixed=fitted.models$m.fixed, sd.fixed=fitted.models$sd.fixed,
                m.y=m.y, s.y=s.y
            )
        },
        mc.cores=no.cores
    )

    # Keep the highest-ranking models, matching the behaviour of prems().
    ML.new <- unlist(parallel::mclapply(
        fitted.models$fitted.models[[k+1]], getElement, 'ML', mc.cores=no.cores
    ))
    s <- order(ML.new, decreasing=FALSE)
    n.keep <- min(max.s, length(s))
    keep <- s[seq_len(n.keep)]

    fitted.models$fitted.models[[k+1]] <- fitted.models$fitted.models[[k+1]][keep]
    fitted.models$model.indicator[[k+1]] <- fitted.models$model.indicator[[k+1]][keep,,drop=FALSE]

    return(fitted.models)
}

getModelFit <- function( fitted.models, size=1, rank=1, no.cores=10, criteria='ML' ){
    ptr.covs.use <- which( fitted.models$sd!=0 )

    ptr <- order( unlist(parallel::mclapply( fitted.models$fitted.models[[size]], getElement, criteria, mc.cores=no.cores ) ))[rank]

    model.fit <- list()
    ptr1 <- fitted.models$model.indicator[[size]][ptr,,drop=FALSE]
    for( i in 1:length(ptr) ){
        model.fit[[i]] <- fitted.models$fitted.models[[size]][[ptr[i]]]
        if( fitted.models$family=='binomial' | fitted.models$family=='gaussian' ){
            nmes <- c( 'I', fitted.models$cnames.fixed, fitted.models$cnames[ptr.covs.use[ptr1[i,]]] )
        }else if( fitted.models$family=='cox' ){
            nmes <- c( fitted.models$cnames.fixed, fitted.models$cnames[ptr.covs.use[ptr1[i,]]] )
        }
        names(model.fit[[i]]$beta) <- nmes
    }

    if( length(ptr)==1 ){
        model.fit <- model.fit[[1]]
    }
    return(model.fit)
}

thin.prems <- function( fit, size, rank ){
    fitted.models <- list()
    model.indicator <- list()
    for( i in seq_len(size) ){
        s <- order(unlist(sapply(fit$fitted.models[[i]], getElement, 'ML')))
        keep <- s[seq_len(min(rank, length(s)))]
        fitted.models[[i]] <- fit$fitted.models[[i]][keep]
        model.indicator[[i]] <- fit$model.indicator[[i]][keep,,drop=FALSE]
    }
    ret <- list(fit$null, fitted.models, model.indicator, fit$cnames.fixed, fit$cnames,
                fit$m, fit$sd, fit$m.fixed, fit$sd.fixed, fit$tau, fit$standardize, fit$family)
    names(ret) <- c('null','fitted.models','model.indicator', 'cnames.fixed', 'cnames',
                    'm', 'sd', 'm.fixed', 'sd.fixed', 'tau', 'standardize', 'family')
    return(ret)
}

predict.prems <- function( fitted.models, newx, newx.fixed=NULL, size=1, rank=1,
                          no.cores=10, criteria='ML' ){
    ptr.covs.use <- which( fitted.models$sd!=0 )
    best.fit <- order( unlist(parallel::mclapply( fitted.models$fitted.models[[size]], getElement, criteria, mc.cores=no.cores ) ))[rank]
    best.fit.model <- fitted.models$fitted.models[[size]][[best.fit]]$beta
    ptr <- fitted.models$model.indicator[[size]][best.fit,]

    ptr1 <- match( fitted.models$cnames[ptr.covs.use[ptr]], colnames(newx) )

    I = 1
    if( fitted.models$family=='cox' )
        I = NULL
    if( is.null(newx.fixed) ){
        X <- as.matrix(cbind( I, newx[,ptr1,drop=FALSE]) )
        pred <- X %*% best.fit.model
    }else{
        X <- as.matrix(cbind( I, newx.fixed, newx[,ptr1,drop=FALSE]) )
        pred <- X %*% best.fit.model
    }

    if( fitted.models$family=='binomial' ){
        pred <- 1 / ( 1 + exp(-pred) )
    }

    return(pred)
}

getICs <- function( fitted.models, k.min=1 ){
    ll <- length(fitted.models$model.indicator)
    res <- matrix( ncol=3, nrow=ll+2-k.min )
    for( ii in k.min:ll ){
        k <- ncol(fitted.models$model.indicator[[ii]])
        if( !is.null(k) ){
            aic <- min(sapply( fitted.models$fitted.models[[ii]], getElement, 'aic' ),na.rm=TRUE)
            ml <- min(sapply( fitted.models$fitted.models[[ii]], getElement, 'ML' ),na.rm=TRUE)
            res[(ii+2-k.min),] <- c( k, aic, ml )
        }
    }
    res[1,] <- c( 0, fitted.models$null$aic, fitted.models$null$ML )
    colnames(res) <- c('k', 'aic', 'ml' )
    return(res)
}

cv.prems <- function( y, x, x.fixed=NULL, no.cores=10, k.min=1, k.max, tau.i=NULL,
                      max.s=50, max2way='all', standardize=TRUE, nfolds=NULL, foldid=NULL,
                      lasso.factor=1, criteria='ML', family='binomial', verbose=TRUE ){
    n <- NROW(y)

    if( is.null(foldid) & is.null(nfolds) ){
        nfolds <- n
        foldid <- seq_len(nfolds)
    }
    if( !is.null(foldid) & is.null(nfolds) ){
        nfolds <- length(unique(foldid))
    }
    if( is.null(foldid) ){
        yy <- y
        if( family=="cox" )
            yy <- y[,ncol(y)]
        if( family=="cox" | family=="binomial" ){
            foldid <- make.folds2(yy, folds=nfolds)
            if( verbose ) print(table(yy, foldid))
        }
        if( family=="gaussian" ){
            foldid <- make.folds.continuous(n, nfolds)
            if( verbose ) print(table(foldid))
        }
    }

    pwll <- matrix(nrow=nfolds, ncol=(k.max-k.min+1))
    if( verbose ){
        print(paste(nfolds,'fold cross-validation'))
    }
    selected.coef <- vector("list", k.max)
    for( k in k.min:k.max ){
        selected.coef[[k]] <- matrix(nrow=nfolds, ncol=k)
    }

    subset_y <- function(y, idx) {
        if( family == "cox" ){
            y[idx,,drop=FALSE]
        }else{
            y[idx]
        }
    }

    for( i in seq_len(nfolds) ){
        train <- which(foldid!=i)
        test <- which(foldid==i)
        y.train <- subset_y(y, train)
        y.test <- subset_y(y, test)
        x.fixed.train <- if( is.null(x.fixed) ) NULL else x.fixed[train,,drop=FALSE]
        x.fixed.test <- if( is.null(x.fixed) ) NULL else x.fixed[test,,drop=FALSE]

        if( is.null(tau.i) ){
            tauest <- NULL
            attempt <- 1
            max.attempts <- 20
            tau.nfolds <- min(10, NROW(y.train))

            while( attempt <= max.attempts ){
                tauest.try <- try(
                    TauEst(y=y.train, x=x[train,,drop=FALSE], x.fixed=x.fixed.train,
                           family=family, nfolds=tau.nfolds, parallel=FALSE),
                    silent=TRUE
                )
                ok <- !inherits(tauest.try, "try-error") &&
                    !is.null(tauest.try$tau.opt) &&
                    length(tauest.try$tau.opt) == 1 &&
                    is.numeric(tauest.try$tau.opt) &&
                    is.finite(tauest.try$tau.opt) &&
                    tauest.try$tau.opt > 0
                if( ok ){
                    tauest <- tauest.try
                    break
                }
                attempt <- attempt + 1
            }
            if( is.null(tauest) ){
                stop("Unable to estimate tau in cross-validation fold ", i,
                     ". Supply tau.i explicitly or inspect the glmnet fit.")
            }
            tau <- tauest$tau.opt * lasso.factor
            if( verbose ) print(paste0("tau=",tau))
        }else{
            tau <- tau.i
        }

        my.fit <- prems(y=y.train, x=x[train,,drop=FALSE], x.fixed=x.fixed.train,
                        family=family, tau=tau, k.max=k.max, max.s=max.s,
                        standardize=standardize, max2way=max2way,
                        no.cores=no.cores, verbose=FALSE)

        for( k in k.min:k.max ){
            kk <- k - k.min + 1
            pred <- predict.prems(my.fit,
                                  newx=x[test,,drop=FALSE],
                                  newx.fixed=x.fixed.test,
                                  size=k, criteria=criteria)

            if( family=='binomial' ){
                lp1 <- log(pred)
                lp1 <- ifelse(is.finite(lp1), lp1, -1000)
                lp0 <- log(1-pred)
                lp0 <- ifelse(is.finite(lp0), lp0, -1000)
                pwll[i,kk] <- sum(y.test*lp1 + (1-y.test)*lp0)
            }else if( family=='gaussian' ){
                pwll[i,kk] <- sum(-(y.test-pred)^2)
            }else if( family=='cox' ){
                pwll[i,kk] <- cox_pl_left_trunc(y.test, pred)
            }

            coef.names <- names(getModelFit(my.fit, size=k, rank=1, criteria=criteria)$beta)
            selected <- setdiff(coef.names, c("I", my.fit$cnames.fixed))
            selected.coef[[k]][i,] <- selected
        }

        if( verbose ){
            print(paste('Fold',i,'complete.'))
        }
    }

    cvm <- apply(pwll, 2, mean)
    cvsd <- apply(pwll, 2, sd)/sqrt(nfolds)
    names(cvm) <- k.min:k.max
    names(cvsd) <- k.min:k.max

    sizes <- (k.min:k.max)[prems.optim(cvm, cvsd)]

    ret <- list(sizes[1], sizes[2], cvm, cvsd, selected.coef)
    names(ret) <- c('best','one.se','cvm','cvsd', 'selected.coef')
    return(ret)
}

prems.optim <- function( cvm, cvsd ){
    best <- order(cvm, decreasing=TRUE)[1]
    if( best == 1 ){
        one.se <- 1
    }else{
        ptr2 <- which((cvm+cvsd[best])[seq_len(best-1)] > cvm[best])
        if( length(ptr2)>0 ){
            one.se <- min(ptr2)
        }else{
            one.se <- best
        }
    }
    return(c(best,one.se))
}

TauEst <- function( y, x, x.fixed=NULL, family='binomial', standardize=TRUE,
                    n.coef=1, fit=NULL, nfolds=NULL, parallel=FALSE ){
    n <- NROW(y)
    if( is.null(nfolds) ){
        nfolds <- n
    }
    if( is.null(x.fixed) ){
        x.fixed <- matrix(nrow=n, ncol=0)
    }
    if( is.null(fit) ){
        lambda.factor <- c(rep(0,ncol(x.fixed)), rep(1,ncol(x)))
        fit <- glmnet::cv.glmnet(x=as.matrix(cbind(x.fixed,x)), y=y,
                      penalty.factor=lambda.factor,
                      family=family, alpha=1, nfolds=nfolds,
                      type.measure='deviance', grouped=FALSE, standardize=standardize,
                      parallel=parallel)
    }

    ncol.fixed <- ncol(x.fixed)

    remove.fixed <- function(beta) {
        if( ncol.fixed == 0 ) return(beta)
        fixed.names <- colnames(x.fixed)
        if( !is.null(fixed.names) ){
            return(beta[!names(beta) %in% fixed.names])
        }
        beta[-seq_len(min(ncol.fixed, length(beta)))]
    }

    lambda.min <- fit$lambda.min
    beta <- remove.fixed(getCoefGlmnet(fit, s=lambda.min))
    if( length(beta) == 0 && length(fit$lambda) >= 2 ){
        lambda.min <- fit$lambda[2]
        beta <- remove.fixed(getCoefGlmnet(fit, s=lambda.min))
    }
    if( length(beta) == 0 ){
        stop("No penalised predictors were selected by Lasso; tau cannot be estimated.")
    }

    s <- rep(1,length(beta))
    if( standardize ){
        ptr <- match(names(beta), colnames(x))
        if( anyNA(ptr) ){
            stop("Predictor names in the glmnet fit do not match colnames(x).")
        }
        s <- apply(x[,ptr,drop=FALSE], 2, sd)
    }
    lambda <- lambda.min * n
    beta1 <- sort(abs(beta*s), decreasing=TRUE)
    n.use <- min(n.coef, length(beta1))
    ptr <- seq_len(n.use)
    tau.opt <- lambda * sum(beta1[ptr]) / sum(beta1[ptr]^2)

    beta <- remove.fixed(getCoefGlmnet(fit, s='lambda.1se'))
    if( length(beta) > 0 ){
        s <- rep(1,length(beta))
        if( standardize ){
            ptr <- match(names(beta), colnames(x))
            if( anyNA(ptr) ){
                stop("Predictor names in the glmnet fit do not match colnames(x).")
            }
            s <- apply(x[,ptr,drop=FALSE], 2, sd)
        }
        lambda <- fit$lambda.1se * n
        beta1 <- sort(abs(beta*s), decreasing=TRUE)
        n.use <- min(n.coef, length(beta1))
        ptr <- seq_len(n.use)
        tau.1se <- lambda * sum(beta1[ptr]) / sum(beta1[ptr]^2)
    }else{
        tau.1se <- NA
    }

    ret <- list(tau.opt, tau.1se, fit)
    names(ret) <- c('tau.opt', 'tau.1se', 'fit.lasso')
    return(ret)
}

my.auc <- function( my.fit, sizes, X, y, rank=1, criteria='ML' ){
    my.pred <- matrix(ncol=length(sizes),nrow=length(y))
    for( i in sizes ){
        ii <- i - min(sizes) + 1
        my.pred[,ii] <- predict.prems( my.fit, as.matrix(X), size=i, rank=rank, no.cores=10, criteria=criteria )
    }
#    r <- matrix(ncol=3,nrow=length(sizes))
    r <- list()
    for( i in 1:length(sizes) ){
        r[[i]] <- pROC::roc(y, my.pred[,i], ci=TRUE, quiet=TRUE)
#        r[i,] <- as.numeric(pROC::roc(y, my.pred[,i], ci=TRUE, quiet=TRUE)$ci)
    }
    names(r) <- sizes
    return(r)
}

make.folds2 <- function(strata, folds=5, seed=NULL) {

    if( !is.null(seed) ) set.seed(seed)

    if( length(folds)!=1 || folds<2 || folds!=as.integer(folds) ){
        stop("folds must be a single integer >= 2")
    }

    if( any(is.na(strata)) ){
        stop("strata contains missing values")
    }

    strata <- as.factor(strata)
    n <- length(strata)

    if( folds > n ){
        stop("number of folds cannot exceed number of samples")
    }

    foldid <- rep(NA_integer_, n)

    tab <- table(strata)
    if( any(tab < folds) ){
        warning("At least one stratum has fewer samples than folds; exact stratification is impossible")
    }

    ## Process smaller strata first, so rare strata are spread as evenly as possible
    strata.levels <- names(sort(tab))

    fold.counts <- rep(0, folds)

    for( s in strata.levels ){

        ptr <- which(strata==s)
        ptr <- sample(ptr)

        ## Prefer currently smaller folds, but randomise ties
        fold.order <- order(fold.counts, runif(folds))

        ## Repeat fold order until all samples in this stratum are assigned
        assign.folds <- rep(fold.order, length.out=length(ptr))

        foldid[ptr] <- assign.folds

        fold.counts <- tabulate(foldid[!is.na(foldid)], nbins=folds)
    }

    return(foldid)
}

make.folds.continuous <- function( n, folds ) {
  if (folds <= 0) stop("folds must be positive.")
  if (n < 0) stop("n must be nonnegative.")

  q <- n %/% folds
  r <- n %% folds

  counts <- rep(q, folds)
  if (r > 0) {
    counts[1:r] <- counts[1:r] + 1
  }

  x <- rep(seq_len(folds), times = counts)
  sample(x)
}
