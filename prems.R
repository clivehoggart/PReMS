library(parallel)
library(glmnet)
library(pROC)
#library(BayesLogit)
library(gplots)
library(MASS)
library(exvatools)

cox_pl_left_trunc <- function( Surv_obj, eta ) {
    if( length(eta) != nrow(Surv_obj) ) stop("eta and surv must be same length")

    start  <- Surv_obj[, 1]
    stop   <- Surv_obj[, 2]
    status <- Surv_obj[, 3]

  # order by stop time
    ord <- order(stop)
    start  <- start[ord]
    stop   <- stop[ord]
    status <- status[ord]
    eta    <- eta[ord]

    loglik <- 0
    for( i in which(status == 1) ){
        t_i <- stop[i]
        riskset <- which( start <= t_i & stop >= t_i )
        loglik <- loglik + eta[i] - log(sum(exp(eta[riskset])))
    }
    return(loglik)
}

my.rpg <- function(z){
    omega <- rpg( num=1, h=1, z )
    return(omega)
}

Bayes_logit <- function( y=y, X, m0, P0, samp, burn ){
    beta.samp <- matrix( ncol=length(m0), nrow=(samp) )
    beta <- rep(0,length(m0))
    kappa <- y-0.5
    omega <- vector(length=length(y))
    tX <- t(X)
    P0.m0 <- P0 %*% m0
    tX.kappa <- tX %*% kappa
    m1 <- tX.kappa + P0.m0
    for( k in 1:(samp+burn) ){
        eta <- X %*% beta
        omega <- sapply( as.vector(eta), my.rpg )
#        V.omega <- solve(t(X) %*% diag(omega) %*% X + P0)
        V.omega <- solve( multd(tX,diag(omega)) %*% X + P0)
        m.omega <- V.omega %*% m1
        beta <- mvrnorm( n=1, mu=m.omega, Sigma=V.omega )
        if( k>burn ){
            beta.samp[(k-burn),] <- beta
        }
    }
    return(beta.samp)
}

plot.cv.prems <- function( cv.fit, ylim=NULL, cex=1 ){
    if( is.null(ylim) ){
        ylim <- range(c( (cv.fit$cvm-cv.fit$cvsd), (cv.fit$cvm+cv.fit$cvsd) ))
    }
    k <- as.numeric(names(cv.fit$cvm))
    plotCI( k, cv.fit$cvm, col='red', xlab='Model size', ylab='', uiw=cv.fit$cvsd, barcol='dimgrey', ylim=ylim,pch=19, cex.lab=cex, cex.axis=cex )
    mtext('Predictive log-likelihood',side=2,line=3,cex=cex)
    abline(v=(k)[cv.fit$one.se], lty=2)
    abline(v=(k)[cv.fit$best], lty=2)
}

getCoefGlmnet <- function( fit, s="lambda.min" ){
    tmp2 <- coef( fit, s=s )
    ptr <- which( tmp2!=0 )
    tmp <- tmp2[ptr]
    names(tmp) <- rownames(tmp2)[ptr]
    beta <- tmp[-1]
    return(beta)
}

cv.auc <- function( y, pred, folds ){
    r <- vector()
    for( i in 1:max(folds) ){
        r[i] <- roc( y[folds==i], pred[folds==i] )$auc
    }
    return(mean(r))
}

null.sim <- function( fitted.models, X, y, from=from, samples, no.cores=10 ){
    waic.null <- vector()
    for( i in 1:samples ){
        selected.covs <- names(getModelFit( fitted.models, from )$beta)[-1]
        ptr <- match( selected.covs, colnames(X) )
        x.null <- X[sample(1:nrow(X)),]
        x.null[,ptr] <- X[,ptr]

        tmp <- ModelSearchIncrease( fitted.models=fitted.models,  X=x.null, y=y, no.cores=no.cores, from=from )
        tmp <- fillICs( fitted.models=tmp, y=y, X=x.null, n.waic=5000, model.sizes=from+1, no.cores=no.cores )
        tmp2 <- getICs(tmp)
        waic.null[i] <- tmp2[ from+2, 4 ]
    }
    return(waic.null)
}

getEta <- function( y, x, beta ){
  eta <- y * x %*% beta
  return(eta)
}

stepUP <- function( old.modelS, P, old.ml, max.s=10 ){
  s <- order( old.ml, decreasing=FALSE )
  new.modelS <- expand.model( old.modelS[s[1],], P )
  for( k in 2:max.s ){
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

getHessian2 <- function( x, beta, tau ){# In matrix notation
  eta <- x %*% beta
  theta <- 1 / ( 1 + exp(-eta) )
  R <- diag(x=as.vector(theta*(1-theta)))
  hess <- t(x) %*% R %*% x + diag(tau,nrow=length(tau))
  return(hess)
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

getMargLikelihood.clean <- function( x.select=NULL, x.fixed=NULL, y, tau=1, delta=1, family='gaussian', beta.tilde0=NULL, m1, sd1, m.fixed, sd.fixed, n.waic=0, burn=10 ){
    n <- length(y)
    x1 <- cbind( rep(1,n), x.fixed, x.select )
    k <- ncol(x1) # total covs + intercept
    k1 <- ifelse( is.null(x.select), 0, ncol(x.select) ) # no. of selected covs
    k2 <- ifelse( is.null(x.fixed), 0, ncol(x.fixed) ) # no. fixed covs
    w.aic <- NA
    lppd <- NA
    p.waic <- NA
    p.waic2 <- NA
    p.anna <- NA
    dic <- NA
    beta.bar <- NA

    if( is.null(beta.tilde0) ){
        beta.tilde0 <- rep(0,k)
    }

    if( family=='cox' ){
        fit <- coxph( y ~ x.fixed + ridge( x.select, theta = tau, scale=FALSE ) )
        beta.tilde <- fit$coef
        l.gamma1 <- fit$loglik[2]
    }
    if( family=='gaussian' | family=='binomial' ){
        fit <- glm( y ~ x.fixed + ridge( x.select, theta = tau, scale=FALSE ), family=family )
        beta.tilde <- fit$coef
        l.gamma1 <- -fit$deviance / 2
        aic <- fit$aic
    }
    if( k2>0 ){
        for( ii in 1:k2 ){
            i <- ii + 1
            beta.tilde[i] <- beta.tilde[i] / sd.fixed[ii]
            beta.tilde[1] <- beta.tilde[1] - m.fixed[ii]*beta.tilde[i]
            beta.bar[i] <- beta.bar[i] / sd.fixed[ii]
            beta.bar[1] <- beta.bar[1] - m.fixed[ii]*beta.bar[i]
        }
    }
    if( k1>0 ){
        for( ii in 1:k1 ){
            i <- ii + k2 + 1
            beta.tilde[i] <- beta.tilde[i] / sd1[ii]
            beta.tilde[1] <- beta.tilde[1] - m1[ii]*beta.tilde[i]
            beta.bar[i] <- beta.bar[i] / sd1[ii]
            beta.bar[1] <- beta.bar[1] - m1[ii]*beta.bar[i]
        }
    }
    ret <- list( -l.gamma1, aic, w.aic, beta.tilde, beta.bar, lppd, p.waic )
######################################################################
# Returning MINUS log-posterior to be consistent with other measures #
# of model fit, ie ICs, which are minimised for best fit             #
######################################################################
    names(ret) <- c( 'ML', 'aic', 'waic', 'beta', 'beta.bar', 'lppd', 'p.waic' )
    return( ret )
}

getMargLikelihood2 <- function( x.select=NULL, x.fixed=NULL, y, tau=1, delta=1, family='gaussian', beta.tilde0=NULL, m1, sd1, m.fixed, sd.fixed, n.waic=0, burn=10 ){
    n <- length(y)
    x1 <- cbind( rep(1,n), x.fixed, x.select )
    k <- ncol(x1) # total covs + intercept
    k1 <- ifelse( is.null(x.select), 0, ncol(x.select) ) # no. of selected covs
    k2 <- ifelse( is.null(x.fixed), 0, ncol(x.fixed) ) # no. fixed covs
    w.aic <- NA
    lppd <- NA
    p.waic <- NA
    p.waic2 <- NA
    p.anna <- NA
    dic <- NA
    beta.bar <- NA

    if( is.null(beta.tilde0) ){
        beta.tilde0 <- rep(0,k)
    }

    if( family=='cox' ){
        if( k1==0 & k2==0 ){
            fit <- coxph( y ~ 1 )
        }else if( k1==0 & k2!=0 ){
            fit <- coxph( y ~ x.fixed )
        }else if( k1!=0 & k2==0 ){
            fit <- coxph( y ~ ridge( x.select, theta = tau, scale=FALSE ) )
        }else{
            fit <- coxph( y ~ x.fixed + ridge( x.select, theta = tau, scale=FALSE ) )
        }
        beta.tilde <- fit$coef
        l.gamma1 <- fit$loglik[2]
        aic = -2*fit$loglik[2] + 2*length(beta.tilde)
    }
    if( family=='gaussian' ){
        norm.const <- sd(y)
        y <- y/norm.const
        y.y <- sum(y*y)

        M1 <- diag(tau,nrow=ncol(x1)) + t(x1) %*% x1
        Inv.M1 <- solve(M1)
        q1 <- y.y - (t(y) %*% x1 %*% Inv.M1 %*% t(x1) %*% y)

        l.gamma1 <- lgamma( 0.5*(n+delta+k1) ) - (n-k1)*log(tau)/2 - (n+delta+k1)*log(1+q1/tau)/2 - lgamma(0.5*(delta+k1)) - 0.5*log(det(M1))
        beta.tilde <- (Inv.M1 %*% t(x1) %*% y) * norm.const
        aic <- 2* (sum(y - x1 %*% beta.tilde)^2 + length(beta.tilde) - 1)
        if(n.waic!=0){
            S <- t(y - x1 %*% beta.tilde) %*% y
            post.samples <- rmvt( delta=beta.tilde, sigma=S, type="shifted" )
        }
    }
    if( family=='binomial' ){
        tau1 <- c( rep(1e-12,k2+1), rep(tau,k1) )
        yy <- 2*y-1
        tmp <- optim( beta.tilde0, fn=getLogPost, gr=getDLogPost, y=yy, x=x1, tau=tau1, method="L-BFGS" )
        beta.tilde <- tmp$par

        # Log-posterior is NEGATIVE of value which is returned by optim -- by default optim minimises
        logPost <- -tmp$value
#        hess <- getHessian2( x1, beta.tilde, tau1 )
        ###########################
        # Recording log-posterior #
        ###########################
        l.gamma1 <- logPost# - 0.5*log(det(hess))

        eta <- getEta( yy, x1, beta.tilde )
        aic <- 2 * (sum(log(1 + exp(-eta))) + length(beta.tilde) - 1)
        if(n.waic!=0){
            beta.post <- Bayes_logit( y=y, X=x1, m0=rep(0,k), P0=diag(tau1,k), samp=n.waic, burn=burn )
#            beta.post <- post.samples$beta
            ll <- matrix( ncol=n, nrow=n.waic )
            for( i in 1:n.waic ){
                eta <- getEta( yy, x1, beta.post[i,] )
                ll[i,] <- 1/(1 + exp(-eta))
            }
            lppd <- sum(log(apply( ll, 2, mean )))
            p.waic <- sum(apply( log(ll), 2, var ))
#            p.waic2 <- 2 * sum( log(apply( ll, 2, mean )) - apply( log(ll), 2,  mean ) )
            w.aic <- -lppd + p.waic

            beta.bar <- apply( beta.post, 2, mean )
        }
    }
    if( family=='binomial' | family=='gaussian' ){
        if( k2>0 ){
            for( ii in 1:k2 ){
                i <- ii + 1
                beta.tilde[i] <- beta.tilde[i] / sd.fixed[ii]
                beta.tilde[1] <- beta.tilde[1] - m.fixed[ii]*beta.tilde[i]
                beta.bar[i] <- beta.bar[i] / sd.fixed[ii]
                beta.bar[1] <- beta.bar[1] - m.fixed[ii]*beta.bar[i]
            }
        }
        if( k1>0 ){
            for( ii in 1:k1 ){
                i <- ii + k2 + 1
                beta.tilde[i] <- beta.tilde[i] / sd1[ii]
                beta.tilde[1] <- beta.tilde[1] - m1[ii]*beta.tilde[i]
                beta.bar[i] <- beta.bar[i] / sd1[ii]
                beta.bar[1] <- beta.bar[1] - m1[ii]*beta.bar[i]
            }
        }
    }else if( family=='cox' ){
        if( k2>0 ){
            for( i in 1:k2 ){
                beta.tilde[i] <- beta.tilde[i] / sd.fixed[i]
                beta.bar[i] <- beta.bar[i] / sd.fixed[i]
            }
        }
        if( k1>0 ){
            for( ii in 1:k1 ){
                i <- ii + k2
                beta.tilde[i] <- beta.tilde[i] / sd1[ii]
                beta.bar[i] <- beta.bar[i] / sd1[ii]
            }
        }
    }
    ret <- list( -l.gamma1, aic, w.aic, beta.tilde, beta.bar, lppd, p.waic )
######################################################################
# Returning MINUS log-posterior to be consistent with other measures #
# of model fit, ie ICs, which are minimised for best fit             #
######################################################################
    names(ret) <- c( 'ML', 'aic', 'waic', 'beta', 'beta.bar', 'lppd', 'p.waic' )
    return( ret )
}

############# Public functions below #############

prems <- function( y, x, x.fixed=NULL, max2way="all", k.max=5, omega=0.5,
                  family='gaussian', tau=1, delta=1, max.s=10, no.cores=10,
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

    ptr.covs.use <- which( s1!=0 )
    if( length(ptr.covs.use) != ncol(x) & verbose ){
        print( paste('WARNING: Variables', setdiff( 1:ncol(x), ptr.covs.use ), 'are monomorphic.') )
    }
    Ncov <- length(ptr.covs.use)

    null <- getMargLikelihood2( y=y, x.fixed=x.fixed, family=family, tau=tau, delta=delta,
                               m1=vector(length=0), sd1=vector(length=0),
                               m.fixed=m.fixed, sd.fixed=s.fixed,
                               n.waic=0, burn=10 )

    model.indicator[[1]] <- cbind(1:Ncov)
    if( verbose ){
        print( paste('Searching',Ncov,'1D models (all possible)') )
    }

    fitted.models[[1]] <- mclapply(1:Ncov, function(ptr)
    {getMargLikelihood2( x.select=x[,ptr.covs.use[ptr],drop=FALSE], x.fixed=x.fixed,
                        y=y, family=family, tau=tau, delta=delta,
                        m1=m1[ptr.covs.use[ptr]], sd1=s1[ptr.covs.use[ptr]],
                        m.fixed=m.fixed, sd.fixed=s.fixed )},
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
        ML <- unlist(mclapply( fitted.models[[1]], getElement, 'aic', mc.cores=no.cores ))
        model.indicator[[2]] <- stepUP( model.indicator[[1]], Ncov, ML, max.s=max2way )
        if( verbose ){
            print( paste('Searching',nrow(model.indicator[[2]]),'2D models') )
        }
    }
    k <- 2
    fitted.models[[k]] <- mclapply( 1:nrow(model.indicator[[k]]), function(i)
    {getMargLikelihood2( x.select=x[,ptr.covs.use[model.indicator[[k]][i,]]],
                        x.fixed=x.fixed, y=y,
                        family=family, tau=tau, delta=delta,
                        m1=m1[ptr.covs.use[model.indicator[[k]][i,]]],
                        sd1=s1[ptr.covs.use[model.indicator[[k]][i,]]],
                        m.fixed=m.fixed, sd.fixed=s.fixed )}, mc.cores=no.cores)
    if( verbose ){
        print("Finished 2D models")
    }
    ML <- unlist(mclapply( fitted.models[[2]], getElement, 'aic', mc.cores=no.cores ))

    s <- order( ML, decreasing=FALSE )
    tmp.fits <- list()
    tmp.indicator <- matrix( ncol=k, nrow=max.s )
    iML <- vector()
    for( j in 1:max.s ){
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
            fitted.models[[k]] <- mclapply(1:nrow(model.indicator[[k]]), function(i)
            {getMargLikelihood2( x.select=x[,ptr.covs.use[model.indicator[[k]][i,]]],
                                x.fixed=x.fixed, y=y,
                                family=family, tau=tau, delta=delta,
                                m1=m1[ptr.covs.use[model.indicator[[k]][i,]]],
                                sd1=s1[ptr.covs.use[model.indicator[[k]][i,]]],
                                m.fixed=m.fixed, sd.fixed=s.fixed )},
            mc.cores=no.cores)
            if( verbose ){
                print( paste('Finished ',k,'D models',sep='') )
            }
            ML <- unlist(mclapply( fitted.models[[k]], getElement, 'aic', mc.cores=no.cores ))

            s <- order( ML, decreasing=FALSE )
            tmp.fits <- list()
            tmp.indicator <- matrix( ncol=k, nrow=max.s )
            iML <- vector()
            for( j in 1:max.s ){
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

ModelSearchIncrease <- function( fitted.models, y, x, x.fixed=NULL, no.cores=10, max.s=max.s ){
    x <- t(t(x)-fitted.models$m)
    x <- t(t(x)/fitted.models$sd)

    if( !is.null(x.fixed) ){
        x.fixed <- t(t(x.fixed)-fitted.models$m.fixed)
        x.fixed <- t(t(x.fixed)/fitted.models$sd.fixed)
    }

    k <- length(fitted.models$fitted.models)
    ptr.covs.use <- which( fitted.models$sd!=0 )
    Ncov <- length(ptr.covs.use)

    ML <- unlist(mclapply( fitted.models$fitted.models[[k]], getElement, 'aic', mc.cores=no.cores ))

    fitted.models$model.indicator[[k+1]] <- stepUP( fitted.models$model.indicator[[k]], Ncov, ML, max.s=max.s )

    fitted.models$fitted.models[[k+1]] <- mclapply(1:nrow(fitted.models$model.indicator[[k+1]]) , function(i)
    {getMargLikelihood2( x.select=x[,ptr.covs.use[fitted.models$model.indicator[[k+1]][i,]]],
                        x.fixed=x.fixed,
                        y=y, family=fitted.models$family,
                        tau=fitted.models$tau, delta=fitted.models$delta,
                        m1=fitted.models$m[ptr.covs.use[fitted.models$model.indicator[[k+1]][i,]]],
                        sd1=fitted.models$sd[ptr.covs.use[fitted.models$model.indicator[[k+1]][i,]]],
                        m.fixed=fitted.models$m.fixed, sd.fixed=fitted.models$sd.fixed )},
    mc.cores=no.cores)
    return(fitted.models)
}

fillICs <- function( fitted.models, y, x, x.fixed=NULL, n.waic, no.cores=10, n.rank=10, model.sizes, verbose=TRUE ){
    x <- t(t(x)-fitted.models$m)
    x <- t(t(x)/fitted.models$sd)

    if( !is.null(x.fixed) ){
        x.fixed <- t(t(x.fixed)-fitted.models$m.fixed)
        x.fixed <- t(t(x.fixed)/fitted.models$sd.fixed)
    }

    ptr.covs.use <- which( fitted.models$sd!=0 )
    for( size in model.sizes ){
        ic <- sapply( fitted.models$fitted.models[[size]], getElement, 'aic' )
        s <- order(ic)
        ptr <- fitted.models$model.indicator[[size]][s[1:n.rank],, drop=FALSE]
        tmp <-  mclapply( 1:n.rank,
                         function(ii){getMargLikelihood2(
                                          x.select=x[ ,ptr.covs.use[ptr[ii,]], drop=FALSE],
                                          x.fixed=x.fixed,
                                          y=y, tau=fitted.models$tau,
                                          family=fitted.models$family,
                                          n.waic=n.waic,
                                          m1=fitted.models$m[ptr.covs.use[ptr[ii,]]],
                                          sd1=fitted.models$sd[ptr.covs.use[ptr[ii,]]],
                                          m.fixed=fitted.models$m.fixed,
                                          sd.fixed=fitted.models$sd.fixed )},
                         mc.cores=no.cores )
        for( i in 1:n.rank ){
            fitted.models$fitted.models[[size]][[s[i]]] <- tmp[[i]]
        }
        if( verbose ){
            print(paste(size,'D models complete',sep=''))
        }
    }
    return( fitted.models )
}

getModelFit <- function( fitted.models, size=NULL, rank=1, no.cores=10, criteria='aic', family=family ){
    ptr.covs.use <- which( fitted.models$sd!=0 )
    model.sizes <- size
    if( is.null(model.sizes) ){
        model.sizes=1:length(fitted.models)
    }
    if( length(model.sizes)>1 ){
        min.ic <- vector()
        for( i in model.sizes ){
            min.ic[i] <- min(unlist(mclapply( fitted.models$fitted.models[[i]], getElement, criteria, mc.cores=no.cores )),na.rm=TRUE)
        }
        size <- order(min.ic)[1]
    }

    ptr <- order( unlist(mclapply( fitted.models$fitted.models[[size]], getElement, criteria, mc.cores=no.cores ) ))[rank]

    model.fit <- list()
    ptr1 <- fitted.models$model.indicator[[size]][ptr,,drop=FALSE]
    for( i in 1:length(ptr) ){
        model.fit[[i]] <- fitted.models$fitted.models[[size]][[ptr[i]]]
        if( family=='binomial' | family=='gaussian' ){
            nmes <- c( 'I', fitted.models$cnames.fixed, fitted.models$cnames[ptr.covs.use[ptr1[i,]]] )
        }else if( family=='cox' ){
            nmes <- c( fitted.models$cnames.fixed, fitted.models$cnames[ptr.covs.use[ptr1[i,]]] )
        }
        names(model.fit[[i]]$beta) <- nmes
        names(model.fit[[i]]$beta.bar) <- nmes
    }

    if( length(ptr)==1 ){
        model.fit <- model.fit[[1]]
    }
    return(model.fit)
}

predict.prems <- function( fitted.models, newx, newx.fixed=NULL, size=NULL, rank=1, no.cores=10, criteria='waic', fit='mean', family='binomial' ){
    ptr.covs.use <- which( fitted.models$sd!=0 )
    best.fit <- order( unlist(mclapply( fitted.models$fitted.models[[size]], getElement, criteria, mc.cores=no.cores ) ))[rank]
    if( fit=='mode' ){
        best.fit.model <- fitted.models$fitted.models[[size]][[best.fit]]$beta
    }
    if( fit=='mean' ){
        best.fit.model <- fitted.models$fitted.models[[size]][[best.fit]]$beta.bar
    }
    ptr <- fitted.models$model.indicator[[size]][best.fit,]

    ptr1 <- match( fitted.models$cnames[ptr.covs.use[ptr]], colnames(newx) )

    I = 1
    if( family=='cox' )
        I = NULL
    if( is.null(newx.fixed) ){
        X <- as.matrix(cbind( I, newx[,ptr1,drop=FALSE]) )
        pred <- X %*% best.fit.model
    }else{
        X <- as.matrix(cbind( I, newx.fixed, newx[,ptr1,drop=FALSE]) )
        pred <- X %*% best.fit.model
    }

    if( family=='binomial' ){
        pred <- 1 / ( 1 + exp(-pred) )
    }

    return(pred)
}

predict.prems.bayes <- function( fitted.models, newx, y.train, x.train, size=NULL, rank=1, iter=10000, no.cores=10, criteria='waic' ){
    ptr.covs.use <- which( fitted.models$sd!=0 )
    for( i in 1:ncol(x.train) ){
        newx[,i] <- (newx[,i]-fitted.models$m[i]) / fitted.models$sd[i]
        x.train[,i] <- (x.train[,i]-fitted.models$m[i]) / fitted.models$sd[i]
    }

    best.fit <- order( unlist(mclapply( fitted.models$fitted.models[[size]], getElement, criteria ) ))[rank]
    ptr1 <- fitted.models$model.indicator[[size]][best.fit,]
    ptr <- match( fitted.models$cnames[ptr.covs.use[ptr1]], colnames(x.train) )


    tau1 <- c( 1e-12, rep(fitted.models$tau,size) )
    iter.cores <- ceiling(iter/no.cores)
    beta <- mclapply(1:no.cores, function(i) {parallel.BayesLogit( y=y.train, X=cbind(1,x.train[,ptr]), m0=rep(0,size+1), P0=diag(tau1), samp=iter.cores, burn=500, dummy=i )}, mc.cores=no.cores )

#    beta <- lapply( post.samples, getElement, 'beta' )

    ptr.test <- match( colnames(x.train)[ptr], colnames(newx) )
    X <- cbind( 1, newx[,ptr.test,drop=FALSE] )
    pred.full <- matrix(ncol=no.cores,nrow=nrow(newx))
    beta2 <- matrix(ncol=ncol(X),nrow=0)
    for( i in 1:no.cores ){
        tmp <- mclapply(1:iter.cores, function(j){ 1/( 1+exp(-X %*% beta[[i]][j,]))}, mc.cores=no.cores )
        theta <- simplify2array(tmp,higher=FALSE)
        pred.full[,i] <- apply( theta, 1, mean )
        beta2 <- rbind( beta2, beta[[i]] )
    }
    pred.full <- apply( pred.full, 1, mean )
    beta.bar <- apply( beta2, 2, mean )
    eta.hat <- X %*% beta.bar
    pred.hat <- 1/(1 + exp(-eta.hat))

    return( cbind( pred.full, pred.hat) )
}

parallel.BayesLogit <- function( y, X, m0, P0, samp, burn, dummy ){
    post.samples <- Bayes_logit( y=y, X=X, m0=m0, P0=P0, samp=samp, burn=burn )
    return(post.samples)
}

getICs <- function( fitted.models, k.min=1 ){
    ll <- length(fitted.models$model.indicator)
    res <- matrix( ncol=4, nrow=ll+2-k.min )
    for( ii in k.min:ll ){
        k <- ncol(fitted.models$model.indicator[[ii]])
        if( !is.null(k) ){
            aic <- min(sapply( fitted.models$fitted.models[[ii]], getElement, 'aic' ),na.rm=TRUE)
            waic <- min(sapply( fitted.models$fitted.models[[ii]], getElement, 'waic' ),na.rm=TRUE)
            ml <- min(sapply( fitted.models$fitted.models[[ii]], getElement, 'ML' ),na.rm=TRUE)
            res[(ii+2-k.min),] <- c( k, aic, waic, ml )
        }
    }
    res[1,] <- c( 0, fitted.models$null$aic, fitted.models$null$waic, fitted.models$null$ML )
    colnames(res) <- c('k', 'aic', 'waic', 'ml' )
    return(res)
}

prems.clean <- function( all.fits ){
    model.indicator <- list()
    fitted.models <- list()
    for( i in 1:length(all.fits$fitted.models) ){
        ptr <- which(!is.na(sapply( all.fits$fitted.models[[i]], getElement, 'waic' )))
        if( length(ptr)>0 ){
            fitted.models[[i]] <- list()
            model.indicator[[i]] <- matrix(ncol=i,nrow=length(ptr))
            for( j in 1:length(ptr) ){
                fitted.models[[i]][[j]] <- all.fits$fitted.models[[i]][[ptr[j]]]
                tmp <- all.fits$model.indicator[[i]][ptr[j],]
                model.indicator[[i]][j,] <- tmp
            }
        }
    }
    ret <- list( all.fits$null, fitted.models, model.indicator, all.fits$cnames, all.fits$m, all.fits$sd, all.fits$tau, all.fits$standardize, all.fits$family )
    names(ret) <- c('null','fitted.models','model.indicator', 'cnames', 'm', 'sd', 'tau', 'standardize', 'family' )
    return( ret )
}

cv.prems <- function( y, x, x.fixed=NULL, no.cores=10, k.min=1, k.max, tau.i=NULL,
                     max.s=50, max2way='all', standardize=TRUE, nfolds=NULL, foldid=NULL,
                     n.waic=100, n.coef=1, lasso.factor=1,
                     criteria='aic', fit='mode', family='binomial', verbose=TRUE ){
    if( is.null(foldid) & is.null(nfolds) ){
        nfolds <- length(y)
        foldid <- 1:nfolds
    }
    if( !is.null(foldid) & is.null(nfolds) ){
        nfolds <- length(unique(foldid))
    }
    if( is.null(foldid) ){
        yy <- y
        if( family=="cox" )
            yy <- y[,3]
        if( family=="cox" | family=="binomial" )
            foldid <- make.folds2( yy, folds=nfolds )
        if( family=="gaussian" )
            foldid <- make.folds.continuous( length(y), nfolds  )
    }
    print( table( foldid ) )

    ll <- vector()
#    pred <- matrix(ncol=(k.max-k.min+1),nrow=length(y))
    pwll <- matrix( nrow=nfolds, ncol=(k.max-k.min+1) )
    if( verbose ){
        print(paste(nfolds,'fold cross-validation'))
    }
    for( i in 1:nfolds ){
        train <- which( foldid!=i )
        test <- which( foldid==i )
        if( is.null(tau.i) ){
            tau.est <- TauEst( y[train], x=x[train,], x.fixed=x.fixed[train,], family=family,
                          standardize=standardize, n.coef=n.coef, nfolds=10 )
#        tau <- tau.est$tau.1se + lasso.penalty*( tau.est$tau.opt - tau.est$tau.1se )
            tau <- tau.est$tau.opt * lasso.factor
            print(paste0("tau=",tau))
        }else{
            tau=tau.i
        }
        my.fit <- prems( y=y[train], x=x[train,], x.fixed=x.fixed[train,,drop=FALSE],
                        family=family, tau=tau, k.max=k.max, max.s=max.s,
                        standardize=standardize, max2way=max2way,
                        no.cores=no.cores, verbose=FALSE )
#        if( criteria=='waic' | fit=='mode' ){
        if( criteria=='waic'  | fit=='mean' ){
            my.fit <- fillICs( fitted.models=my.fit, y=y[train], x=x[train,],
                              n.waic=n.waic, model.sizes=k.min:k.max, no.cores=no.cores, verbose=FALSE )
        }
        for( k in k.min:k.max ){
            kk <- k - k.min + 1
            pred <- predict.prems( my.fit,
                                  newx=x[test,,drop=FALSE],
                                  newx.fixed=x.fixed[test,,drop=FALSE],
                                  size=k, family=family, criteria=criteria, fit=fit )
            if( family=='binomial' ){
                lp1 <- log(pred)
                lp1 <- ifelse( is.finite(lp1), lp1, -1000 )
                lp0 <- log(1-pred)
                lp0 <- ifelse( is.finite(lp0), lp0, -1000 )
                pwll[i,kk] <- sum( y[test]*lp1 + (1-y[test])*lp0 )
            }else if( family=='gaussian' ){
                pwll[i,kk] <- sum ( -( y[test] - pred )^2 )
            }else if( family=='cox' ){
                pwll[i,kk] <- cox_pl_left_trunc( y[test], pred )
            }
        }
        if( verbose ){
            print(paste('Fold',i,'complete.'))
        }
    }
    cvm <- apply( pwll, 2, mean )
    cvsd <- apply( pwll, 2, sd )/sqrt(nfolds)
    names(cvm) <- k.min:k.max
    names(cvsd) <- k.min:k.max

    sizes <- (k.min:k.max)[prems.optim( cvm, cvsd )]

    ret <- list( sizes[1], sizes[2], cvm, cvsd )
    names(ret) <- c('best','one.se','cvm','cvsd' )
    return( ret )
}

prems.optim <- function( cvm, cvsd ){
    k <- 1:length(cvm)
    best <- order(cvm,decreasing=TRUE)[1]
    ptr2 <- which((cvm+cvsd[best])[1:(best-1)]>cvm[best])
    if( length(ptr2)>0 ){
        one.se <- min(ptr2)
    }
    if( length(ptr2)==0 ){
        one.se <- best
    }
    return(c(best,one.se))
}

new.optim <- function( cvm, cvsd ){
    k <- 1:length(cvm)
    one.se <- vector()
    best <- order( cvm-cvsd, decreasing=TRUE )[1]
    for( i in 1:best ){
        ptr2 <- which((cvm+cvsd[i])[1:(i-1)]>cvm[i])
        if( length(ptr2)>0 ){
            one.se[i] <- k[min(ptr2)]
        }
        if( length(ptr2)==0 ){
            one.se[i] <- i
        }
    }
    best <- max(which(one.se==max(one.se)))
    return(c(best,one.se[best]))
}

TauEst <- function( y, x, x.fixed=NULL, family='binomial', standardize=TRUE,
                   n.coef=1, fit=NULL, nfolds=NULL ){
    if( is.null(nfolds) ){
        nfolds <- length(y)
    }
    if( is.null(fit) ){
        if( !is.null(x.fixed) ){
            lambda.factor <- c( rep(0,ncol(x.fixed)), rep(1,ncol(x)) )
        }else{
            lambda.factor <- rep(1,ncol(x))
            x.fixed <- matrix( ncol=0, nrow=length(y) )
        }
        fit <- cv.glmnet( x=as.matrix(cbind(x.fixed,x)), y=y, penalty.factor=lambda.factor,
                     family=family, alpha=1, nfolds=nfolds,
                     type.measure='deviance', grouped=FALSE, standardize=standardize )
    }

    beta <- getCoefGlmnet( fit, s='lambda.min' )[-(1:(1+ncol(x.fixed)))]
    if( length(beta) > 0 ){
        s <- rep(1,ncol(x))
        if( standardize ){
            ptr <- match( names(beta), colnames(x) )
            s <- apply( x[,ptr,drop=FALSE], 2, sd )
        }
        lambda <- fit$lambda.min * length(y)
        beta1 <- sort(abs(beta*s),decreasing=TRUE)
        n.coef <- ifelse( n.coef>length(beta), length(beta), n.coef )
        ptr <- 1:n.coef
        tau.opt <- lambda * sum(beta1[ptr]) / sum(beta1[ptr]^2)
    }else{
        tau.opt <- NA
    }

    beta <- getCoefGlmnet( fit, s='lambda.1se' )[-(1:(1+ncol(x.fixed)))]
    if( length(beta) > 0 ){
        s <- rep(1,length(beta))
        if( standardize ){
            ptr <- match( names(beta), colnames(x) )
            s <- apply( x[,ptr,drop=FALSE], 2, sd )
        }
        lambda <- fit$lambda.1se * length(y)
        beta1 <- sort(abs(beta*s),decreasing=TRUE)
        n.coef <- ifelse( n.coef>length(beta), length(beta), n.coef )
        ptr <- 1:n.coef
        tau.1se <- lambda * sum(beta1[ptr]) / sum(beta1[ptr]^2)
    }else{
        tau.1se <- NA
    }

    ret <- list( tau.opt, tau.1se, fit )
    names(ret) <- c('tau.opt', 'tau.1se', 'fit.lasso')
    return(ret)
}

my.auc <- function( my.fit, sizes, X, y, rank=1, criteria='aic' ){
    my.pred <- matrix(ncol=length(sizes),nrow=length(y))
    for( i in sizes ){
        ii <- i - min(sizes) + 1
        my.pred[,ii] <- predict.prems( my.fit, as.matrix(X), size=i, rank=rank, no.cores=10, criteria=criteria, fit='mode', family='binomial' )
    }
#    r <- matrix(ncol=3,nrow=length(sizes))
    r <- list()
    for( i in 1:length(sizes) ){
        r[[i]] <- roc( y, my.pred[,i], ci=TRUE )
#        r[i,] <- as.numeric(roc( y, my.pred[,i], ci=TRUE )$ci)
    }
    names(r) <- sizes
    return(r)
}

make.folds2 <- function( strata, folds ){
    which.fold <- rep(NA,length=sum(!is.na(strata)))
    strata <- as.factor(strata)
    for( ll in levels(strata) ){
        ptr.ll <- which( strata==ll )
        fold <- c( rep( 1:folds, times=floor(length(ptr.ll)/folds) ), sample( 1:folds, length(ptr.ll)%%folds) )
        which.fold[sample( ptr.ll, replace=FALSE )] <- fold
    }
    return(which.fold)
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
