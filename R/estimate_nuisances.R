.estimate.conditional.survival <- function(Y, Delta, A, W, newW, fit.times, fit.treat, event.SL.library, cens.SL.library, survSL.control, survSL.cvControl, verbose, save.fit) {
    ret <- list(fit.times=fit.times)
    AW <- cbind(A, W)
    if(0 %in% fit.treat & 1 %in% fit.treat) {
        newAW <- rbind(cbind(A=0, newW), cbind(A=1, newW))
    } else {
        newAW <- cbind(A=fit.treat, newW)
    }
    res <- require(survSuperLearner)
    if(!res) stop("Please install the package survSuperLearner via:\n devtools::install_github('tedwestling/survSuperLearner')")
    if(is.null(survSL.control)) survSL.control <- list(saveFitLibrary = save.fit)

    fit <- survSuperLearner(time = Y, event = Delta,  X = AW, newX = newAW, new.times = fit.times, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose=verbose, control = survSL.control, cvControl = survSL.cvControl)
    if(save.fit) ret$surv.fit <- fit
    if(0 %in% fit.treat) {
        ret$event.pred.0 <- fit$event.SL.predict[1:nrow(newW),]
        if(any(ret$event.pred.0 == 0)) ret$event.pred.0[ret$event.pred.0 == 0] <- min(ret$event.pred.0[ret$event.pred.0 > 0])
        ret$cens.pred.0 <- fit$cens.SL.predict[1:nrow(newW),]
        if(any(ret$cens.pred.0 == 0)) ret$cens.pred.0[ret$cens.pred.0 == 0] <- min(ret$cens.pred.0[ret$cens.pred.0 > 0])
        if(1 %in% fit.treat) {
            ret$event.pred.1 <- fit$event.SL.predict[-(1:nrow(newW)),]
            if(any(ret$event.pred.1 == 0)) ret$event.pred.1[ret$event.pred.1 == 0] <- min(ret$event.pred.1[ret$event.pred.1 > 0])

            ret$cens.pred.1 <- fit$cens.SL.predict[-(1:nrow(newW)),]
            if(any(ret$cens.pred.1 == 0)) ret$cens.pred.1[ret$cens.pred.1 == 0] <- min(ret$cens.pred.1[ret$cens.pred.1 > 0])
        }
    } else {
        ret$event.pred.1 <- fit$event.SL.predict
        if(any(ret$event.pred.1 == 0)) ret$event.pred.1[ret$event.pred.1 == 0] <- min(ret$event.pred.1[ret$event.pred.1 > 0])
        ret$cens.pred.1 <- fit$cens.SL.predict
        if(any(ret$cens.pred.1 == 0)) ret$cens.pred.1[ret$cens.pred.1 == 0] <- min(ret$cens.pred.1[ret$cens.pred.1 > 0])
    }
    ret$event.coef <- fit$event.coef
    ret$cens.coef <- fit$cens.coef
    return(ret)
}

.estimate.propensity <- function(A, W, newW, SL.library, save.fit, verbose) {
    ret <- list()
    library(SuperLearner)
    if (length(SL.library) == 1) {
        if(length(unlist(SL.library)) == 2 & ncol(W) > 1) {
            screen <- get(SL.library[[1]][2])
            whichScreen <- screen(Y = A, X = W, family = 'binomial')
        } else {
            whichScreen <- rep(TRUE, ncol(W))
        }
        learner <- get(SL.library[[1]][1])
        prop.fit <- learner(Y = A, X = W[,whichScreen, drop=FALSE], newX = newW[,whichScreen, drop=FALSE], family='binomial')
        ret$prop.pred <- prop.fit$pred
        if(save.fit) {
            ret$prop.fit <- list(whichScreen = whichScreen, pred.alg = prop.fit)
        }
    } else {
        prop.fit <- SuperLearner(Y=A, X=W, newX=newW, family='binomial',
                                 SL.library=SL.library, method = "method.NNloglik", verbose = verbose)
        ret$prop.pred <- c(prop.fit$SL.predict)
        if(save.fit) ret$prop.fit <- prop.fit
    }
    return(ret)
}
