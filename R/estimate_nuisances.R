.estimate.conditional.survival <- function(Y, Delta, A, fit.times, fit.treat, method, W, ...) {
    ret <- list(times=fit.times)
    if(method == "coxph") {
        AW <- cbind(A, W)
        library(survival)
        df <- as.data.frame(cbind(Y=Y, Delta=Delta, A=A, W))
        fit <- coxph(Surv(Y, Delta) ~ ., data=df)
        if(0 %in% fit.treat) {
            AW0 <- as.data.frame(cbind(A=0, W))
            ret$S.hats.0 <- t(summary(survfit(fit, newdata=AW0,  se.fit = FALSE, conf.int = FALSE), times=fit.times)$surv)
        }
        if(1 %in% fit.treat) {
            AW1 <- as.data.frame(cbind(A=1, W))
            ret$S.hats.1 <- t(summary(survfit(fit, newdata=AW1,  se.fit = FALSE, conf.int = FALSE), times=fit.times)$surv)
        }
    }
    if(method == "randomForestSRC") {
        library(randomForestSRC)
        library(survival)
        data <- data.frame(Y, Delta, A)
        data <- cbind(data, W)
        fit <- rfsrc(Surv(Y, Delta) ~ ., data=data, ...)
        if(0 %in% fit.treat) {
            AW0 <- data
            AW0$A <- 0
            survs <- predict(fit, newdata=AW0, importance='none')$survival
            ret$S.hats.0 <- t(sapply(1:nrow(survs), function(i) {
                approx(c(0,fit$time.interest), c(0,survs[i,]), method='linear', xout = fit.times)$y
            }))
        }
        if(1 %in% fit.treat) {
            AW1 <- data
            AW1$A <- 1
            survs <- predict(fit, newdata=AW1, importance='none')$survival
            ret$S.hats.1 <- t(sapply(1:nrow(survs), function(i) {
                approx(c(0,fit$time.interest), c(0,survs[i,]), method='linear', xout = fit.times)$y
            }))
        }


    }
    return(ret)
}

.estimate.propensity <- function(A, W.propensity, method, SL.library=NULL) {
    ret <- list(method=method, SL.library=SL.library)
    if (method == "glm") {
        prop.fit <- glm(A ~ W.propensity, family='binomial')
        ret$g.hats <- prop.fit$fitted.values
        ret$propensity.coef <- prop.fit$coef
    }
    if (method == "SuperLearner") {
        library(SuperLearner)
        sl.fit <- SuperLearner(Y=A, X=as.data.frame(W.propensity), family='binomial',
                               SL.library=SL.library, method = "method.NNloglik")
        ret$g.hats <- sl.fit$SL.predict
        ret$SL.coef <- sl.fit$coef
    }
    return(ret)
}
