.surv.difference <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1, conf.band=TRUE, band.end.pts=c(0,Inf), conf.level=.95) {
    logit <- function(x) log(x / (1-x))
    logit.prime <- function(x) 1/(x * (1-x))
    expit <- function(x) 1/(1 + exp(-x))

    n <- nrow(IF.vals.0)
    quant <- qt(1-(1-conf.level)/2, n - 1)#qnorm(1-(1-conf.level)/2)
    df <- data.frame(time=c(0,fit.times), surv.diff=c(0,surv.1-surv.0))

    diff <- surv.1 - surv.0
    IF.diff <- IF.vals.1 - IF.vals.0
    se.diff <- sqrt(colMeans(IF.diff^2))
    se.diff[se.diff == 0] <- NA

    in.bounds <- -1 < diff & diff < 1
    logit.diff <- logit((diff + 1)/2)
    IF.diff.logit <- IF.diff
    for(j in 1:length(diff)) IF.diff.logit[,j] <- IF.diff[,j] * logit.prime((diff[j] + 1) / 2) / 2

    se.diff.logit <- sqrt(colMeans(IF.diff.logit^2))
    se.diff.logit[se.diff.logit == 0] <- NA

    df$se <- c(0,se.diff)/sqrt(n)
    df$se.logit <- c(0, se.diff.logit / sqrt(n))
    ll <- 2 * expit(logit.diff - quant * se.diff.logit / sqrt(n)) - 1
    ul <- 2 * expit(logit.diff + quant * se.diff.logit / sqrt(n)) - 1
    # ll[!in.bounds] <- pmin(approx(fit.times[in.bounds], ll[in.bounds], xout=fit.times[!in.bounds], rule=2)$y, diff[!in.bounds])
    # ul[!in.bounds] <- pmax(approx(fit.times[in.bounds], ul[in.bounds], xout=fit.times[!in.bounds], rule=2)$y, diff[!in.bounds])
    df$ptwise.lower <- c(0, ll)#pmax(est - quant * res$se, 0)
    df$ptwise.upper <- c(0, ul) #pmin(est + quant * res$se, 1)

    # df$ptwise.lower <- df$surv.diff - quant * df$se
    # df$ptwise.upper <- df$surv.diff + quant * df$se
    df$ptwise.pval <- c(1, pchisq(( logit.diff / (se.diff.logit / sqrt(n)) )^2, df=1, lower.tail = FALSE))


    if(conf.band) {
        unif.info <- .estimate.uniform.quantile(IF.diff.logit[,!is.na(se.diff.logit) & fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2]], conf.level)
        unif.quant <- unif.info$quantile

        ll <- 2 * expit(logit.diff - unif.quant * se.diff.logit / sqrt(n)) - 1
        ul <- 2 * expit(logit.diff + unif.quant * se.diff.logit / sqrt(n)) - 1
        # ll[!in.bounds] <- pmin(approx(fit.times[in.bounds], ll[in.bounds], xout=fit.times[!in.bounds], rule=2)$y, diff[!in.bounds])
        # ul[!in.bounds] <- pmax(approx(fit.times[in.bounds], ul[in.bounds], xout=fit.times[!in.bounds], rule=2)$y, diff[!in.bounds])
        df$unif.lower <- c(0, ll)#pmax(est - quant * res$se, 0)
        df$unif.upper <- c(0, ul) #pmin(est + quant * res$se, 1)

        df$unif.lower[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA
        df$unif.upper[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA

        #df$unif.pval <- c(1, pchisq((diff/( unif.quant * se.logit))^2, df=1, lower.tail = FALSE))
    }

    res <- list(surv.diff.df=df)
    if(conf.band) {
        res$surv.diff.unif.quant <- unif.quant
        res$surv.diff.sim.maxes <- unif.info$sim.maxes
    }
    return(res)
}

.surv.ratio <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1, conf.band=TRUE, band.end.pts=c(0, Inf),conf.level=.95) {
    n <- nrow(IF.vals.0)
    quant <- qt(1-(1-conf.level)/2, df = n - 1)
    df <- data.frame(time=c(0,fit.times), log.surv.ratio=c(0,log(surv.1) - log(surv.0)), surv.ratio=c(1,surv.1 / surv.0))

    IF.log.ratio <- IF.vals.1 / matrix(surv.1, nrow=nrow(IF.vals.1), ncol=ncol(IF.vals.1), byrow=TRUE) -
        IF.vals.0 / matrix(surv.0, nrow=nrow(IF.vals.0), ncol=ncol(IF.vals.0), byrow=TRUE)
    se.log.ratio <- sqrt(colMeans(IF.log.ratio^2))
    se.log.ratio[se.log.ratio == 0] <- NA
    se.log.ratio[is.infinite(se.log.ratio)] <- NA
    df$se.log.ratio <- c(0,se.log.ratio)/sqrt(n)
    df$ptwise.lower.log <- df$log.surv.ratio - quant * df$se.log.ratio
    df$ptwise.upper.log <- df$log.surv.ratio + quant * df$se.log.ratio
    df$ptwise.lower <- exp(df$ptwise.lower.log)
    df$ptwise.upper <- exp(df$ptwise.upper.log)
    df$ptwise.pval <- c(1,pchisq((df$log.surv.ratio[-1]/(df$se.log.ratio[-1]))^2, df=1, lower.tail = FALSE))


    if(conf.band) {
        log.unif.info <- .estimate.uniform.quantile(IF.log.ratio[,!is.na(se.log.ratio) & fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2]], conf.level)
        log.unif.quant <- log.unif.info$quantile
        df$unif.lower.log <- df$log.surv.ratio - log.unif.quant * df$se.log.ratio
        df$unif.upper.log <- df$log.surv.ratio + log.unif.quant * df$se.log.ratio
        df$unif.lower <- exp(df$unif.lower.log)
        df$unif.upper <- exp(df$unif.upper)

        df$unif.lower[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA
        df$unif.upper[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA

      #  df$unif.pval <- c(1, pchisq((df$log.surv.ratio[-1]/( log.unif.quant * df$se.log.ratio[-1]))^2, df=1, lower.tail = FALSE))
    }

    res <- list(surv.ratio.df=df)
    if(conf.band) {
        res$log.surv.ratio.unif.quant <- log.unif.quant
        res$log.surv.ratio.sim.maxes <- log.unif.info$sim.maxes
    }

    return(res)
}

.risk.ratio <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1, conf.band=TRUE, band.end.pts=c(0,Inf), conf.level=.95) {
    n <- nrow(IF.vals.0)
    quant <- qt(1-(1-conf.level)/2, df = n - 1)
    risk.0 <- 1-surv.0
    risk.1 <- 1-surv.1
    IF.vals.0 <- -IF.vals.0
    IF.vals.1 <- -IF.vals.1
    df <- data.frame(time=c(0,fit.times), log.risk.ratio=c(NA,log(risk.1) - log(risk.0)), risk.ratio=c(NA,risk.1 / risk.0))

    IF.log.ratio <- IF.vals.1 / matrix(risk.1, nrow=nrow(IF.vals.1), ncol=ncol(IF.vals.1), byrow=TRUE) -
        IF.vals.0 / matrix(risk.0, nrow=nrow(IF.vals.0), ncol=ncol(IF.vals.0), byrow=TRUE)
    se.log.ratio <- sqrt(colMeans(IF.log.ratio^2))
    se.log.ratio[se.log.ratio == 0] <- NA
    se.log.ratio[is.infinite(se.log.ratio)] <- NA
    df$se.log.ratio <- c(NA,se.log.ratio)/sqrt(n)
    df$ptwise.lower.log <- df$log.risk.ratio - quant * df$se.log.ratio
    df$ptwise.upper.log <- df$log.risk.ratio + quant * df$se.log.ratio
    df$ptwise.lower <- exp(df$ptwise.lower.log)
    df$ptwise.upper <- exp(df$ptwise.upper.log)
    df$ptwise.pval <- c(1,pchisq((df$log.risk.ratio[-1]/(df$se.log.ratio[-1]))^2, df=1, lower.tail = FALSE))


    if(conf.band) {
        log.unif.info <- .estimate.uniform.quantile(IF.log.ratio[,!is.na(se.log.ratio) & fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2]], conf.level)
        log.unif.quant <- log.unif.info$quantile
        df$unif.lower.log <- df$log.risk.ratio - log.unif.quant * df$se.log.ratio
        df$unif.upper.log <- df$log.risk.ratio + log.unif.quant * df$se.log.ratio
        df$unif.lower <- exp(df$unif.lower.log)
        df$unif.upper <- exp(df$unif.upper)

        df$unif.lower[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA
        df$unif.upper[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA

       # df$unif.pval <- c(1, pchisq((df$log.risk.ratio[-1]/( log.unif.quant * df$se.log.ratio[-1]))^2, df=1, lower.tail = FALSE))
    }

    res <- list(risk.ratio.df=df)
    if(conf.band) {
        res$log.risk.ratio.unif.quant <- log.unif.quant
        res$log.risk.ratio.sim.maxes <- log.unif.info$sim.maxes
    }

    return(res)
}

.nnt <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1, conf.band=TRUE, band.end.pts=c(0, Inf), conf.level=.95) {
    n <- nrow(IF.vals.0)
    quant <- qt(1-(1-conf.level)/2, df = n - 1)
    df <- data.frame(time=c(0,fit.times), nnt=c(NA,1/(surv.1-surv.0)))#, log.nnt=c(NA, -log(surv.1-surv.0)))

    # IF.log.nnt <- -(IF.vals.1 - IF.vals.0)/matrix(surv.1 - surv.0, nrow=nrow(IF.vals.1), ncol=ncol(IF.vals.1), byrow=TRUE)
    # se.log.nnt <- sqrt(colMeans(IF.log.nnt^2))
    # se.log.nnt[se.log.nnt == 0] <- NA
    # se.log.nnt[is.infinite(se.log.nnt)] <- NA
    # df$se.log.nnt <- c(0,se.log.nnt)/sqrt(n)
    # df$ptwise.lower.log <- df$log.nnt - quant * df$se.log.nnt
    # df$ptwise.upper.log <- df$log.nnt + quant * df$se.log.nnt
    IF.nnt <- (IF.vals.0 - IF.vals.1) / matrix((surv.1 - surv.0)^2, nrow=nrow(IF.vals.0), ncol=ncol(IF.vals.0), byrow=TRUE)
    se.nnt <- sqrt(colMeans(IF.nnt^2))
    se.nnt[se.nnt == 0] <- NA

    df$se <- c(0,se.nnt)/sqrt(n)
    df$ptwise.lower <- df$nnt - quant * df$se
    df$ptwise.upper <- df$nnt + quant * df$se

    if(conf.band) {
        unif.info <- .estimate.uniform.quantile(IF.nnt[,!is.na(se.nnt) & !is.nan(se.nnt) & !is.infinite(se.nnt) & fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2]], conf.level)
        unif.quant <- unif.info$quantile
        df$unif.lower <- df$nnt - unif.quant * df$se
        df$unif.upper <- df$nnt + unif.quant * df$se
        df$unif.lower[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA
        df$unif.upper[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA

    }

    res <- list(nnt.df=df)
    if(conf.band) {
        res$nnt.unif.quant <- unif.quant
        res$nnt.sim.maxes <- unif.info$sim.maxes
    }
    return(res)
}
