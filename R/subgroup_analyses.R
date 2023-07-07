#' Compute subgroup counterfactual survival functions
#'
#' This function estimates subgroup-specific counterfactual survival functions and contrasts from right-censored data subject to potential confounding using the output from the \code{\link{CFsurvfit}} function.
#'
#' See the documentation for \code{\link{CFsurvfit}} for an explanation of counterfactual survival functions and how they are estimated.
#'
#' @param fit Previously fitted CF survival curves from the \code{CFsurvfit}
#' @param subgroup.inds Indices of the original data set corresponding to the subgroup for which the counterfactual survivals and contrasts are desired.
#' @param conf.band Logical indicating whether to compute simultaneous confidence bands.
#' @param conf.level Desired coverage of confidence intervals/bands.
#' @param surv.diffs Logical indicating whether to return an estimate of the difference in the survival functions, along with confidence intervals and tests.
#' @param surv.ratios Logical indicating whether to return an estimate of the ratio in the survival functions, along with confidence intervals and tests.
#' @param risk.ratios Logical indicating whether to return an estimate of the difference in the survival functions, along with confidence intervals and tests.
#' @param nnt Logical indicating whether to return an estimate of the number needed to treat (nnt), along with confidence intervals.
#' @param verbose Logical indicating whether progress should be printed.
#' @return \code{subgroup.CFsurvfit} returns a list with the same structure as the output from \code{\link{CFsurvfit}} -- see the documentation of that function for additional details.
#' @examples
#' # Define parameters
#' n <- 300
#' expit <- function(x) 1/(1 + exp(-x))
#' betaT <- 2; lambdaT <- 20; betaC <- 2; lambdaC <- 15
#'
#' # Simulate data
#' covar <- runif(n, min=-1, max=1)
#' g0s <- expit(.2 - covar)
#' rx <- rbinom(n, size=1, prob=g0s)
#' event.time <- rweibull(n, shape = betaT, scale = lambdaT * exp(-covar - 1 + rx))
#' cens.time <- rweibull(n, shape = betaC, scale = lambdaC * exp(-covar/5 - rx/5))
#' cens.time[cens.time > 15] <- 15
#' obs.time <- pmin(event.time, cens.time)
#' obs.event <- as.numeric(event.time <= cens.time)
#'
#' # Estimate the CF survivals
#' fit <- CFsurvfit(time=obs.time, event=obs.event, treat=rx, cond.surv.method = "coxph", confounders = data.frame(covar), propensity.method = "glm", surv.diffs=TRUE, surv.ratios=TRUE, risk.ratios=TRUE, nnt=TRUE, verbose=TRUE)
#'
#' # Suppose we want the counterfactual survivals & contrasts among patients with covariate
#' #     values > 0
#' subgp <- which(covar > 0)
#' sub.fit <- subgroup.CFsurvfit(fit, subgp, verbose=TRUE)
#'
##' # Define the true conditional survival functions of the control (S0) and treatment (S1) groups among units with covar > 0
#' S0 <- function(t) sapply(t, function(t0) integrate(function(w) pweibull(t0, shape=betaT, scale=lambdaT * exp(-w-1), lower.tail = FALSE), lower=0, upper=1)$value)
#' S1 <- function(t) sapply(t, function(t0) integrate(function(w) pweibull(t0, shape=betaT, scale=lambdaT * exp(-w), lower.tail = FALSE), lower=0, upper=1)$value)
#'
#'
#' # Plot the output
#' \dontrun{
#' library(ggplot2)
#'
#' # First plot the survival curves + conf intervals + conf bands
#' sub.fit$surv.df$true.surv <- c(S1(c(0, sub.fit$fit.times)), S0(c(0, sub.fit$fit.times)))
#' ggplot(sub.fit$surv.df) +
#'     geom_line(aes(time, true.surv, group=trt), color='black') +
#'     geom_line(aes(time, surv, color=as.factor(trt), group=trt)) +
#'     geom_line(aes(time, ptwise.lower, color=as.factor(trt), group=trt), linetype=2) +
#'     geom_line(aes(time, ptwise.upper, color=as.factor(trt), group=trt), linetype=2) +
#'     geom_line(aes(time, unif.logit.lower, color=as.factor(trt), group=trt), linetype=3) +
#'     geom_line(aes(time, unif.logit.upper, color=as.factor(trt), group=trt), linetype=3) +
#'     scale_color_discrete("Treatment") +
#'     xlab("Time") +
#'     ylab("Survival") +
#'     coord_cartesian(xlim=c(0,15), ylim=c(0,1))
#'
#' # Next plot the survival difference
#' sub.fit$surv.diff.df$true.surv.diff <- c(S1(sub.fit$surv.diff.df$time) - S0(sub.fit$surv.diff.df$time))
#' ggplot(sub.fit$surv.diff.df) +
#'     geom_line(aes(time, true.surv.diff), color='red') +
#'     geom_line(aes(time, surv.diff)) +
#'     geom_line(aes(time, ptwise.lower), linetype=2) +
#'     geom_line(aes(time, ptwise.upper), linetype=2) +
#'     geom_line(aes(time, unif.lower), linetype=3) +
#'     geom_line(aes(time, unif.upper), linetype=3) +
#'     xlab("Time") +
#'     ylab("Survival difference (treatment - control)") +
#'     coord_cartesian(xlim=c(0,15), ylim=c(0,1))
#'
# Survival ratio
#' sub.fit$surv.ratio.df$true.surv.ratio <- c(S1(sub.fit$surv.ratio.df$time) / S0(sub.fit$surv.ratio.df$time))
#' ggplot(sub.fit$surv.ratio.df) +
#'     geom_line(aes(time, true.surv.ratio), color='red') +
#'     geom_line(aes(time, surv.ratio)) +
#'     geom_line(aes(time, ptwise.lower), linetype=2) +
#'     geom_line(aes(time, ptwise.upper), linetype=2) +
#'     geom_line(aes(time, unif.lower), linetype=3) +
#'     geom_line(aes(time, unif.upper), linetype=3) +
#'     xlab("Time") +
#'     ylab("Survival ratio (treatment / control)") +
#'     coord_cartesian(xlim=c(0,15), ylim=c(0,10))}


subgroup.CFsurvfit <- function(fit, subgroup.inds, conf.band=TRUE, conf.level=.95, surv.diffs=TRUE, surv.ratios=TRUE, risk.ratios=FALSE, nnt=FALSE, verbose=FALSE) {
    fit.times <- fit$fit.times
    time <- fit$data$time[subgroup.inds]
    event <- fit$data$event[subgroup.inds]
    treat <- fit$data$treat[subgroup.inds]
    fit.treat <- fit$fit.treat
    g.hats <- fit$g.hats[subgroup.inds]

    #### ESTIMATE CF SURVIVALS ####
    surv.df <- data.frame()
    result <- list(fit.times=fit.times, fit.treat=fit.treat, subgroup.inds=subgroup.inds, surv.df=surv.df)
    if(verbose) message("Computing counterfactual survivals...")
    if(1 %in% fit.treat) {
        S.hats.1 <-fit$S.hats.1[subgroup.inds,]
        G.hats.1 <-fit$G.hats.1[subgroup.inds,]
        surv.1 <- .get.survival(Y=time, Delta=event, A=treat, times=fit.times, S.hats=S.hats.1, G.hats=G.hats.1, g.hats=g.hats)
        surv.df.1 <- data.frame(time=c(0,fit.times), trt=1, surv=c(1, surv.1$surv))
        result$IF.vals.1 <- surv.1$IF.vals

        c.int <- .surv.confints(fit.times, surv.1$surv, surv.1$IF.vals, conf.band = conf.band, conf.level=conf.level)
        surv.df.1$se <- c(0,c.int$res$se)
        surv.df.1$ptwise.lower <- c(1,c.int$res$ptwise.lower)
        surv.df.1$ptwise.upper <- c(1,c.int$res$ptwise.upper)

        if(conf.band) {
            surv.df.1$unif.lower <- c(1,c.int$res$unif.lower)
            surv.df.1$unif.upper <- c(1,c.int$res$unif.upper)
            result$surv.1.unif.quant <- c.int$unif.quant
            result$surv.1.sim.maxes <- c.int$sim.maxes
        }
        result$surv.df <- rbind(result$surv.df, surv.df.1)
    }

    if(0 %in% fit.treat) {
        S.hats.0 <-fit$S.hats.0[subgroup.inds,]
        G.hats.0 <-fit$G.hats.0[subgroup.inds,]
        surv.0 <- .get.survival(Y=time, Delta=event, A=1-treat, times=fit.times, S.hats=S.hats.0, G.hats=G.hats.0, g.hats=1-g.hats)
        surv.df.0 <- data.frame(time=c(0,fit.times), trt=0, surv=c(1, surv.0$surv))
        result$IF.vals.0 <- surv.0$IF.vals

        c.int <- .surv.confints(fit.times, surv.0$surv, surv.0$IF.vals, conf.band = conf.band, conf.level=conf.level)
        surv.df.0$se <- c(0,c.int$res$se)
        surv.df.0$ptwise.lower <- c(1,c.int$res$ptwise.lower)
        surv.df.0$ptwise.upper <- c(1,c.int$res$ptwise.upper)

        if(conf.band) {
            surv.df.0$unif.lower <- c(1,c.int$res$unif.lower)
            surv.df.0$unif.upper <- c(1,c.int$res$unif.upper)
            result$surv.0.unif.quant <- c.int$unif.quant
            result$surv.0.sim.maxes <- c.int$sim.maxes
        }
        result$surv.df <- rbind(result$surv.df, surv.df.0)
    }

    #### ESTIMATE CF CONTRASTS ####
    if(surv.diffs & identical(sort(fit.treat), c(0,1))) {
        if(verbose) message("Computing survival differences...")
        surv.diff <- .surv.difference(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, conf.level=conf.level)
        result <- c(result, surv.diff)
    }

    if(surv.ratios & identical(sort(fit.treat), c(0,1))) {
        if(verbose) message("Computing survival ratios...")
        surv.ratio <- .surv.ratio(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, conf.level=conf.level)
        result <- c(result, surv.ratio)
    }

    if(risk.ratios & identical(sort(fit.treat), c(0,1))) {
        if(verbose) message("Computing risk ratios...")
        risk.ratio <- .risk.ratio(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, conf.level=conf.level)
        result <- c(result, risk.ratio)
    }

    if(nnt & identical(sort(fit.treat), c(0,1))) {
        if(verbose) message("Computing number needed to treat...")
        nnt <- .nnt(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, conf.level=conf.level)
        result <- c(result, nnt)
    }

    result$g.hats <- g.hats

    if(1 %in% fit.treat) {
        result$S.hats.1 <- S.hats.1
        result$G.hats.1 <- G.hats.1
    }
    if(0 %in% fit.treat) {
        result$S.hats.0 <- S.hats.0
        result$G.hats.0 <- G.hats.0
    }

    return(result)

}
