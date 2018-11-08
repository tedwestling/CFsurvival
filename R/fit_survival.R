#' Estimate counterfactual survival functions
#'
#' This function estimates counterfactual survival functions and contrasts from right-censored data subject to potential confounding.
#'
#' For a binary treatment measured at baseline, the counterfactual survival through time \code{t} for treatment level \code{a} is defined as the probability of survival up to time \code{t} were all units in the population assigned to treatment \code{a}. This quantity is estimated using an augmented inverse probability of treatment/censoring weighted (AIPTW) estimator. Our estimator of the counterfactual survival function requires further estimation of three nuisance parameters: the conditional probability of receiving treatment given a set of potential confounders (i.e. the treatment propensity), the conditional survival function of the event given confounders, and the conditional survival function of censoring given confounders. At the moment, we provide two ways of estimating each of these nuisance parameters. The treatment propensity can be estimated either using a logistic regression or SuperLearner (van der Laan and Polley, 2006). The conditional survival functions of event and censoring can be estimated using Cox proportional hazard models or survival random forests (Ishwaran et al.). Additionally, the estimates of any of these parameters can be passed in directly if another method is desired.
#'
#'  \code{CFsurvival} by default returns 95\% pointwise confidence intervals and 95\% uniform confidence bands, both of which are based on the influence function of the estimator. Contrasts of the counterfactual survival under treatment and control can also be obtained, including the survival difference, survival ratio, risk ratio, and number needed to treat. Pointwise confidence intervals and uniform confidence bands for these contrasts will also be computed.
#'
#' @param time n x 1 numeric vector of observed right-censored follow-up times; i.e. the minimum of the event and censoring times.
#' @param event n x 1 numeric vector of status indicators of whether an event was observed.
#' @param treat n x 1 numeric vector of observed binary treatment/exposure group; either 0 or 1.
#' @param fit.times k x 1 numeric vector of time grid at which the counterfactual survival function is desired. Only values > 0 are allowed. Defauls to all unique positive values in \code{time}.
#' @param fit.treat Optional subset of \code{c(0,1)} for which the counterfactual survival curves are desired. For example, if \code{fit.treat=c(0,1)} (default behavior) then both placebo and treatment curves will be estimated, and if \code{fit.treat=0} then only the placebo curve will be estimated.
#' @param cond.surv.method The method that should be used to estimate the conditional survival and censoring probabilities. One of \code{randomForestSRC}, \code{coxph}, or \code{NULL}. If \code{randomForestSRC} is specified, the \code{randomForestSRC} package will be used to estimate random survival forests, and must be installed. If \code{coxph} is specified, the \code{coxph} function from the \code{survival} package will be used. If \code{cond.surv.method = NULL}, then \code{S.hats.0}, \code{G.hats.0}, and/or \code{S.hats.1} and \code{G.hats.1} must be specified.
#' @param confounders n x p numeric matrix of potential confounders to use when estimating the conditional survival probabilities. Missing values are not allowed
#' @param cens.subset Optional numeric vector of the columns of \code{confounders} that should be used to estimate the conditional censoring probabilities. Defaults to using all availalbe confounders to estimate the censoring probabilities.
#' @param propensity.method The method that should be used to estimate the propensity (i.e. probability of treatment) model. One of \code{glm}, \code{SuperLearner}, or \code{NULL}. If \code{glm} is specified, then a logistic regression will be used; if \code{SuperLearner} is specified, then the \code{SuperLearner} package will be used (and needs to be installed). If \code{propensity.method = NULL}, then \code{g.hats} must be provided.
#' @param treat.subset Optional numeric vector of the columns of \code{confounders} that should be used to estimate the propensities. Defaults to using all available confounders to estimate the propensities.
#' @param SL.library If \code{propensity.method = SuperLearner}, the library to use for the SuperLearner.
#' @param S.hats.0 Optional n x k matrix of estimates of the conditional survival of the event given treatmet = 0 and confounders. If \code{propensity.method = NULL} and \code{0 \%in\% fit.treat}, then \code{S.hats.0} must be specified.
#' @param S.hats.1 Optional n x k matrix of estimates of the conditional survival of the event given treatmet = 1 and confounders. If \code{propensity.method = NULL} and \code{1 \%in\% fit.treat}, then \code{S.hats.1} must be specified.
#' @param G.hats.0 Optional n x k matrix of estimates of the conditional survival of censoring given treatmet = 0 and confounders. If \code{propensity.method = NULL} and \code{0 \%in\% fit.treat}, then \code{G.hats.0} must be specified.
#' @param G.hats.1 Optional n x k matrix of estimates of the conditional survival of censoring given treatmet = 1 and confounders. If \code{propensity.method = NULL} and \code{1 \%in\% fit.treat}, then \code{G.hats.1} must be specified.
#' @param g.hats Optional n x 1 numeric vector of estimated propensities. If \code{propensity.method = NULL}, then \code{g.hats} must be specified.
#' @param conf.band Logical indicating whether to compute simultaneous confidence bands.
#' @param conf.level Desired coverage of confidence intervals/bands.
#' @param surv.diffs Logical indicating whether to return an estimate of the difference in the survival functions, along with confidence intervals and tests.
#' @param surv.ratios Logical indicating whether to return an estimate of the ratio in the survival functions, along with confidence intervals and tests.
#' @param risk.ratios Logical indicating whether to return an estimate of the difference in the survival functions, along with confidence intervals and tests.
#' @param nnt Logical indicating whether to return an estimate of the number needed to treat (nnt), along with confidence intervals.
#' @param verbose Logical indicating whether progress should be printed to the command line.
#' @param ... Additional arguments to be passed on to \code{randomForestSRC::rfsrc} (only used if \code{cond.surv.method = randomForestSRC}).
#' @return \code{CFsurvfit} returns a named list with the following elements:
#' \item{fit.times}{The time points at which the counterfactual survival curves (and contrasts) were fit.}
#' \item{surv.df}{A data frame with the estimated counterfactual survival functions and CIs.}
#' \item{IF.vals.0, IF.vals.1}{n x k matrices with the influence values for the counterfactual survival functions. Rows index observations, columns index time points.}
#' \item{surv.diff.df}{A data frame with the estimated counterfactual survival difference (treatment survival minus control survival), as well as confidence intervals and tests of the null hypothesis that the difference equals zero. Only returned if \code{surv.diffs=TRUE}.}
#' \item{surv.ratio.df}{A data frame with the estimated counterfactual survival ratio (treatment survival divided by control survival), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{surv.ratios=TRUE}.}
#' \item{risk.ratio.df}{A data frame with the estimated counterfactual risk ratio (treatment risk divided by control risk), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{risk.ratios=TRUE}.}
#' \item{nnt.df}{A data frame with the estimated counterfactual number needed to treat (one divided by the survival difference), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{nnt=TRUE}.}
#' \item{surv.0.sim.maxes, surv.1.sim.maxes, surv.diff.}{Samples from the approximate distribution of the maximum of the limiting Gaussian processes of the counterfactual survival functions.}
#' \item{g.hats}{The estimated treatment propensities.}
#' \item{S.hats.0, S.hats.1}{The estimated conditional survival functions of the event.}
#' \item{G.hats.0, G.hats.1}{The estimated conditional survival functions of censoring.}
#' \item{data}{The original time, event, and treatment data supplied to the function.}
#' @references Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. \emph{The Annals of Applied Statistics}, 2(3), 841-860.
#' @references van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). Super learner. \emph{Statistical Applications in Genetics and Molecular Biology}, 6(1).
#' @examples
#' # Define parameters
#' n <- 300
#' expit <- function(x) 1/(1 + exp(-x))
#' betaT <- 2; lambdaT <- 20; betaC <- 2; lambdaC <- 15
#'
#' # Define the true survival functions of the control (S0) and treatment (S1) groups
#' S0 <- function(t) sapply(t, function(t0) integrate(function(w) .5 * pweibull(t0, shape=betaT, scale=lambdaT * exp(-w-1), lower.tail = FALSE), lower=-1, upper=1)$value)
#' S1 <- function(t) sapply(t, function(t0) integrate(function(w) .5 * pweibull(t0, shape=betaT, scale=lambdaT * exp(-w), lower.tail = FALSE), lower=-1, upper=1)$value)
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
#' # It is a good idea to check the min/max of the propensity estimates and the min of the censoring estimates
#' # If they are very small, there may be positivity violations, which may result in invalid inference.
#' min(fit$g.hats)
#' max(fit$g.hats)
#' min(fit$G.hats.0)
#' min(fit$G.hats.1)
#'
#' # Plot the output
#' \dontrun{
#' library(ggplot2)
#'
#' # First plot the survival curves + conf intervals + conf bands
#' fit$surv.df$true.surv <- c(S1(c(0, fit$fit.times)), S0(c(0, fit$fit.times)))
#' ggplot(fit$surv.df) +
#'     geom_line(aes(time, true.surv, group=trt), color='black') +
#'     geom_line(aes(time, surv, color=as.factor(trt), group=trt)) +
#'     geom_line(aes(time, ptwise.lower, color=as.factor(trt), group=trt), linetype=2) +
#'     geom_line(aes(time, ptwise.upper, color=as.factor(trt), group=trt), linetype=2) +
#'     geom_line(aes(time, unif.lower, color=as.factor(trt), group=trt), linetype=3) +
#'     geom_line(aes(time, unif.upper, color=as.factor(trt), group=trt), linetype=3) +
#'     scale_color_discrete("Treatment") +
#'     xlab("Time") +
#'     ylab("Survival") +
#'     coord_cartesian(xlim=c(0,15), ylim=c(0,1))
#'
#' # Next plot the survival difference
#' fit$surv.diff.df$true.surv.diff <- c(S1(fit$surv.diff.df$time) - S0(fit$surv.diff.df$time))
#' ggplot(fit$surv.diff.df) +
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
#' fit$surv.ratio.df$true.surv.ratio <- c(S1(fit$surv.ratio.df$time) / S0(fit$surv.ratio.df$time))
#' ggplot(fit$surv.ratio.df) +
#'     geom_line(aes(time, true.surv.ratio), color='red') +
#'     geom_line(aes(time, surv.ratio)) +
#'     geom_line(aes(time, ptwise.lower), linetype=2) +
#'     geom_line(aes(time, ptwise.upper), linetype=2) +
#'     geom_line(aes(time, unif.lower), linetype=3) +
#'     geom_line(aes(time, unif.upper), linetype=3) +
#'     xlab("Time") +
#'     ylab("Survival ratio (treatment / control)") +
#'     coord_cartesian(xlim=c(0,15), ylim=c(0,10))
#'
#' # Risk ratio
#' fit$risk.ratio.df$true.risk.ratio <- (1 - S1(fit$risk.ratio.df$time)) / (1 - S0(fit$risk.ratio.df$time))
#' ggplot(fit$risk.ratio.df) +
#'     geom_line(aes(time, true.risk.ratio), color='red') +
#'     geom_line(aes(time, risk.ratio)) +
#'     geom_line(aes(time, ptwise.lower), linetype=2) +
#'     geom_line(aes(time, ptwise.upper), linetype=2) +
#'     geom_line(aes(time, unif.lower), linetype=3) +
#'     geom_line(aes(time, unif.upper), linetype=3) +
#'     coord_cartesian(xlim=c(0,15), ylim=c(0,1)) +
#'     xlab("Time") +
#'     ylab("Risk ratio (treatment / control)")
#'
#' # Number needed to treat
#' fit$nnt.df$true.nnt <- 1/(S1(fit$nnt.df$time) - S0(fit$nnt.df$time))
#' ggplot(fit$nnt.df) +
#'     geom_line(aes(time, true.nnt), color='red') +
#'     geom_line(aes(time, nnt)) +
#'     geom_line(aes(time, ptwise.lower), linetype=2) +
#'     geom_line(aes(time, ptwise.upper), linetype=2) +
#'     geom_line(aes(time, unif.lower), linetype=3) +
#'     geom_line(aes(time, unif.upper), linetype=3) +
#'     coord_cartesian(xlim=c(0,15), ylim=c(0, 10)) +
#'     xlab("Time") +
#'     ylab("Number needed to treat (NNT)")}



CFsurvfit <- function(time, event, treat, fit.times=sort(unique(time[time > 0 & time < max(time)])), fit.treat=c(0,1), cond.surv.method=c("randomForestSRC", "coxph", NULL), confounders=NULL, cens.subset=NULL, propensity.method=c("SuperLearner", "glm", NULL), treat.subset=NULL, SL.library=NULL, S.hats.0=NULL, G.hats.0=NULL, S.hats.1=NULL, G.hats.1=NULL, g.hats=NULL, conf.band=TRUE, conf.level=.95, surv.diffs=TRUE, surv.ratios=TRUE, risk.ratios=FALSE, nnt=FALSE, verbose=FALSE, ...) {
    .args <- mget(names(formals()),sys.frame(sys.nframe()))
    do.call(.check.input, .args)

    if(any(fit.times <= 0)) {
        fit.times <- fit.times[fit.times > 0]
        warning("fit.times <= 0 removed.")
    }
    if(any(fit.times > max(time))) {
        fit.times <- fit.times[fit.times <= max(time)]
        warning("fit.times > max(time) removed.")
    }

    #### ESTIMATE PROPENSITY ####

    if(is.null(g.hats)) {
        if(verbose) message("Estimating propensities...")
        if(is.null(treat.subset)) treat.subset <- 1:ncol(confounders)
        g.fit <- .estimate.propensity(A=treat, W.propensity=confounders[,treat.subset], method=propensity.method, SL.library=SL.library)
        g.hats <- g.fit$g.hats
    }

    #### ESTIMATE EVENT SURVIVALS ####

    if((1 %in% fit.treat & is.null(S.hats.1)) | (0 %in% fit.treat & is.null(S.hats.0))) {
        if(verbose) message("Estimating conditional event survivals...")
        S.hats <- .estimate.conditional.survival(Y=time, Delta=event, A=treat, fit.times=fit.times, fit.treat=fit.treat, method=cond.surv.method, W=confounders, ...)
        if(1 %in% fit.treat & is.null(S.hats.1)) {
            S.hats.1 <- S.hats$S.hats.1
        }
        if(0 %in% fit.treat & is.null(S.hats.0)) {
            S.hats.0 <- S.hats$S.hats.0
        }
    }

    #### ESTIMATE CENSORING SURVIVALS ####

    if((1 %in% fit.treat & is.null(G.hats.1)) | (0 %in% fit.treat & is.null(G.hats.0))) {
        if(verbose) message("Estimating conditional censoring survivals...")
        G.hats <- .estimate.conditional.survival(Y=time, Delta=1-event, A=treat, fit.times=fit.times, fit.treat=fit.treat, method=cond.surv.method, W=confounders[,cens.subset], ...)
        if(1 %in% fit.treat & is.null(G.hats.1)) {
            G.hats.1 <- G.hats$S.hats.1
        }
        if(0 %in% fit.treat & is.null(G.hats.0)) {
            G.hats.0 <- G.hats$S.hats.0
        }
    }

    #### ESTIMATE CF SURVIVALS ####

    surv.df <- data.frame()
    result <- list(fit.times=fit.times, fit.treat=fit.treat, surv.df=surv.df)
    if(verbose) message("Computing counterfactual survivals...")
    if(1 %in% fit.treat) {
        surv.1 <- .get.survival(Y=time, Delta=event, A=treat, times=fit.times, S.hats=S.hats.1, G.hats=G.hats.1, g.hats=g.hats)
        surv.df.1 <- data.frame(time=c(0,fit.times), trt=1, surv=c(1, surv.1$surv.iso))
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
        surv.0 <- .get.survival(Y=time, Delta=event, A=1-treat, times=fit.times, S.hats=S.hats.0, G.hats=G.hats.0, g.hats=1-g.hats)
        surv.df.0 <- data.frame(time=c(0,fit.times), trt=0, surv=c(1, surv.0$surv.iso))
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

    result$data <- data.frame(time, event, treat)

    return(result)

}

.check.input <- function(time, event, treat, fit.times, fit.treat, cond.surv.method, confounders, cens.subset, propensity.method, treat.subset, SL.library, S.hats.0, G.hats.0, S.hats.1, G.hats.1, g.hats, conf.band, conf.level, surv.diffs, surv.ratios, risk.ratios, nnt, verbose, ...) {
    if(any(time < 0)) stop("Only positive event/censoring times allowed!")
    if(any(time == 0 & event == 1)) stop("Events at time zero not allowed.")
    if(any(!(event %in% c(0,1)))) stop("Event must be binary.")
    if(any(!(treat %in% c(0,1)))) stop("Treatment must be binary.")
    if(length(time) != length(event) | length(time) != length(treat)) stop("time, event, and treat must be n x 1 vectors")
    if(any(!(fit.treat %in% c(0,1)))) stop("fit.treat must be a subset of c(0,1).")
    if(!is.null(S.hats.1) && dim(S.hats.1) != c(nrow(time), length(fit.times))) {
        stop("S.hats must be an n x k matrix (n is number of observations, k is length of fit.times)")
    }
    if(!is.null(S.hats.0) && dim(S.hats.0) != c(nrow(time), length(fit.times))) {
        stop("S.hats must be an n x k matrix (n is number of observations, k is length of fit.times)")
    }
    if(!is.null(G.hats.1) && dim(G.hats.1) != c(nrow(time), length(fit.times))) {
        stop("G.hats must be an n x k matrix (n is number of observations, k is length of fit.times)")
    }
    if(!is.null(G.hats.0) && dim(G.hats.0) != c(nrow(time), length(fit.times))) {
        stop("G.hats must be an n x k matrix (n is number of observations, k is length of fit.times)")
    }
    if(!is.null(cond.surv.method)) {
        if(!(cond.surv.method %in% c("coxph", "randomForestSRC"))) {
            stop("cond.surv.method must be one of coxph or randomForestSRC if not NULL")
        }
        if(!is.null(S.hats.0) | !is.null(G.hats.0) | !is.null(S.hats.1) | !is.null(G.hats.1)) {
            warning("cond.surv.method not NULL but S.hats or G.hats provided.")
        }
    } else {
        if(0 %in% fit.treat & (is.null(S.hats.0) | is.null(G.hats.0))) {
            stop("S.hats.0 and G.hats.0 if cond.surv.method is not specified and 0 is a treatment of interest.")
        }
        if(1 %in% fit.treat & (is.null(S.hats.1) | is.null(G.hats.1))) {
            stop("S.hats.1 and G.hats.1 if cond.surv.method is not specified and 1 is a treatment of interest.")
        }
    }
    if((surv.diffs | surv.ratios | risk.ratios | nnt ) & !identical(sort(fit.treat), c(0,1))) {
        warning("surv.diffs, surv.ratios, or risk ratios specified but both treatment regimens not requested -- contrasts will not be provided. Re-run with fit.treat = c(0,1) for survival contrasts.")
    }

}



.get.survival <- function(Y, Delta, A, times, S.hats, G.hats, g.hats, isotonize=TRUE) {
    times <- times[times > 0]
    n <- length(Y)
    ord <- order(times)
    times <- times[ord]
    S.hats <- S.hats[,ord]
    G.hats <- G.hats[,ord]

    int.vals <- t(sapply(1:n,function(i) {
        vals <- diff(c(1,1/S.hats[i,]))* 1/ G.hats[i,]
        vals[times >= Y[i]] <- 0
        cumsum(vals)
    }))
    S.hats.Y <- sapply(1:n, function(i) stepfun(times, c(1,S.hats[i,]))(Y[i]))
    G.hats.Y <- sapply(1:n, function(i) stepfun(times, c(1,G.hats[i,]))(Y[i]))
    IF.vals <- matrix(NA, nrow=n, ncol=length(times))
    surv <- rep(NA, length(times))
    for(k in 1:length(times)) {
        t0 <- times[k]
        S.hats.t0 <- S.hats[,k]
        inner.func.1 <- ifelse(Y <= t0 & Delta == 1, 1/(S.hats.Y * G.hats.Y), 0 )
        inner.func.2 <- int.vals[,k]
        if.func <- as.numeric(A == 1) * S.hats.t0 * ( -inner.func.1 + inner.func.2) / g.hats + S.hats.t0
        surv[k] <- mean(if.func)
        IF.vals[,k] <- if.func - surv[k]
    }
    res <- list(times=times, surv=pmin(1,pmax(0,surv)), IF.vals=IF.vals)
    if(isotonize) {
        res$surv.iso <- NA
        res$surv.iso[!is.na(res$surv)] <- 1 - isoreg(res$times[!is.na(res$surv)], 1-res$surv[!is.na(res$surv)])$yf
    }

    return(res)
}

.surv.confints <- function(times, est, IF.vals, isotonize=TRUE, conf.band=TRUE, conf.level=.95) {
    n <- nrow(IF.vals)
    res <- NULL
    res$se <- sqrt(colMeans(IF.vals^2)) / sqrt(n)
    res$se[res$se == 0] <- NA
    quant <- qnorm(1-(1-conf.level)/2)

    # Raw intervals and bands based on un-isotonized survival ests
    res$ptwise.lower <- pmax(est - quant * res$se, 0)
    res$ptwise.upper <- pmin(est + quant * res$se, 1)

    # Isotonized intervals and bands
    if(isotonize) {
        res$ptwise.lower <- NA
        res$ptwise.lower[!is.na(res$ptwise.lower)] <- 1 - isoreg(times[!is.na(res$ptwise.lower)], 1-res$ptwise.lower[!is.na(res$ptwise.lower)])$yf
        res$ptwise.upper <- NA
        res$ptwise.upper[!is.na(res$ptwise.upper)] <- 1 - isoreg(times[!is.na(res$ptwise.upper)], 1-res$ptwise.upper[!is.na(res$ptwise.upper)])$yf
    }
    out <- NULL
    if(conf.band) {
        unif.vals <- .estimate.uniform.quantile(IF.vals[,!is.na(res$se)], conf.level)
        unif.quant <- unif.vals$quantile
        out$sim.maxes <- unif.vals$maxes
        out$unif.quant <- unif.quant
        res$unif.lower <- pmax(est - unif.quant * res$se, 0)
        res$unif.upper <- pmin(est + unif.quant * res$se, 1)
        if(isotonize) {
            res$unif.lower <- 1 - isoreg(times[!is.na(res$unif.lower)], 1-res$unif.lower[!is.na(res$unif.lower)])$yf
            res$unif.upper <- 1 - isoreg(times[!is.na(res$unif.upper)], 1-res$unif.upper[!is.na(res$unif.upper)])$yf
        }
    }
    out$res <- res
    return(out)
}

.estimate.uniform.quantile <- function(IF.vals, conf.level=.95) {
    n <- nrow(IF.vals)
    IF.vals <- scale(IF.vals)
    maxes <- replicate(1e4, max(abs(rbind(rnorm(n)/sqrt(n)) %*% IF.vals)))
    return(list(quantile=quantile(maxes, 1 - (1 - conf.level)/2), maxes=maxes))
}

# conditional.cum.inc <- function(surv0, IF0, survs, IF.vals, conf.level=.95) {
#     n <- nrow(IF.vals)
#     dlog <- log(survs) - log(surv0)
#     logIF0 <- IF0 / surv0
#     logIF.vals <- IF.vals / matrix(survs, nrow=nrow(IF.vals), ncol=ncol(IF.vals), byrow=TRUE)
#     dlog.IF <- logIF.vals - logIF0
#
#     SEs <- apply(dlog.IF, 2, sd)
#     SEs[SEs == 0] <- NA
#     quant <- qnorm(1-(1-conf.level)/2)
#     nas <- is.na(SEs)
#     unif.info <- estimate.uniform.quantile(dlog.IF[,!nas], conf.level)
#     unif.quant <- unif.info$quantile
#
#     ptwise.lower <- dlog - quant * SEs / sqrt(n)
#     ptwise.upper <- dlog + quant * SEs / sqrt(n)
#     unif.lower <- dlog - unif.quant * SEs / sqrt(n)
#     unif.upper <- dlog + unif.quant * SEs / sqrt(n)
#
#     ptwise.lower.iso <- 1 - isoreg(1:sum(!nas), 1-ptwise.lower[!nas])$yf
#     ptwise.upper.iso <- 1 - isoreg(1:sum(!nas), 1-ptwise.upper[!nas])$yf
#
#     unif.lower.iso <- 1 - isoreg(1:sum(!nas), 1-unif.lower[!nas])$yf
#     unif.upper.iso <- 1 - isoreg(1:sum(!nas), 1-unif.upper[!nas])$yf
#
#     res <- data.frame(cond.cum.inc = pmin(pmax(1 - exp(dlog),0),1))
#     res$ptwise.lower <- res$ptwise.upper <- res$unif.lower <- res$unif.upper <- rep(NA, ncol(IF.vals))
#     res$ptwise.lower[!nas] <- pmax(1 - exp(ptwise.upper.iso), 0)
#     res$ptwise.upper[!nas] <- pmin(1 - exp(ptwise.lower.iso), 1)
#     res$unif.lower[!nas] <- pmax(1 - exp(unif.upper.iso), 0)
#     res$unif.upper[!nas] <- pmin(1 - exp(unif.lower.iso), 1)
#     res$log.ptwise.se <- quant * SEs / sqrt(n)
#     res$log.unif.se <- unif.quant * SEs / sqrt(n)
#     res$SE.log <- SEs
#
#     IF.vals <- -dlog.IF * matrix(exp(dlog), nrow=nrow(dlog.IF), ncol=ncol(dlog.IF), byrow=TRUE)
#     return(list(df=res, IF.vals=IF.vals, maxes=unif.info$maxes))
# }


