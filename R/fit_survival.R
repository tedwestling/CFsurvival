# TODO:
#  CIs/CBs for contrasts

#' Estimate counterfactual survival functions
#'
#' This function estimates counterfactual survival functions and contrasts from right-censored data subject to potential confounding.
#'
#' For a binary treatment measured at baseline, the counterfactual survival through time \code{t} for treatment level \code{a} is defined as the probability of survival up to time \code{t} were all units in the population assigned to treatment \code{a}. This quantity is estimated using an augmented inverse probability of treatment/censoring weighted (AIPTW) estimator. Our estimator of the counterfactual survival function requires further estimation of three nuisance parameters: the conditional probability of receiving treatment given a set of potential confounders (i.e. the treatment propensity), the conditional survival function of the event given confounders, and the conditional survival function of censoring given confounders. At the moment, we provide two ways of estimating each of these nuisance parameters. The treatment propensity can be estimated either using a logistic regression or SuperLearner (van der Laan and Polley, 2006). The conditional survival functions of event and censoring can be estimated using Cox proportional hazard models or survival random forests (Ishwaran et al.). Additionally, the estimates of any of these parameters can be passed in directly if another method is desired.
#'
#'  \code{CFsurvival} by default returns 95\% pointwise confidence intervals and 95\% uniform confidence bands, both of which are based on the influence function of the estimator. Contrasts of the counterfactual survival under treatment and control can also be obtained, including the survival difference, survival ratio, risk ratio, and number needed to treat. Pointwise confidence intervals and uniform confidence bands for these contrasts will also be computed.
#'
#' @param time \code{n x 1} numeric vector of observed right-censored follow-up times; i.e. the minimum of the event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of whether an event was observed.
#' @param treat \code{n x 1} numeric vector of observed binary treatment/exposure group; either 0 or 1.
#' @param confounders \code{n x p} numeric matrix of potential confounders to use when estimating the conditional survival probabilities. Missing values are not allowed
#' @param fit.times \code{k x 1} numeric vector of time grid at which the counterfactual survival function should be computed. Only values > 0 are allowed. Note that even if the survival is only desired at one timepoint (e.g. a final study time), fit.times should be a relatively fine grid between 0 and this final time in order for the computation to work. Defauls to all unique positive values in \code{time}.
#' @param fit.treat Optional subset of \code{c(0,1)} for which the counterfactual survival curves are desired. For example, if \code{fit.treat=c(0,1)} (default behavior) then both curves will be estimated, and if \code{fit.treat=0} then only the curve under non-treatment will be estimated.
#' @param nuisance.options List of options for nuisance parameter specification. See \code{\link{nuisance.options}} for details. Defaults to SuperLearners with rich libraries.
#' @param conf.band Logical indicating whether to compute simultaneous confidence bands.
#' @param conf.level Desired coverage of confidence intervals/bands.
#' @param contrasts Character vector indicating which (if any) contrasts of the survival functions are desired. Can include \code{"surv.diff"} (difference of survival functions), \code{"surv.ratio"} (ratio of survival functions), \code{"risk.ratio"} (ratio of the risk functions), or \code{"nnt"} (number needed to treat; i.e. recipricol of survival difference). If not \code{NULL}, then \code{fit.treat} must be \code{c(0,1)}. Defaults to \code{c("surv.diff", "surv.ratio")}.
#' @param verbose Logical indicating whether progress should be printed to the command line.
#' @return \code{CFsurvival} returns a named list with the following elements:
#' \item{fit.times}{The time points at which the counterfactual survival curves (and contrasts) were fit.}
#' \item{surv.df}{A data frame with the estimated counterfactual survival functions and CIs.}
#' \item{IF.vals.0, IF.vals.1}{n x k matrices with the influence values for the counterfactual survival functions. Rows index observations, columns index time points.}
#' \item{surv.diff.df}{A data frame with the estimated counterfactual survival difference (treatment survival minus control survival), as well as confidence intervals and tests of the null hypothesis that the difference equals zero. Only returned if \code{surv.diffs=TRUE}.}
#' \item{surv.ratio.df}{A data frame with the estimated counterfactual survival ratio (treatment survival divided by control survival), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{surv.ratios=TRUE}.}
#' \item{risk.ratio.df}{A data frame with the estimated counterfactual risk ratio (treatment risk divided by control risk), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{risk.ratios=TRUE}.}
#' \item{nnt.df}{A data frame with the estimated counterfactual number needed to treat (one divided by the survival difference), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{nnt=TRUE}.}
#' \item{surv.0.sim.maxes, surv.1.sim.maxes, surv.diff.}{Samples from the approximate distribution of the maximum of the limiting Gaussian processes of the counterfactual survival functions.}
#' \item{prop.pred}{The estimated treatment propensities.}
#' \item{event.pred.0, event.pred.1}{The estimated conditional survival functions of the event.}
#' \item{cens.pred.0, censd.pred.1}{The estimated conditional survival functions of censoring.}
#' \item{data}{The original time, event, and treatment data supplied to the function.}
#' @details \code{surv.df$ptwise.lower} and \code{surv.df$ptwise.upper} are lower and upper endpoints of back-transformed pointwise confidence intervals on the logit scale. Pointwise confidence intervals are only provided for values of \code{fit.times} that are at least as large as the smallest observed event time in the corresponding treatment cohort. \code{surv.df$unif.ew.lower} and \code{surv.df$unif.ew.upper} are lower and upper endpoints of an equal-width confidence band. \code{surv.df$unif.logit.lower} and \code{surv.df$unif.logit.upper} are lower and upper endpoints of an equi-precision confidence band on the logit scale, and trasformed back to the survival scale. This confidence band only covers values of \code{fit.times} such that the corresponding estimated survival function is \code{<= .99} and \code{>= .01}.
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
#' fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx, confounders =  data.frame(covar), contrasts = NULL, verbose = TRUE, fit.times = seq(0, 14.5, by=.1), fit.treat = c(0,1),
#' nuisance.options = list(prop.SL.library = c("SL.mean", "SL.bayesglm"),
#'                        event.SL.library = c("survSL.km", "survSL.coxph", "survSL.weibreg", "survSL.expreg"),
#'                        cens.SL.library = c("survSL.km", "survSL.coxph", "survSL.weibreg", "survSL.expreg"), #' cross.fit =TRUE, V = 5, save.nuis.fits = FALSE), conf.band = TRUE, conf.level = .95)

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
#'     geom_step(aes(time, surv, color=as.factor(trt), group=trt)) +
#'     geom_step(aes(time, ptwise.lower, color=as.factor(trt), group=trt), linetype=2) +
#'     geom_step(aes(time, ptwise.upper, color=as.factor(trt), group=trt), linetype=2) +
#'     geom_step(aes(time, unif.lower, color=as.factor(trt), group=trt), linetype=3) +
#'     geom_step(aes(time, unif.upper, color=as.factor(trt), group=trt), linetype=3) +
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
#'     coord_cartesian(xlim=c(0,15))
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
#'     geom_hline(yintercept=0, color='blue', linetype=4) +
#'     coord_cartesian(xlim=c(0,15), ylim=c(1,10))
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


CFsurvival <- function(time, event, treat, confounders, fit.times=sort(unique(time[time > 0 & time < max(time[event == 1])])), fit.treat=c(0,1), nuisance.options = list(), conf.band=TRUE, conf.level=.95, contrasts = c("surv.diff", "surv.ratio"), verbose=FALSE) {
    .args <- mget(names(formals()), sys.frame(sys.nframe()))

    # if(is.null(contrasts)) contrasts <- c("surv.diff", "surv.ratio")

    nuis <- do.call("CFsurvival.nuisance.options", nuisance.options)

    n <- length(time)

    if(sum(1-event) == 0) {
        message("No censored events; unanticipated errors may occur.")
        nuis$cens.pred.0 <- nuis$cens.pred.1 <- matrix(1, nrow=n, ncol=k)
    }
    if(sum(event) == 0) {
        stop("No uncensored events; cannot perform estimation.")
    }

    if(!nuis$cross.fit) {
        if(!is.null(nuis$folds) | nuis$V > 1) {
            message("nuisance.options$cross.fit set to FALSE, but V > 1 or folds specified. Cross-fitting will not be performed.")
        }
        nuis$V <- 1
        nuis$folds <- rep(1, n)
    }
    if(!is.null(nuis$folds)) {
        nuis$folds <- as.numeric(factor(nuis$folds))
        nuis$V <- length(unique(nuis$folds))
    }
    if(nuis$cross.fit & is.null(nuis$folds)) {
        while (sum(event) < nuis$V | sum(1-event) < nuis$V) {
            if (nuis$V == 1) {
                message("Cannot perform cross-fitting with one uncensored event.")
                nuis$cross.fit <- FALSE
                nuis$V <- 1
                break
            }
            message(paste0("Number of (event = 1) or (event = 0) is less than the number of folds; reducing number of folds to ", nuis$V-1))
            nuis$V <- nuis$V-1
        }
        event.0 <- which(event == 0)
        event.1 <- which(event == 1)
        folds.0 <- sample(rep(1:nuis$V, length = length(event.0)))
        folds.1 <- sample(rep(1:nuis$V, length = length(event.1)))
        nuis$folds <- rep(NA, n)
        nuis$folds[event.0] <- folds.0
        nuis$folds[event.1] <- folds.1
    }


    contrasts <- unique(tolower(contrasts))
    contrasts <- contrasts[contrasts %in% c("surv.diff", "surv.ratio", "risk.ratio", "nnt")]
    .check.input(time=time, event=event, treat=treat, confounders=confounders, fit.times=fit.times, fit.treat=fit.treat, nuisance.options=nuis, conf.band=conf.band, conf.level=conf.level, contrasts=contrasts, verbose=verbose)

    if(any(fit.times < 0)) {
        fit.times <- fit.times[fit.times > 0]
        message("fit.times < 0 removed.")
    }
    if(any(fit.times == 0)) fit.times <- fit.times[fit.times > 0]
    if(any(fit.times > max(time[event == 1]))) {
        fit.times <- fit.times[fit.times <= max(time[event == 1])]
        message("fit.times > max(time[event == 1]) removed.")
    }


    if(is.null(nuis$eval.times)) nuis$eval.times <- sort(unique(c(0,time[time > 0 & time <= max(fit.times)], max(fit.times))))
    k <- length(nuis$eval.times)

    confounders <- as.data.frame(confounders)

    surv.df <- data.frame()
    result <- list(fit.times=fit.times, fit.treat=fit.treat, surv.df=surv.df)

    #### ESTIMATE PROPENSITY ####

    if(is.null(nuis$prop.pred)) {
        if(verbose) message("Estimating propensities...")
        nuis$prop.pred <- rep(NA, n)
        if(nuis$V > 1) {
            if(nuis$save.nuis.fits) result$prop.fits <- vector(mode='list',length=nuis$V)
            for(v in 1:nuis$V) {
                if(verbose) message(paste("Fold ", v, "..."))
                train <- nuis$fold != v
                test <- nuis$fold == v
                prop.fit <- .estimate.propensity(A=treat[train], W=confounders[train,,drop=FALSE], newW=confounders[test,, drop=FALSE], SL.library=nuis$prop.SL.library, save.fit = nuis$save.nuis.fits, verbose = FALSE)
                nuis$prop.pred[test] <- prop.fit$prop.pred
                if(nuis$save.nuis.fits) result$prop.fits[[v]] <- prop.fit$prop.fit
            }
        } else {
            prop.fit <- .estimate.propensity(A=treat, W.propensity=confounders[,nuis$prop.subset, drop=FALSE], newW=confounders[,nuis$prop.subset, drop=FALSE], SL.library=nuis$prop.SL.library, save.fit = nuis$save.nuis.fits, verbose = FALSE)
            nuis$prop.pred[test] <- prop.fit$prop.pred
            if(nuis$save.nuis.fits) {
                result <- c(result, prop.fit = prop.fit$prop.fit)
            }
        }
    }

    #### ESTIMATE CONDITIONAL SURVIVALS ####

    if((1 %in% fit.treat & (is.null(nuis$event.pred.1) | is.null(nuis$cens.pred.1))) | (0 %in% fit.treat & (is.null(nuis$event.pred.0) | is.null(nuis$cens.pred.0)))) {
        if(verbose) message("Estimating conditional survivals...")
        if(0 %in% fit.treat & is.null(nuis$event.pred.0)) {
            do.event.pred.0 <- TRUE
            nuis$event.pred.0 <- matrix(NA, nrow=n, ncol=k)
        } else do.event.pred.0 <- FALSE
        if(1 %in% fit.treat & is.null(nuis$event.pred.1)) {
            do.event.pred.1 <- TRUE
            nuis$event.pred.1 <- matrix(NA, nrow=n, ncol=k)
        } else do.event.pred.1 <- FALSE
        if(0 %in% fit.treat & is.null(nuis$cens.pred.0)) {
            do.cens.pred.0 <- TRUE
            nuis$cens.pred.0 <- matrix(NA, nrow=n, ncol=k)
        } else do.cens.pred.0 <- FALSE
        if(1 %in% fit.treat & is.null(nuis$cens.pred.1)) {
            do.cens.pred.1 <- TRUE
            nuis$cens.pred.1 <- matrix(NA, nrow=n, ncol=k)
        } else do.cens.pred.1 <- FALSE
        if(nuis$V > 1) {
            if(nuis$save.nuis.fits) result$surv.fits <- vector(mode='list',length=nuis$V)
            for(v in 1:nuis$V) {
                if(verbose) message(paste("Fold ", v, "..."))
                train <- nuis$fold != v
                test <- nuis$fold == v
                surv.fit <- .estimate.conditional.survival(Y=time[train], Delta=event[train], A=treat[train], W=confounders[train,, drop=FALSE], newW=confounders[test,, drop=FALSE], event.SL.library=nuis$event.SL.library, fit.times=nuis$eval.times, fit.treat=fit.treat, cens.SL.library=nuis$cens.SL.library, save.fit = nuis$save.nuis.fits, verbose = FALSE)
                if(do.event.pred.0) nuis$event.pred.0[test,] <- surv.fit$event.pred.0
                if(do.event.pred.1) nuis$event.pred.1[test,] <- surv.fit$event.pred.1
                if(do.cens.pred.0) nuis$cens.pred.0[test,] <- surv.fit$cens.pred.0
                if(do.cens.pred.1) nuis$cens.pred.1[test,] <- surv.fit$cens.pred.1
                if(nuis$save.nuis.fits) result$surv.fits[[v]] <- surv.fit$surv.fit
            }
        } else {
            surv.fit <- .estimate.conditional.survival(Y=time, Delta=event, A=treat, W=confounders, newW=confounders, event.SL.library=nuis$event.SL.library, fit.times=nuis$eval.times, fit.treat=fit.treat, cens.SL.library=nuis$cens.SL.library, save.fit = nuis$save.nuis.fits, verbose = FALSE)
            if(do.event.pred.0) nuis$event.pred.0 <- surv.fit$event.pred.0
            if(do.event.pred.1) nuis$event.pred.1 <- surv.fit$event.pred.1
            if(do.cens.pred.0) nuis$cens.pred.0 <- surv.fit$cens.pred.0
            if(do.cens.pred.1) nuis$cens.pred.1 <- surv.fit$cens.pred.1
            if(nuis$save.nuis.fits) result$surv.fit <- surv.fit$surv.fit
        }
    }

    #### ESTIMATE CF SURVIVALS ####
    if(verbose) message("Computing counterfactual survivals...")
    if(0 %in% fit.treat) {
        surv.0 <- .get.survival(Y=time, Delta=event, A=1-treat, fit.times=fit.times, eval.times=nuis$eval.times, S.hats=nuis$event.pred.0, G.hats=nuis$cens.pred.0, g.hats=1-nuis$prop.pred)
        surv.df.0 <- data.frame(time=c(0,fit.times), trt=0, surv=c(1, surv.0$surv.iso))
        result$IF.vals.0 <- surv.0$IF.vals

        q.01 <- quantile(time[event == 1 & treat == 0], .01)#min(fit.times[surv.0$surv <= .99])
        q.99 <- max(fit.times[!is.na(surv.0$surv.iso) && surv.0$surv.iso >= .01])
        min.obs.time <- min(time[event == 1 & treat == 0])

        c.int <- .surv.confints(fit.times, surv.0$surv, surv.0$IF.vals, conf.band = conf.band, band.end.pts = c(q.01, q.99), conf.level=conf.level)
        surv.df.0$se <- c(0,c.int$res$se)
        surv.df.0$se.logit <- c(0,c.int$res$se.logit)
        surv.df.0$ptwise.lower <- c(1,c.int$res$ptwise.lower)
        surv.df.0$ptwise.upper <- c(1,c.int$res$ptwise.upper)
        surv.df.0$ptwise.logit.lower <- c(1,c.int$res$ptwise.logit.lower)
        surv.df.0$ptwise.logit.upper <- c(1,c.int$res$ptwise.logit.upper)

        #surv.df.0$se[surv.df.0$time < min.obs.time & surv.df.0$time > 0] <- surv.df.0$se.logit[surv.df.0$time < min.obs.time & surv.df.0$time > 0] <- surv.df.0$ptwise.lower[surv.df.0$time < min.obs.time & surv.df.0$time > 0] <- surv.df.0$ptwise.upper[surv.df.0$time < min.obs.time & surv.df.0$time > 0] <- NA

        if(conf.band) {
            surv.df.0$unif.ew.lower <- c(1,c.int$res$unif.ew.lower)
            surv.df.0$unif.ew.upper <- c(1,c.int$res$unif.ew.upper)
            result$surv.0.unif.ew.quant <- c.int$unif.ew.quant
            result$surv.0.ew.sim.maxes <- c.int$ew.sim.maxes

            surv.df.0$unif.logit.lower <- c(1,c.int$res$unif.logit.lower)
            surv.df.0$unif.logit.upper <- c(1,c.int$res$unif.logit.upper)
            result$surv.0.unif.logit.quant <- c.int$unif.logit.quant
            result$surv.0.logit.sim.maxes <- c.int$logit.sim.maxes
        }
        result$surv.df <- rbind(result$surv.df, surv.df.0)
    }

    if(1 %in% fit.treat) {
        surv.1 <- .get.survival(Y=time, Delta=event, A=treat, fit.times=fit.times, eval.times=nuis$eval.times, S.hats=nuis$event.pred.1, G.hats=nuis$cens.pred.1, g.hats=nuis$prop.pred)
        surv.df.1 <- data.frame(time=c(0,fit.times), trt=1, surv=c(1, surv.1$surv.iso))
        result$IF.vals.1 <- surv.1$IF.vals

        q.01 <- quantile(time[event == 1 & treat == 1], .01)#min(fit.times[surv.1$surv <= .99])
        q.99 <- max(fit.times[surv.1$surv.iso >= .01])
        min.obs.time <- min(time[event == 1 & treat == 1])

        c.int <- .surv.confints(fit.times, surv.1$surv, surv.1$IF.vals, conf.band = conf.band, band.end.pts = c(q.01, q.99), conf.level=conf.level)
        surv.df.1$se <- c(0,c.int$res$se)
        surv.df.1$se.logit <- c(0,c.int$res$se.logit)
        surv.df.1$ptwise.lower <- c(1,c.int$res$ptwise.lower)
        surv.df.1$ptwise.upper <- c(1,c.int$res$ptwise.upper)
        surv.df.1$ptwise.logit.lower <- c(1,c.int$res$ptwise.logit.lower)
        surv.df.1$ptwise.logit.upper <- c(1,c.int$res$ptwise.logit.upper)

        #surv.df.1$se[surv.df.1$time < min.obs.time & surv.df.1$time > 0] <- surv.df.1$se.logit[surv.df.1$time < min.obs.time & surv.df.1$time > 0] <- surv.df.1$ptwise.lower[surv.df.1$time < min.obs.time & surv.df.1$time > 0] <- surv.df.1$ptwise.upper[surv.df.1$time < min.obs.time & surv.df.1$time > 0] <- NA

        if(conf.band) {
            surv.df.1$unif.ew.lower <- c(1,c.int$res$unif.ew.lower)
            surv.df.1$unif.ew.upper <- c(1,c.int$res$unif.ew.upper)
            result$surv.1.unif.ew.quant <- c.int$unif.ew.quant
            result$surv.1.ew.sim.maxes <- c.int$ew.sim.maxes

            surv.df.1$unif.logit.lower <- c(1,c.int$res$unif.logit.lower)
            surv.df.1$unif.logit.upper <- c(1,c.int$res$unif.logit.upper)
            result$surv.1.unif.logit.quant <- c.int$unif.logit.quant
            result$surv.1.logit.sim.maxes <- c.int$logit.sim.maxes
        }
        result$surv.df <- rbind(result$surv.df, surv.df.1)
    }

    #### ESTIMATE CF CONTRASTS ####
    if(identical(sort(fit.treat), c(0,1))) {
        q.01 <- quantile(time[event == 1], .01)#min(fit.times[surv.0$surv <= .99])
        q.99 <- max(fit.times[surv.0$surv.iso >= .01 | surv.1$surv.iso >= .01])
    }

    if('surv.diff' %in% contrasts & identical(sort(fit.treat), c(0,1))) {
        if(verbose) message("Computing survival differences...")
        surv.diff <- .surv.difference(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, band.end.pts=c(q.01, q.99), conf.level=conf.level)
        result <- c(result, surv.diff)
    }

    if('surv.ratio' %in% contrasts & identical(sort(fit.treat), c(0,1))) {
        if(verbose) message("Computing survival ratios...")
        surv.ratio <- .surv.ratio(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, band.end.pts=c(q.01, q.99),conf.level=conf.level)
        result <- c(result, surv.ratio)
    }

    if('risk.ratio' %in% contrasts  & identical(sort(fit.treat), c(0,1))) {
        if(verbose) message("Computing risk ratios...")
        risk.ratio <- .risk.ratio(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, band.end.pts=c(q.01, q.99), conf.level=conf.level)
        result <- c(result, risk.ratio)
    }

    if('nnt' %in% contrasts  & identical(sort(fit.treat), c(0,1))) {
        if(verbose) message("Computing number needed to treat...")
        nnt <- .nnt(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, band.end.pts=c(q.01, q.99), conf.level=conf.level)
        result <- c(result, nnt)
    }

    result$data <- data.frame(time, event, treat)
    result$nuisance <- nuis

    return(result)
}


#' Initialize options for nuisance estimation
#'
#' This function initializes the options for nuisance parameter estimation (i.e. conditional survivals of event and censoring and treatment propensity) for use in \code{CFsurvival}. The conditional survivals of event and censoring are estimated using \code{\link[survSuperLearner]{survSuperLearner}}, or any of the individual learners therein. Treatment propensity is estimated using the \code{\link[SuperLearner]{SuperLearner}}, or any of the individual learners therein. Alternatively, the estimation process can be overriden by providing predictions from pre-fit nuisance estimators.
#'
#' @param cross.fit Logical indicating whether to cross-fit nuisance parameters. Defaults to \code{TRUE}.
#' @param V Positive integer number of folds for cross-fitting. Defaults to 10.
#' @param folds Optional \code{n x 1} vector indicating which fold each observation is in. If \code{NULL}, folds will be randomly assigned in such a way to balance observed events across folds.
#' @param eval.times Grid of time values on which to perform the estimation procedure. Defaults to the sorted unique values in \code{time}.
#' @param save.nuis.fits Logical indicating whether to save the fitted nuisance objects.
#' @param event.SL.library The library of candidate learners to be passed on to \code{\link[survSuperLearner]{survSuperLearner}} to estimate the conditional survival of the event. If only a single learner is provided(e.g. \code{survSL.km} for Kaplan-Meier estimator, \code{survSL.coxph} for Cox model, or \code{survSL.rfsrc} for survival random forest), then just this learner will be used, and no super learning will be performed. Defaults to a full library with screening and all algorithms currently implemented in \code{survSuperLearner}. If \code{event.SL.library = NULL}, then \code{event.pred.0} and/or \code{event.pred.1} must be specified.
#' @param event.pred.0 Optional \code{n x k} matrix of estimates of the conditional survival of the event given treatment = 0 and confounders. If \code{event.SL.library = NULL} and \code{0 \%in\% fit.treat}, then \code{event.pred.0} must be specified. If \code{event.SL.library} is not \code{NULL}, then \code{event.pred.0 } is ignored.
#' @param event.pred.1 Optional \code{n x k} matrix of estimates of the conditional survival of the event given treatment = 1 and confounders. If \code{event.SL.library = NULL} and \code{1 \%in\% fit.treat}, then \code{event.pred.1} must be specified. If \code{event.SL.library} is not \code{NULL}, then \code{event.pred.1 } is ignored.
#' @param cens.SL.library The library of candidate learners to estimate the conditional survival of censoring. As with \code{event.SL.library}, single learners can be specified. Defaults to a full library with screening and all algorithms currently implemented in \code{survSuperLearner}. If \code{cens.SL.library = NULL}, then \code{cens.pred.0} and/or \code{cens.pred.1} must be specified.
#' @param cens.pred.0 Optional \code{n x k} matrix of estimates of the conditional survival of censoring given treatment = 0 and confounders. If \code{cens.SL.library = NULL} and \code{0 \%in\% fit.treat}, then \code{cens.pred.0} must be specified. If \code{cens.SL.library} is not \code{NULL}, then \code{cens.pred.0 } is ignored.
#' @param cens.pred.1 Optional \code{n x k} matrix of estimates of the conditional survival of censoring given treatment = 1 and confounders. If \code{cens.SL.library = NULL} and \code{1 \%in\% fit.treat}, then \code{cens.pred.1} must be specified. If \code{cens.SL.library} is not \code{NULL}, then \code{cens.pred.1} is ignored.
#' @param prop.SL.library The library to use for estimation of the treatment propensities using \code{\link[SuperLearner]{SuperLearner}}. If only a single learner is provided(e.g. \code{SL.mean} for marginal mean, \code{SL.glm} for logistic regression, or \code{SL.ranger} for random forest), then just this learner will be used, and no super learning will be performed. If \code{prop.SL.library = NULL}, then \code{prop.pred} must be provided.
#' @param prop.pred Optional \code{n x 1} numeric vector of estimated probabilities that \code{treat = 1} given the confounders. If \code{prop.SL.library = NULL}, then \code{prop.pred} must be specified, otherwise it is ignored.
#' @return Named list containing the nuisance options.
CFsurvival.nuisance.options <- function(cross.fit = TRUE, V = 10, folds = NULL, eval.times = NULL, save.nuis.fits = FALSE,
                                        event.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "surv.screen.glmnet", "surv.screen.marg", "All") ),  event.pred.0 = NULL, event.pred.1 = NULL,
                                        cens.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "surv.screen.glmnet", "surv.screen.marg", "All") ),  cens.pred.0 = NULL, cens.pred.1 = NULL,
                                        prop.SL.library = lapply(c("SL.mean", "SL.bayesglm", "SL.gam", "SL.earth", "SL.ranger", "SL.xgboost"), function(alg) c(alg, "screen.glmnet", "screen.corRank", "All") ), prop.pred = NULL) {
    list(cross.fit = cross.fit, V = V, folds = folds, eval.times = eval.times, save.nuis.fits = save.nuis.fits,
         event.SL.library = event.SL.library, event.pred.0 = event.pred.0, event.pred.1 = event.pred.1,
         cens.SL.library = cens.SL.library, cens.pred.0 = cens.pred.0, cens.pred.1 = cens.pred.1,
         prop.SL.library = prop.SL.library,  prop.pred = prop.pred)
}

.check.input <- function(time, event, treat, confounders, fit.times, fit.treat, nuisance.options, conf.band, conf.level, contrasts, verbose, ...) {
    if(any(time < 0)) stop("Only non-negative event/censoring times allowed!")
    if(any(time == 0 & event == 1)) stop("Events at time zero not allowed.")
    if(any(!(event %in% c(0,1)))) stop("Event must be binary.")
    if(any(!(treat %in% c(0,1)))) stop("Treatment must be binary.")
    if(length(time) != length(event) | length(time) != length(treat)) stop("time, event, and treat must be n x 1 vectors")
    if(any(!(fit.treat %in% c(0,1)))) stop("fit.treat must be a subset of c(0,1).")
    if(!is.null(nuisance.options$event.pred.1) && dim(nuisance.options$event.pred.1) != c(nrow(time), length(fit.times))) {
        stop("event.pred must be an n x k matrix (n is number of observations, k is length of fit.times)")
    }
    if(!is.null(nuisance.options$event.pred.0) && dim(nuisance.options$event.pred.0) != c(nrow(time), length(fit.times))) {
        stop("event.pred must be an n x k matrix (n is number of observations, k is length of fit.times)")
    }
    if(!is.null(nuisance.options$cens.pred.1) && dim(nuisance.options$cens.pred.1) != c(nrow(time), length(fit.times))) {
        stop("cens.pred must be an n x k matrix (n is number of observations, k is length of fit.times)")
    }
    if(!is.null(nuisance.options$cens.pred.0) && dim(nuisance.options$cens.pred.0) != c(nrow(time), length(fit.times))) {
        stop("cens.pred must be an n x k matrix (n is number of observations, k is length of fit.times)")
    }
    if(is.null(nuisance.options$event.SL.library)) {
        if(0 %in% fit.treat & is.null(nuisance.options$event.pred.0)) {
            stop("event.pred.0 must be provided if event.SL.library is not specified and 0 is a treatment of interest.")
        }
        if(1 %in% fit.treat & is.null(nuisance.options$event.pred.1)) {
            stop("event.pred.1 must be provided if event.SL.library is not specified and 1 is a treatment of interest.")
        }
    }
    if(is.null(nuisance.options$cens.SL.library)) {
        if(0 %in% fit.treat & is.null(nuisance.options$cens.pred.0)) {
            stop("cens.pred.0 must be provided if cens.SL.library is not specified and 0 is a treatment of interest.")
        }
        if(1 %in% fit.treat & is.null(nuisance.options$cens.pred.1)) {
            stop("cens.pred.1 must be provided if cens.SL.library is not specified and 1 is a treatment of interest.")
        }
    }
    if(length(contrasts) > 0 & !identical(sort(fit.treat), c(0,1))) {
        warning("contrast specified but both treatment regimens not requested -- contrasts will not be provided. Re-run with fit.treat = c(0,1) for survival contrasts.")
    }
    if(any(is.na(time) | is.na(event))) {
        stop("Missing time or event detected; missing data not allowed.")
    }
    if(any(is.na(confounders) | is.na(treat))) {
        warning("Missing confounders or exposures detected; missing data not allowed.")
    }
}

.get.survival <- function(Y, Delta, A, fit.times, eval.times, S.hats, G.hats, g.hats, isotonize=TRUE) {
    fit.times <- fit.times[fit.times > 0]
    n <- length(Y)
    ord <- order(eval.times)
    eval.times <- eval.times[ord]
    S.hats <- S.hats[,ord]
    G.hats <- G.hats[,ord]

    int.vals <- t(sapply(1:n, function(i) {
        vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
        if(any(eval.times[-1] > Y[i])) vals[eval.times[-1] > Y[i]] <- 0
        c(0,cumsum(vals))
    }))
    S.hats.Y <- sapply(1:n, function(i) stepfun(eval.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
    G.hats.Y <- sapply(1:n, function(i) stepfun(eval.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
    IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
    surv <- rep(NA, length(fit.times))
    for(t0 in fit.times) {
        k <- min(which(eval.times >= t0))
        S.hats.t0 <- S.hats[,k]
        inner.func.1 <- ifelse(Y <= t0 & Delta == 1, 1/(S.hats.Y * G.hats.Y), 0 )
        inner.func.2 <- int.vals[,k]
        if.func <- as.numeric(A == 1) * S.hats.t0 * ( -inner.func.1 + inner.func.2) / g.hats + S.hats.t0
        k1 <- which(fit.times == t0)
        surv[k1] <- mean(if.func)
        IF.vals[,k1] <- if.func - surv[k1]
    }
    res <- list(times=fit.times, surv=pmin(1,pmax(0,surv)), IF.vals=IF.vals)
    if(isotonize) {
        res$surv.iso <- NA
        res$surv.iso[!is.na(res$surv)] <- 1 - isoreg(res$times[!is.na(res$surv)], 1-res$surv[!is.na(res$surv)])$yf
    }

    return(res)
}

.surv.confints <- function(times, est, IF.vals, isotonize=TRUE, conf.band=TRUE, band.end.pts=c(0,Inf), conf.level=.95) {
    logit <- function(x) log(x / (1-x))
    logit.prime <- function(x) 1/(x * (1-x))
    expit <- function(x) 1/(1 + exp(-x))

    n <- nrow(IF.vals)
    in.bounds <- est > 0 & est < 1
    IF.vals.logit <- IF.vals
    for(j in 1:length(est)) IF.vals.logit[,j] <- IF.vals[,j] * logit.prime(est[j])

    res <- NULL
    res$se <- sqrt(colMeans(IF.vals^2)) / sqrt(n)
    res$se.logit <- sqrt(colMeans(IF.vals.logit^2)) / sqrt(n)

    res$se[res$se == 0] <- NA
    res$se.logit[is.infinite(res$se.logit)] <- NA
    quant <- qt(1-(1-conf.level)/2, n - 1) #qnorm(1-(1-conf.level)/2)

    #Intervals and bands based on un-isotonized survival ests
    res$ptwise.lower <- pmax(est - quant * res$se, 0)
    res$ptwise.upper <- pmin(est + quant * res$se, 1)
    res$ptwise.logit.lower <- expit(logit(est) - quant * res$se.logit) #pmax(est - quant * res$se, 0)
    res$ptwise.logit.upper <- expit(logit(est) + quant * res$se.logit) #pmin(est + quant * res$se, 1)
    #
    # F.vals <-  isoreg(times, 1-est)$yf
    # F.inv.hat <- stepfun(F.vals[-1], times)
    #
    # ln.q.ptwise <- as.numeric(F.inv.hat(expit(logit(F.vals) - quant * res$se.logit )))
    # un.q.ptwise <- as.numeric(F.inv.hat(expit(logit(F.vals) + quant * res$se.logit )))

    # res$ptwise.inv.upper <- 1-sapply(times, function(t) {
    #     if(any(un.q.ptwise <= t,na.rm=TRUE)) F.vals[max(which(un.q.ptwise <= t), na.rm=TRUE)]
    #     else 0
    # })
    # res$ptwise.inv.lower <- 1-sapply(times, function(t) {
    #     if(any(ln.q.ptwise >= t,na.rm=TRUE)) F.vals[min(which(ln.q.ptwise >= t),na.rm=TRUE)]
    #     else min(ln.q.ptwise, na.rm=TRUE)
    # })

    # Isotonized intervals and bands
    out <- NULL
    if(conf.band) {
        unif.vals <- .estimate.uniform.quantile(IF.vals[,!is.na(res$se)], conf.level, scale = FALSE)
        unif.quant <- unif.vals$quantile
        out$ew.sim.maxes <- unif.vals$maxes
        out$unif.ew.quant <- unif.quant
        res$unif.ew.lower <- pmax(est - unif.quant / sqrt(n), 0)
        res$unif.ew.upper <- pmin(est + unif.quant / sqrt(n), 1)

        if(isotonize) {
            res$unif.ew.lower[!is.na(res$unif.ew.lower)] <- 1 - isoreg(times[!is.na(res$unif.ew.lower)], 1-res$unif.ew.lower[!is.na(res$unif.ew.lower)])$yf
            res$unif.ew.upper[!is.na(res$unif.ew.upper)] <- 1 - isoreg(times[!is.na(res$unif.ew.upper)], 1-res$unif.ew.upper[!is.na(res$unif.ew.upper)])$yf
        }

        unif.logit.vals <- .estimate.uniform.quantile(IF.vals.logit[,!is.na(res$se.logit) & times >= band.end.pts[1] & times <= band.end.pts[2]], conf.level)
        unif.logit.quant <- unif.logit.vals$quantile
        out$logit.sim.maxes <- unif.logit.vals$maxes
        out$unif.logit.quant <- unif.logit.quant
        res$unif.logit.lower <- expit(logit(est) - unif.logit.quant * res$se.logit) #pmax(est - unif.quant * res$se, 0)
        res$unif.logit.upper <- expit(logit(est) + unif.logit.quant * res$se.logit) # pmin(est + unif.quant * res$se, 1)
        res$unif.logit.lower[times < band.end.pts[1] | times > band.end.pts[2]] <- NA
        res$unif.logit.upper[times < band.end.pts[1] | times > band.end.pts[2]] <- NA
        # res$unif.lower[!in.bounds] <- pmin(approx(times[in.bounds], res$unif.lower[in.bounds], xout=times[!in.bounds], rule=2)$y, est[!in.bounds])
        # res$unif.upper[!in.bounds] <- pmax(approx(times[in.bounds], res$unif.lower[in.bounds], xout=times[!in.bounds], rule=2)$y, est[!in.bounds])
        if(isotonize) {
            res$unif.logit.lower[!is.na(res$unif.logit.lower)] <- 1 - isoreg(times[!is.na(res$unif.logit.lower)], 1-res$unif.logit.lower[!is.na(res$unif.logit.lower)])$yf
            res$unif.logit.upper[!is.na(res$unif.logit.upper)] <- 1 - isoreg(times[!is.na(res$unif.logit.upper)], 1-res$unif.logit.upper[!is.na(res$unif.logit.upper)])$yf
        }
    }
    out$res <- res
    return(out)
}

.estimate.uniform.quantile <- function(IF.vals, conf.level=.95, scale=TRUE) {
    n <- nrow(IF.vals)
    if(scale) IF.vals <- scale(IF.vals)
    maxes <- replicate(1e4, max(abs(rbind(rt(n, df = n - 1)/sqrt(n)) %*% IF.vals)))
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


