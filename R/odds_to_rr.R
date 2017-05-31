#' @title Get relative risks estimates from logistic regressions or odds ratio values
#' @name odds_to_rr
#'
#' @description \code{odds_to_rr()} converts odds ratios from a logistic regression
#'                model (including mixed models) into relative risks; \code{or_to_rr()}
#'                converts a single odds ratio estimate into a relative risk estimate.
#'
#' @param fit A fitted binomial generalized linear (mixed) model with logit-link function
#'          (logistic (multilevel) regression model).
#' @param or Numeric, an odds ratio estimate.
#' @param p0 Numeric, proportion of the incidence in the outcome variable
#'             (base line risk).
#'
#' @return A data frame with relative risks and lower/upper confidence interval for
#'           the relative risks estimates; for \code{or_to_rr()}, the risk ratio
#'           estimate.
#'
#' @references Zhang J, Yu KF. 1998. What's the Relative Risk? A Method of Correcting the Odds Ratio in Cohort Studies of Common Outcomes. JAMA; 280(19): 1690-1. \doi{10.1001/jama.280.19.1690}
#'             \cr \cr
#'             Grant RL. 2014. Converting an odds ratio to a range of plausible relative risks for better communication of research findings. BMJ 348:f7450. \doi{10.1136/bmj.f7450}
#'
#' @details This function extracts the odds ratios (exponentiated model coefficients)
#'            from logistic regressions (fitted with \code{glm} or \code{glmer})
#'            and their related confidence intervals, and transforms these values
#'            into relative risks (and their related confidence intervals).
#'            \cr \cr
#'            The formula for transformation is based on Zhang and Yu (1998)
#'            and Grant (2014):
#'            \code{RR <- OR / (1 - P0 + (P0 * OR))}, where \code{OR} is the odds
#'            ratio and \code{P0} indicates the proportion of the incidence in
#'            the outcome variable.
#'
#' @examples
#' library(sjmisc)
#' library(lme4)
#' # create binary response
#' sleepstudy$Reaction.dicho <- dicho(sleepstudy$Reaction, dich.by = "median")
#' # fit model
#' fit <- glmer(Reaction.dicho ~ Days + (Days | Subject),
#'              data = sleepstudy, family = binomial("logit"))
#' # convert to relative risks
#' odds_to_rr(fit)
#'
#'
#' data(efc)
#' # create binary response
#' y <- ifelse(efc$neg_c_7 < median(na.omit(efc$neg_c_7)), 0, 1)
#' # create data frame for fitted model
#' mydf <- data.frame(y = as.factor(y),
#'                    sex = efc$c161sex,
#'                    dep = to_factor(efc$e42dep),
#'                    barthel = efc$barthtot,
#'                    education = to_factor(efc$c172code))
#' # fit model
#' fit <- glm(y ~., data = mydf, family = binomial(link = "logit"))
#' # convert to relative risks
#' odds_to_rr(fit)
#'
#' # replicate OR/RR for coefficient "sex" from above regression
#' or_to_rr(1.913887, .5516)
#'
#' @export
odds_to_rr <- function(fit) {
  # check model family
  fitinfo <- get_glm_family(fit)
  # no binomial model with logit-link?
  if (!fitinfo$is_bin && !fitinfo$is_logit)
    stop("`fit` must be a binomial model with logit-link (logistic regression).", call. = F)
  # get model estimates
  est <- exp(stats::coef(summary(fit))[, 1])
  # get confidence intervals
  if (is_merMod(fit))
    ci <- stats::confint(fit, method = "Wald", parm = "beta_")
  else
    ci <- stats::confint(fit)
  # bind to data frame
  or.dat <- data.frame(est, exp(ci))
  colnames(or.dat) <- c("OR", "lower.ci", "upper.ci")
  # get P0, i.e. the incidence ratio of the outcome for the
  # non-exposed group
  modfram <- stats::model.frame(fit)
  # make sure that outcome is 0/1-numeric, so we can simply
  # compute the mean to get the ratio
  P0 <- mean(sjlabelled::as_numeric(modfram[[1]], start.at = 0, keep.labels = F), na.rm = T)
  # compute relative risks for estimate and confidence intervals
  rr.dat <- or.dat / ((1 - P0) + (P0 * or.dat))
  colnames(rr.dat) <- c("RR", "lower.ci", "upper.ci")
  rr.dat
}


#' @rdname odds_to_rr
#' @export
or_to_rr <- function(or, p0) {
  or / (1 - p0 + (p0 * or))
}
