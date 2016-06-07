#' @title Get relative risks estimates from logistic regressions
#' @name odds_to_rr
#'
#' @description This function converts odds ratios from a logistic regression
#'                model (including mixed models) into relative risks.
#'
#' @param fit A fitted binomial generalized linear (mixed) model with logit-link function
#'          (logistic (multilevel) regression model).
#'
#' @return A data frame with relative risks and lower/upper confidence interval for
#'           the relative risks estimates.
#'
#' @references Zhang J, Yu KF. 1998. What's the Relative Risk? A Method of Correcting the Odds Ratio in Cohort Studies of Common Outcomes. JAMA; 280(19): 1690-1. \doi{10.1001/jama.280.19.1690}
#'
#' @details This function extracts the odds ratios (exponentiated model coefficients)
#'            from logistic regressions (fitted with \code{glm} or \code{glmer})
#'            and their related confidence intervals, and transforms these values
#'            into relative risks (and their related confidence intervals).
#'            \cr \cr
#'            The formula for transformation is based on Zhang and Yu (1998):
#'            \code{RR <- OR / ((1 - P0) + (P0 * OR))}, where \code{OR} is the odds
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
  P0 <- mean(sjmisc::to_value(modfram[[1]], start.at = 0, keep.labels = F), na.rm = T)
  # compute relative risks for estimate and confidence intervals
  rr.dat <- or.dat / ((1 - P0) + (P0 * or.dat))
  colnames(rr.dat) <- c("RR", "lower.ci", "upper.ci")
  rr.dat
}
