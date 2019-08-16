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
#' @param p0 Numeric, the risk of having a positive outcome in the control or
#'   unexposed group (reference group), i.e. the number of outcome or "successes"
#'   in the control divided by the total number of observations in the control
#'   group.
#'
#' @return A data frame with relative risks and lower/upper confidence interval for
#'           the relative risks estimates; for \code{or_to_rr()}, the risk ratio
#'           estimate.
#'
#' @references
#'   Grant RL. 2014. Converting an odds ratio to a range of plausible relative risks for better communication of research findings. BMJ 348:f7450. \doi{10.1136/bmj.f7450}
#'   \cr \cr
#'   Wang Z. 2013. Converting Odds Ratio to Relative Risk in Cohort Studies with Partial Data Information. J Stat Soft 2013;55. \doi{10.18637/jss.v055.i05}
#'   \cr \cr
#'   Zhang J, Yu KF. 1998. What's the Relative Risk? A Method of Correcting the Odds Ratio in Cohort Studies of Common Outcomes. JAMA; 280(19): 1690-1. \doi{10.1001/jama.280.19.1690}
#'
#'
#' @details This function extracts the odds ratios (exponentiated model coefficients)
#'            from logistic regressions (fitted with \code{glm} or \code{glmer})
#'            and their related confidence intervals, and transforms these values
#'            into relative risks (and their related confidence intervals).
#'            \cr \cr
#'            The formula for transformation is based on Zhang and Yu (1998),
#'            Wang (2013) and Grant (2014):
#'            \code{RR <- OR / (1 - P0 + (P0 * OR))}, where \code{OR} is the odds
#'            ratio and \code{P0} indicates the proportion of the incidence in
#'            the outcome variable for the control group (reference group).
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
#' mydf <- data.frame(
#'   y = as.factor(y),
#'   sex = to_factor(efc$c161sex),
#'   dep = to_factor(efc$e42dep),
#'   barthel = efc$barthtot,
#'   education = to_factor(efc$c172code)
#' )
#' # fit model
#' fit <- glm(y ~., data = mydf, family = binomial(link = "logit"))
#' # convert to relative risks
#' odds_to_rr(fit)
#'
#' # replicate OR/RR for coefficient "sex" from above regression
#' # p0 ~ .44, or ~ 1.914
#' prop.table(table(mydf$y, mydf$sex))
#' or_to_rr(1.914, 0.1055 / (.1324 + .1055))
#'
#' @importFrom stats coef confint model.frame
#' @importFrom sjlabelled as_numeric
#' @export
odds_to_rr <- function(fit) {
  # check model family
  fitinfo <- get_glm_family(fit)

  # no binomial model with logit-link?
  if (!fitinfo$is_bin && !fitinfo$is_logit)
    stop("`fit` must be a binomial model with logit-link (logistic regression).", call. = F)

  # get model estimates
  est <- insight::get_parameters(fit)
  est$estimate <- exp(est$estimate)

  # get confidence intervals
  if (is_merMod(fit))
    ci <- stats::confint(fit, method = "Wald", parm = "beta_")
  else
    ci <- stats::confint(fit)

  # bind to data frame
  or.dat <- data.frame(est, exp(ci))
  colnames(or.dat) <- c("Parameter", "OR", "CI_low", "CI_high")

  # get P0, i.e. the incidence ratio of the outcome for the
  # non-exposed group
  modfram <- insight::get_data(fit)

  # make sure that outcome is 0/1-numeric, so we can simply
  # compute the mean to get the ratio
  outcome <- sjmisc::recode_to(sjlabelled::as_numeric(insight::get_response(fit)))

  P0 <- c()
  for (i in 1:nrow(est)) {
    P0 <- c(P0, .baseline_risk_for_predictor(modfram, outcome, est$parameter[i]))
  }

  # compute relative risks for estimate and confidence intervals
  rr.dat <- or.dat[, 2:4] / ((1 - P0) + (P0 * or.dat[, 2:4]))
  rr.dat <- cbind(or.dat$Parameter, or.dat$OR, rr.dat)

  colnames(rr.dat) <- c("Parameter", "Odds Ratio", "Risk Ratio", "CI_low", "CI_high")
  rownames(rr.dat) <- NULL

  rr.dat
}


.baseline_risk_for_predictor <- function(data, outcome, parameter) {
  if (parameter == "(Intercept)") return(mean(outcome))

  if (!(parameter %in% colnames(data))) {
    find.factors <- lapply(colnames(data), function(.i) {
      v <- data[[.i]]
      if (is.factor(v)) {
        return(paste0(.i, levels(v)))
      }
      return(.i)
    })
    names(find.factors) <- colnames(data)
    parameter <- names(find.factors)[which(sapply(find.factors, function(.i) {
      parameter %in% .i
    }))]
  }

  if (is.numeric(data[[parameter]])) {
    mean(outcome)
  } else {
    p <- prop.table(table(data[[parameter]], outcome))
    p[1, 2] / sum(p[1, ])
  }
}


#' @rdname odds_to_rr
#' @export
or_to_rr <- function(or, p0) {
  or / (1 - p0 + (p0 * or))
}
