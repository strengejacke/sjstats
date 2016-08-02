#' @title Tjur's Coefficient of Discrimination
#' @name cod
#'
#' @description This method calculates the Coefficient of Discrimination \code{D}
#'                for generalized linear (mixed) models for binary data. It is
#'                an alternative to other Pseudo-R-squared values
#'                like Nakelkerke's R2 or Cox-Snell R2.
#'
#' @param x Fitted \code{\link{glm}} or \code{\link[lme4]{glmer}} model.
#'
#' @return The \code{D} Coefficient of Discrimination, also known as
#'           Tjur's R-squared value.
#'
#' @note The Coefficient of Discrimination \code{D} can be read like any
#'         other (Pseudo-)R-squared value.
#'
#' @references Tjur T (2009) Coefficients of determination in logistic regression models -
#'               a new proposal: The coefficient of discrimination. The American Statistician,
#'               63(4): 366-372
#'
#' @seealso \code{\link{r2}} for Nagelkerke's and Cox and Snell's pseudo
#'            r-squared coefficients.
#'
#' @examples
#' library(sjmisc)
#' data(efc)
#'
#' # Tjur's R-squared value
#' efc$services <- ifelse(efc$tot_sc_e > 0, 1, 0)
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep,
#'            data = efc, family = binomial(link = "logit"))
#' cod(fit)
#'
#' @importFrom stats predict predict.glm residuals
#' @export
cod <- function(x) {
  # check for valid object class
  if (!any(class(x) == "glmerMod") && !any(class(x) == "glm")) {
    stop("`x` must be an object of class `glm` or `glmerMod`.", call. = F)
  }

  # mixed models (lme4)
  if (any(class(x) == "glmerMod")) {
    # check for package availability
    y <- lme4::getME(x, "y")
    pred <- stats::predict(x, type = "response", re.form = NULL)
  } else {
    y <- x$y
    pred <- stats::predict.glm(x, type = "response")
  }
  # delete pred for cases with missing residuals
  if (anyNA(stats::residuals(x))) pred <- pred[!is.na(stats::residuals(x))]

  categories <- unique(y)
  m1 <- mean(pred[which(y == categories[1])], na.rm = T)
  m2 <- mean(pred[which(y == categories[2])], na.rm = T)

  cod = abs(m2 - m1)
  names(cod) <- "D"

  return(structure(class = "sjstats_r2", list(cod = cod)))
}



#' @title Compute R-squared of (generalized) linear (mixed) models
#' @name r2
#'
#' @description Compute R-squared values of linear (mixed) models, or
#'                pseudo-R-squared values for generalized linear (mixed) models.
#'
#' @param x Fitted model of class \code{lm}, \code{glm}, \code{lmerMod}/\code{lme}
#'            or \code{glmerMod}.
#' @param n Optional, a \code{lmerMod} object, representing the fitted null-model
#'          to \code{x} (unconditional model). If \code{n} is given, the pseudo-r-squared
#'          for random intercept and random slope variances are computed (see Kwok et al. 2008;
#'          see 'Examples' and 'Details').
#'
#' @return \itemize{
#'           \item For linear models, the r-squared and adjusted r-squared values.
#'           \item For linear mixed models, the r-squared and Omega-squared values.
#'           \item For \code{glm} objects, Cox & Snell's and Nagelkerke's pseudo r-squared values.
#'           \item For \code{glmerMod} objects, Tjur's coefficient of determination.
#'         }
#'
#' @details If \code{n} is given, the Pseudo-R2 statistic is the proportion of
#'          explained variance in the random effect after adding co-variates or
#'          predictors to the model, or in short: the proportion of the explained
#'          variance in the random effect of the full (conditional) model \code{x}
#'          compared to the null (unconditional) model \code{n}.
#'
#' @note For linear models, the r-squared and adjusted r-squared value is returned,
#'         as provided by the \code{summary}-function.
#'         \cr \cr
#'         For linear mixed models, an r-squared approximation by computing the
#'         correlation between the fitted and observed values, as suggested by
#'         Byrnes (2008), is returned as well as the Omega-squared value as
#'         suggested by Xu (2003), unless \code{n} is specified. If \code{n}
#'         is given, pseudo r-squared measures based on the variances of random
#'         intercept (tau 00, between-group-variance) and random slope (tau 11,
#'         random-slope-variance) are returned.
#'         \cr \cr
#'         For generalized linear models, Cox & Snell's and Nagelkerke's
#'         pseudo r-squared values are returned.
#'         \cr \cr
#'         For generalized linear mixed models, the coefficient of determination
#'         as suggested by Tjur (2009) (see also \code{\link{cod}}).
#'
#' @references \itemize{
#'               \item \href{http://glmm.wikidot.com/faq}{DRAFT r-sig-mixed-models FAQ}
#'               \item Byrnes, J. 2008. Re: Coefficient of determination (R^2) when using lme() (\url{https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q2/000713.html})
#'               \item Kwok OM, Underhill AT, Berry JW, Luo W, Elliott TR, Yoon M. 2008. Analyzing Longitudinal Data with Multilevel Models: An Example with Individuals Living with Lower Extremity Intra-Articular Fractures. Rehabilitation Psychology 53(3): 370–86. \doi{10.1037/a0012765}
#'               \item Xu, R. 2003. Measuring explained variation in linear mixed effects models. Statist. Med. 22:3527-3541. \doi{10.1002/sim.1572}
#'               \item Tjur T. 2009. Coefficients of determination in logistic regression models - a new proposal: The coefficient of discrimination. The American Statistician, 63(4): 366-372
#'             }
#'
#' @examples
#' library(sjmisc)
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' r2(fit)
#'
#' data(efc)
#' fit <- lm(barthtot ~ c160age + c12hour, data = efc)
#' r2(fit)
#'
#' # Pseudo-R-squared values
#' efc$services <- ifelse(efc$tot_sc_e > 0, 1, 0)
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep,
#'            data = efc, family = binomial(link = "logit"))
#' r2(fit)
#'
#' # Pseudo-R-squared values for random effect variances
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' fit.null <- lmer(Reaction ~ 1 + (Days | Subject), sleepstudy)
#' r2(fit, fit.null)
#'
#'
#' @importFrom stats model.response model.frame fitted var residuals
#' @importFrom sjmisc is_empty
#' @export
r2 <- function(x, n = NULL) {
  rsq <- NULL
  osq <- NULL
  adjr2 <- NULL
  # do we have a glm? if so, report pseudo_r2
  if (any(class(x) == "glm")) {
    return(pseudo_ralt(x))
    # do we have a glmer?
  } else if (any(class(x) == "glmerMod")) {
    return(cod(x))
    # do we have a simple linear model?
  } else if (identical(class(x), "lm")) {
    rsq <- summary(x)$r.squared
    adjr2 <- summary(x)$adj.r.squared
    # name vectors
    names(rsq) <- "R2"
    names(adjr2) <- "adj.R2"
    # return results
    return(structure(class = "sjstats_r2", list(r2 = rsq, adjr2 = adjr2)))
    # else do we have a mixed model?
  } else if (any(class(x) == "plm")) {
    rsq <- summary(x)$r.squared[1]
    adjr2 <- summary(x)$r.squared[2]
    # name vectors
    names(rsq) <- "R2"
    names(adjr2) <- "adj.R2"
    # return results
    return(structure(class = "sjstats_r2", list(r2 = rsq, adjr2 = adjr2)))
  } else if (sjmisc::str_contains(class(x), pattern = c("lmerMod", "lme"),
                                  ignore.case = T, logic = "OR")) {
    # do we have null model?
    if (!is.null(n)) {
      # compute tau for both models
      tau_full <- icc(x)
      tau_null <- icc(n)
      # get taus. tau.00 is the random intercept variance, i.e. for growth models,
      # the difference in the outcome's mean at first time point
      rsq0 <- (attr(tau_null, "tau.00") - attr(tau_full, "tau.00")) / attr(tau_null, "tau.00")
      # tau.11 is the variance of the random slopes, i.e. how model predictors
      # affect the trajectory of subjects over time (for growth models)
      rsq1 <- (attr(tau_null, "tau.11") - attr(tau_full, "tau.11")) / attr(tau_null, "tau.11")
      # if model has no random slope, we need to set this value to NA
      if (is.null(rsq1) || sjmisc::is_empty(rsq1)) rsq1 <- NA
      # name vectors
      names(rsq0) <- "R2(tau-00)"
      names(rsq1) <- "R2(tau-11)"
      # return results
      return(structure(class = "sjstats_r2", list(r2_tau00 = rsq0, r2_tau11 = rsq1)))
    } else {
      # compute "correlation"
      lmfit <-  lm(resp_val(x) ~ stats::fitted(x))
      # get r-squared
      rsq <- summary(lmfit)$r.squared
      # get omega squared
      osq <- 1 - stats::var(stats::residuals(x)) / stats::var(resp_val(x))
      # name vectors
      names(rsq) <- "R2"
      names(osq) <- "O2"
      # return results
      return(structure(class = "sjstats_r2", list(r2 = rsq, o2 = osq)))
    }
  } else {
    warning("`r2` only works on linear (mixed) models of class \"lm\", \"lme\" or \"lmerMod\".", call. = F)
    return(NULL)
  }
}

#' @importFrom stats nobs deviance
pseudo_ralt <- function(x) {
  # get nr of observations
  n <- stats::nobs(x)
  CoxSnell <- 1 - exp((stats::deviance(x) - x$null.deviance) / n)
  Nagelkerke <- CoxSnell / (1 - exp(-x$null.deviance / n))
  names(CoxSnell) <- "CoxSnell"
  names(Nagelkerke) <- "Nagelkerke"
  return(structure(class = "sjstats_r2", list(CoxSnell = CoxSnell, Nagelkerke = Nagelkerke)))
}
