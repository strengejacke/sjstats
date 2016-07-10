#' @title Coefficient of Variation
#' @name cv
#' @description Compute coefficient of variation for single variables
#'                (standard deviation divided by mean) or for fitted
#'                linear (mixed effects) models (root mean squared error
#'                (RMSE) divided by mean of dependent variable).
#'
#' @param x (Numeric) vector or a fitted linear model of class
#'          \code{\link{lm}}, \code{\link[lme4]{merMod}} (\pkg{lme4}) or
#'          \code{\link[nlme]{lme}} (\pkg{nlme}).
#' @param ... More fitted model objects, to compute multiple coefficients of
#'              variation at once.
#' @return The coefficient of variation of \code{x}.
#'
#' @details The advantage of the cv is that it is unitless. This allows
#'            coefficient of variation to be compared to each other in ways
#'            that other measures, like standard deviations or root mean
#'            squared residuals, cannot be.
#'            \cr \cr
#'            "It is interesting to note the differences between a model's CV
#'            and R-squared values. Both are unitless measures that are indicative
#'            of model fit, but they define model fit in two different ways: CV
#'            evaluates the relative closeness of the predictions to the actual
#'            values while R-squared evaluates how much of the variability in the
#'            actual values is explained by the model"
#'            (\href{http://www.ats.ucla.edu/stat/mult_pkg/faq/general/coefficient_of_variation.htm}{source: UCLA-FAQ}).
#'
#' @seealso \code{\link{rmse}}
#'
#' @references \itemize{
#'               \item \href{http://www.ats.ucla.edu/stat/mult_pkg/faq/general/coefficient_of_variation.htm}{UCLA-FAQ: What is the coefficient of variation?}
#'               \item  Everitt, Brian (1998). The Cambridge Dictionary of Statistics. Cambridge, UK New York: Cambridge University Press
#'             }
#'
#' @examples
#' data(efc)
#' cv(efc$e17age)
#'
#' fit <- lm(neg_c_7 ~ e42dep, data = efc)
#' cv(fit)
#'
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' cv(fit)
#'
#' library(nlme)
#' fit <- lme(distance ~ age, data = Orthodont)
#' cv(fit)
#'
#' @importFrom stats sd
#' @export
cv <- function(x, ...) {
  # return value
  cv_ <- cv_helper(x)
  # check if we have multiple parameters
  if (nargs() > 1) {
    # get input list
    params_ <- list(...)
    cv_ <- c(cv_, sapply(params_, cv_helper))
  }
  return(cv_)
}


cv_helper <- function(x) {
  # check if we have a fitted linear model
  if (class(x) == "lm" || any(class(x) == "lmerMod") ||
      any(class(x) == "lme") || any(class(x) == "merModLmerTest")) {
    if (class(x) == "lm") {
      # dependent variable in lm
      dv <- x$model[[1]]
    } else if (any(class(x) == "lmerMod") || any(class(x) == "merModLmerTest")) {
      # dependent variable in lmerMod
      dv <- lme4::getME(x, "y")
    } else if (any(class(x) == "lme")) {
      # dependent variable in lme
      dv <- unname(nlme::getResponse(x))
    }
    # compute mean of dependent variable
    mw <- mean(dv, na.rm = TRUE)
    # check if mean is zero?
    if (mw != 0) {
      # cv = root mean squared error (RMSE) divided by mean of dep. var.
      return(rmse(x) / mw)
    } else {
      warning("Mean of dependent variable is zero. Cannot compute model's coefficient of variation.", call. = F)
    }
  } else {
    # compute mean of variable
    mw <- mean(x, na.rm = TRUE)
    # check if mean is zero?
    if (mw != 0) {
      #  we assume a simple vector
      return(stats::sd(x, na.rm = TRUE) / mw)
    } else {
      warning("Mean of `x` is zero. Cannot compute coefficient of variation.", call. = F)
    }
  }
  return(NULL)
}