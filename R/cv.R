#' @title Compute model quality
#' @name cv
#'
#' @description Compute the coefficient of variation.
#'
#' @param x Fitted linear model of class \code{lm}, \code{merMod} (\pkg{lme4})
#'   or \code{lme} (\pkg{nlme}).
#' @param ... More fitted model objects, to compute multiple coefficients of
#'   variation at once.
#'
#' @details The advantage of the cv is that it is unitless. This allows
#'         coefficient of variation to be compared to each other in ways
#'         that other measures, like standard deviations or root mean
#'         squared residuals, cannot be.
#'
#' @return Numeric, the coefficient of variation.
#'
#' @examples
#' data(efc)
#' fit <- lm(barthtot ~ c160age + c12hour, data = efc)
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

  cv_
}


#' @importFrom performance rmse
#' @importFrom insight get_response
cv_helper <- function(x) {
  # check if we have a fitted linear model
  if (inherits(x, c("lm", "lmerMod", "lme", "merModLmerTest")) && !inherits(x, "glm")) {
    # get response
    dv <- insight::get_response(x)
    mw <- mean(dv, na.rm = TRUE)
    stddev <- performance::rmse(x)
  } else {
    mw <- mean(x, na.rm = TRUE)
    stddev <- stats::sd(x, na.rm = TRUE)
  }

  # check if mean is zero?
  if (mw == 0)
    stop("Mean of dependent variable is zero. Cannot compute model's coefficient of variation.", call. = F)

  stddev / mw
}
