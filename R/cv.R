#' @rdname rmse
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


cv_helper <- function(x) {
  # check if we have a fitted linear model
  if (inherits(x, c("lm", "lmerMod", "lme", "merModLmerTest")) && !inherits(x, "glm")) {
    # get response
    dv <- resp_val(x)
    mw <- mean(dv, na.rm = TRUE)
    stddev <- rmse(x)
  } else {
    mw <- mean(x, na.rm = TRUE)
    stddev <- stats::sd(x, na.rm = TRUE)
  }

  # check if mean is zero?
  if (mw == 0)
    stop("Mean of dependent variable is zero. Cannot compute model's coefficient of variation.", call. = F)

  stddev / mw
}
