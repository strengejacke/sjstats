#' @title Robust standard errors for regression models
#' @name robust
#' @description \code{robust()} computes robust standard error for regression models.
#'                This method wraps the \code{\link[lmtest]{coeftest}}-function with
#'                robust covariance matrix estimators based on the
#'                \code{\link[sandwich]{vcovHC}}-function, and returns the
#'                result as tidy data frame.
#'                \cr \cr
#'                \code{svy()} is intended to compute standard errors for survey
#'                designs (complex samples) fitted with regular \code{lm} or
#'                \code{glm} functions, as alternative to the \pkg{survey}-package.
#'                It simulates sampling weights by adjusting the residual degrees
#'                of freedom based on the precision weights used to fit \code{x},
#'                and then calls \code{robust()} with the adjusted model.
#'
#' @param x A fitted model of any class that is supported by the \code{coeftest()}-function.
#'          For \code{svy()}, \code{x} must be \code{lm} object, fitted with weights.
#' @param vcov Character vector, specifying the estimation type for the
#'          heteroskedasticity-consistent covariance matrix estimation
#'          (see \code{\link[sandwich]{vcovHC}} for details).
#' @param conf.int Logical, \code{TRUE} if confidence intervals based on robust
#'          standard errors should be included.
#' @param exponentiate Logical, whether to exponentiate the coefficient estimates
#'          and confidence intervals (typical for logistic regression).
#'
#' @return A summary of the model, including estimates, robust standard error,
#'           p-value and - optionally - the confidence intervals.
#'
#' @note \code{svy()} simply calls \code{robust()}, but first adjusts the
#'       residual degrees of freedom based on the model weights.
#'       Hence, for \code{svy()}, \code{x} should be fitted with weights.
#'       This simulates \emph{sampling weights} like in survey designs, though
#'       \code{lm} and \code{glm} implement \emph{precision weights}.
#'       The results from \code{svy()} are usually more accurate than simple
#'       weighted standard errors for complex samples. However, results from
#'       the \pkg{survey} package are still more exactly, especially
#'       regarding the estimates.
#'       \cr \cr
#'       \code{vcov} for \code{svy()} defaults to \code{"HC1"}, because
#'       standard errors with this estimation type come closest to the standard
#'       errors from the \pkg{survey}-package.
#'       \cr \cr
#'       Currently, \code{svy()} only works for objects of class \code{lm}.
#'
#' @examples
#' data(efc)
#' fit <- lm(barthtot ~ c160age + c12hour + c161sex + c172code, data = efc)
#' summary(fit)
#' robust(fit)
#'
#' confint(fit)
#' robust(fit, conf.int = TRUE)
#' robust(fit, vcov = "HC1", conf.int = TRUE) # "HC1" should be Stata default
#'
#' library(sjmisc)
#' # dichtomozize service usage by "service usage yes/no"
#' efc$services <- sjmisc::dicho(efc$tot_sc_e, dich.by = 0)
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep,
#'            data = efc, family = binomial(link = "logit"))
#'
#' robust(fit)
#' robust(fit, conf.int = TRUE, exponentiate = TRUE)
#'
#' @importFrom stats qt
#' @importFrom lmtest coeftest
#' @importFrom sandwich vcovHC
#' @importFrom tibble tibble add_column
#' @export
robust <- function(x, vcov = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"), conf.int = FALSE, exponentiate = FALSE) {
  # match arguments
  vcov <- match.arg(vcov)

  # compute robust standard errors
  tmp <- lmtest::coeftest(x, vcov. = sandwich::vcovHC(x, type = vcov))

  # create tidy tibble
  result <- tibble::tibble(
    term = rownames(tmp),
    estimate = tmp[, 1],
    std.error = tmp[, 2],
    statistic = tmp[, 3],
    p.value = tmp[, 4]
  )

  # add CI
  if (conf.int) {
    # denominator df
    dendf <- summary(x)$df[2]
    # add columns with CI
    result <- tibble::add_column(
      result,
      conf.low = result$estimate - (stats::qt(.975, df = dendf) * result$std.error),
      conf.high = result$estimate + stats::qt(.975, df = dendf) * result$std.error
    )
  }

  # exponentiate results?
  if (exponentiate) {
    result$estimate <- exp(result$estimate)
    result$conf.low <- exp(result$conf.low)
    result$conf.high <- exp(result$conf.high)
  }

  result
}


#' @rdname robust
#' @importFrom stats weights
#' @export
svy <- function(x, vcov = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5"), conf.int = FALSE, exponentiate = FALSE) {
  # match arguments
  vcov <- match.arg(vcov)

  # check if we have lm-object
  if (inherits(x, "lm", which = TRUE) == 1) {

    # check if model has weights
    w <- stats::weights(x)

    if (!is.null(w))
      # re-weight residuals
      x$df.residual <- with(x, sum(weights) - length(coefficients))
    else
      warning("Model has no weights. Computing robust standard error for non-weighted model.", call. = F)

  } else {
    # no sampling-weights adjustment for other models than lm right now
    warning("`x` must be of class `lm`. Computing robust standard errors now without adjusting residual df.", call. = F)
  }

  # compute robust se
  suppressWarnings(robust(x, vcov, conf.int, exponentiate))
}
