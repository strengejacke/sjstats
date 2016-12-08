#' @title Robust standard errors for regression models
#' @name robust
#' @description Compute robust standard error for regression models. This method
#'                wraps the \code{\link[lmtest]{coeftest}}-function with
#'                robust covariance matrix estimators based on the
#'                \code{\link[sandwich]{vcovHC}}-function, and returns the
#'                result as tidy data frame.
#'
#' @param x A fitted model of any class that is supported by the \code{coeftest()}-function.
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
  result <- tibble::tibble(term = rownames(tmp),
                           estimate = tmp[, 1],
                           std.error = tmp[, 2],
                           statistic = tmp[, 3],
                           p.value = tmp[, 4])

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
