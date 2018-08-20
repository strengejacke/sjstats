#' @title Robust standard errors for regression models
#' @name robust
#' @description \code{robust()} computes robust standard error for regression models.
#'    This method calls one of the \code{vcov*()}-functions from the
#'    \pkg{sandwich}-package for robust covariance matrix estimators. Results are
#'    returned as tidy data frame.
#'    \cr \cr
#'    \code{svy()} is intended to compute standard errors for survey
#'    designs (complex samples) fitted with regular \code{lm} or
#'    \code{glm} functions, as alternative to the \pkg{survey}-package.
#'    It simulates sampling weights by adjusting the residual degrees
#'    of freedom based on the precision weights used to fit \code{x},
#'    and then calls \code{robust()} with the adjusted model.
#'
#' @param x A fitted model of any class that is supported by the \code{vcov*()}-functions
#'    from the \pkg{sandwich} package. For \code{svy()}, \code{x} must be
#'    \code{lm} object, fitted with weights.
#' @param vcov.fun String, indicating the name of the \code{vcov*()}-function
#'    from the \pkg{sandwich}-package, e.g. \code{vcov.fun = "vcovCL"}.
#' @param vcov.type Character vector, specifying the estimation type for the
#'    robust covariance matrix estimation (see \code{\link[sandwich]{vcovHC}}
#'    for details).
#' @param vcov.args List of named vectors, used as additional arguments that
#'    are passed down to \code{vcov.fun}.
#' @param conf.int Logical, \code{TRUE} if confidence intervals based on robust
#'    standard errors should be included.
#' @param exponentiate Logical, whether to exponentiate the coefficient estimates
#'    and confidence intervals (typical for logistic regression).
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
#'       \code{vcov.type} for \code{svy()} defaults to \code{"HC1"}, because
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
#' robust(fit, vcov.type = "HC1", conf.int = TRUE) # "HC1" should be Stata default
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
#' @importFrom stats qt pt df.residual qnorm pnorm nobs coef
#' @export
robust <- function(x, vcov.fun = "vcovHC", vcov.type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"), vcov.args = NULL, conf.int = FALSE, exponentiate = FALSE) {

  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("Package `sandwich` needed for this function. Please install and try again.")
  }

  # match arguments
  vcov.type <- match.arg(vcov.type)

  # get coefficients
  est <- stats::coef(x)

  # compute robust standard errors based on vcov
  vcov.fun <- get(vcov.fun, asNamespace("sandwich"))
  .vcov <- do.call(vcov.fun, c(list(x = x, type = vcov.type), vcov.args))

  se <- sqrt(diag(.vcov))

  dendf <- tryCatch(
    stats::df.residual(x),
    error = function(x) { NULL },
    warning = function(x) { NULL },
    finally = function(x) { NULL }
  )

  # 2nd try
  if (is.null(dendf)) {
    dendf <- tryCatch(
      summary(x)$df[2],
      error = function(x) { NULL },
      warning = function(x) { NULL },
      finally = function(x) { NULL }
    )
  }

  # 3rd try
  if (is.null(dendf)) {
    dendf <- tryCatch(
      stats::nobs(x) - length(est),
      error = function(x) { NULL },
      warning = function(x) { NULL },
      finally = function(x) { NULL }
    )
  }


  t.stat <- est / se

  if (is.null(dendf)) {
    p.value <- 2 * stats::pnorm(abs(t.stat), lower.tail = FALSE)
    se.factor <- stats::qnorm(.975)
  } else {
    p.value <- 2 * stats::pt(abs(t.stat), df = dendf, lower.tail = FALSE)
    se.factor <- stats::qt(.975, df = dendf)
  }


  # create tidy data frame
  result <- data_frame(
    term = names(est),
    estimate = est,
    std.error = se,
    statistic = t.stat,
    p.value = p.value
  )

  # add CI
  if (conf.int) {
    # add columns with CI
    result <- add_cols(
      result,
      conf.low = result$estimate - se.factor * result$std.error,
      conf.high = result$estimate + se.factor * result$std.error,
      .after = "std.error"
    )
  }

  # exponentiate results?
  if (exponentiate) {
    result$estimate <- exp(result$estimate)
    if (obj_has_name(result, "conf.low")) result$conf.low <- exp(result$conf.low)
    if (obj_has_name(result, "conf.high")) result$conf.high <- exp(result$conf.high)
  }

  result
}


#' @rdname robust
#' @importFrom stats weights
#' @export
svy <- function(x, vcov.fun = "vcovHC", vcov.type = c("HC1", "const", "HC", "HC0", "HC3", "HC2", "HC4", "HC4m", "HC5"), vcov.args = NULL, conf.int = FALSE, exponentiate = FALSE) {
  # match arguments
  vcov.type <- match.arg(vcov.type)

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
  suppressWarnings(robust(
    x,
    vcov.fun = vcov.fun,
    vcov.type = vcov.type,
    vcov.args = vcov.args,
    conf.int = conf.int,
    exponentiate = exponentiate
  ))
}
