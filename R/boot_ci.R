#' @title Standard error and confidence intervals for bootstrapped estimates
#' @name boot_ci
#'
#' @description Compute nonparametric bootstrap estimate, standard error,
#'              confidence intervals and p-value for a vector of bootstrap
#'              replicate estimates.
#'
#' @param data A data frame that containts the vector with bootstrapped
#' estimates, or directly the vector (see 'Examples').
#' @param ci.lvl Numeric, the level of the confidence intervals.
#' @param select Optional, unquoted names of variables (as character vector)
#' with bootstrapped estimates. Required, if either `data` is a data frame
#' (and no vector), and only selected variables from `data` should be processed.
#' @param method Character vector, indicating if confidence intervals should be
#' based on bootstrap standard error, multiplied by the value of the quantile
#' function of the t-distribution (default), or on sample quantiles of the
#' bootstrapped values. See 'Details' in `boot_ci()`. May be abbreviated.
#'
#' @return A data frame with either bootstrap estimate, standard error, the
#' lower and upper confidence intervals or the p-value for all bootstrapped
#' estimates.
#'
#' @details The methods require one or more vectors of bootstrap replicate
#' estimates as input.
#'
#' - `boot_est()`: returns the bootstrapped estimate, simply by computing
#'   the mean value of all bootstrap estimates.
#' - `boot_se()`: computes the nonparametric bootstrap standard error by
#'   calculating the standard deviation of the input vector.
#' - The mean value of the input vector and its standard error is used by
#'   `boot_ci()` to calculate the lower and upper confidence interval,
#'   assuming a t-distribution of bootstrap estimate replicates (for
#'   `method = "dist"`, the default, which is
#'   `mean(x) +/- qt(.975, df = length(x) - 1) * sd(x)`); for
#'   `method = "quantile"`, 95\% sample quantiles are used to compute the
#'   confidence intervals (`quantile(x, probs = c(0.025, 0.975))`). Use
#'   `ci.lvl` to change the level for the confidence interval.
#' - P-values from `boot_p()` are also based on t-statistics, assuming normal
#'   distribution.
#'
#' @references Carpenter J, Bithell J. Bootstrap confdence intervals: when, which, what? A practical guide for medical statisticians. Statist. Med. 2000; 19:1141-1164
#'
#' @seealso []`bootstrap()`] to generate nonparametric bootstrap samples.
#'
#' @examples
#' data(efc)
#' bs <- bootstrap(efc, 100)
#'
#' # now run models for each bootstrapped sample
#' bs$models <- lapply(
#'   bs$strap,
#'   function(.x) lm(neg_c_7 ~ e42dep + c161sex, data = .x)
#' )
#'
#' # extract coefficient "dependency" and "gender" from each model
#' bs$dependency <- vapply(bs$models, function(x) coef(x)[2], numeric(1))
#' bs$gender <- vapply(bs$models, function(x) coef(x)[3], numeric(1))
#'
#' # get bootstrapped confidence intervals
#' boot_ci(bs$dependency)
#'
#' # compare with model fit
#' fit <- lm(neg_c_7 ~ e42dep + c161sex, data = efc)
#' confint(fit)[2, ]
#'
#' # alternative function calls.
#' boot_ci(bs$dependency)
#' boot_ci(bs, "dependency")
#' boot_ci(bs, c("dependency", "gender"))
#' boot_ci(bs, c("dependency", "gender"), method = "q")
#'
#'
#' # compare coefficients
#' mean(bs$dependency)
#' boot_est(bs$dependency)
#' coef(fit)[2]
#' @export
boot_ci <- function(data, select = NULL, method = c("dist", "quantile"), ci.lvl = 0.95) {
  # match arguments
  method <- match.arg(method)

  # evaluate arguments, generate data
  if (is.null(select)) {
    .dat <- as.data.frame(data)
  } else {
    .dat <- data[select]
  }

  # compute confidence intervals for all values
  transform_boot_result(lapply(.dat, function(x) {
    # check if method should be based on t-distribution of
    # bootstrap values or quantiles
    if (method == "dist") {
      # get bootstrap standard error
      bootse <- stats::qt((1 + ci.lvl) / 2, df = length(x) - 1) * stats::sd(x, na.rm = TRUE)
      # lower and upper confidence interval
      ci <- mean(x, na.rm = TRUE) + c(-bootse, bootse)
    } else {
      # CI based on quantiles of bootstrapped values
      ci <- stats::quantile(x, probs = c((1 - ci.lvl) / 2, (1 + ci.lvl) / 2))
    }
    # give proper names
    names(ci) <- c("conf.low", "conf.high")
    ci
  }))
}


#' @rdname boot_ci
#' @export
boot_se <- function(data, select = NULL) {
  # evaluate arguments, generate data
  if (is.null(select)) {
    .dat <- as.data.frame(data)
  } else {
    .dat <- data[select]
  }
  # compute confidence intervalls for all values
  transform_boot_result(lapply(.dat, function(x) {
    # get bootstrap standard error
    se <- stats::sd(x, na.rm = TRUE)
    names(se) <- "std.err"
    se
  }))
}


#' @rdname boot_ci
#' @export
boot_p <- function(data, select = NULL) {
  # evaluate arguments, generate data
  if (is.null(select)) {
    .dat <- as.data.frame(data)
  } else {
    .dat <- data[select]
  }
  # compute confidence intervalls for all values
  transform_boot_result(lapply(.dat, function(x) {
    # compute t-statistic
    t.stat <- mean(x, na.rm = TRUE) / stats::sd(x, na.rm = TRUE)
    # compute p-value
    p <- 2 * stats::pt(abs(t.stat), df = length(x) - 1, lower.tail = FALSE)
    names(p) <- "p.value"
    p
  }))
}


#' @rdname boot_ci
#' @export
boot_est <- function(data, select = NULL) {
  # evaluate arguments, generate data
  if (is.null(select)) {
    .dat <- as.data.frame(data)
  } else {
    .dat <- data[select]
  }
  # compute mean for all values (= bootstrapped estimate)
  transform_boot_result(lapply(.dat, function(x) {
    estimate <- mean(x, na.rm = TRUE)
    names(estimate) <- "estimate"
    estimate
  }))
}


transform_boot_result <- function(res) {
  # transform a bit, so we have each estimate in a row, and ci's as columns...
  rownames_as_column(as.data.frame(t(as.data.frame(res))), var = "term")
}
