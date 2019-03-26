#' @title Standard error and confidence intervals for bootstrapped estimates
#' @name boot_ci
#'
#' @description Compute nonparametric bootstrap estimate, standard error,
#'              confidence intervals and p-value for a vector of bootstrap
#'              replicate estimates.
#'
#' @param data A data frame that containts the vector with bootstrapped
#'          estimates, or directly the vector (see 'Examples').
#' @param ... Optional, unquoted names of variables with bootstrapped estimates.
#'          Required, if either \code{data} is a data frame (and no vector),
#'          and only selected variables from \code{data} should be processed.
#'          You may also use functions like \code{:} or tidyselect's
#'          \code{\link[tidyselect]{select_helpers}}.
#' @param method Character vector, indicating if confidence intervals should be
#'          based on bootstrap standard error, multiplied by the value of the
#'          quantile function of the t-distribution (default), or on sample
#'          quantiles of the bootstrapped values. See 'Details' in \code{boot_ci()}.
#'          May be abbreviated.
#' @inheritParams std_beta
#'
#' @return A \code{\link[tibble]{tibble}} with either bootstrap estimate,
#'         standard error, the lower and upper confidence intervals or the
#'         p-value for all bootstrapped estimates.
#'
#' @details The methods require one or more vectors of bootstrap replicate estimates
#'          as input.
#'          \itemize{
#'            \item{
#'              \code{boot_est()} returns the bootstrapped estimate, simply by
#'              computing the mean value of all bootstrap estimates.
#'            }
#'            \item{
#'              \code{boot_se()} computes the nonparametric bootstrap standard
#'              error by calculating the standard deviation of the input vector.
#'            }
#'            \item{
#'              The mean value of the input vector and its standard error is used
#'              by \code{boot_ci()} to calculate the lower and upper confidence
#'              interval, assuming a t-distribution of bootstrap estimate replicates
#'              (for \code{method = "dist"}, the default, which is
#'              \code{mean(x) +/- qt(.975, df = length(x) - 1) * sd(x)}); for
#'              \code{method = "quantile"}, 95\% sample quantiles are used to compute
#'              the confidence intervals (\code{quantile(x, probs = c(.025, .975))}).
#'              Use \code{ci.lvl} to change the level for the confidence interval.
#'            }
#'            \item{
#'              P-values from \code{boot_p()} are also based on t-statistics,
#'              assuming normal distribution.
#'            }
#'          }
#'
#' @references Carpenter J, Bithell J. Bootstrap confdence intervals: when, which, what? A practical guide for medical statisticians. Statist. Med. 2000; 19:1141-1164
#'
#' @seealso \code{\link{bootstrap}} to generate nonparametric bootstrap samples.
#'
#' @examples
#' library(dplyr)
#' library(purrr)
#' data(efc)
#' bs <- bootstrap(efc, 100)
#'
#' # now run models for each bootstrapped sample
#' bs$models <- map(bs$strap, ~lm(neg_c_7 ~ e42dep + c161sex, data = .x))
#'
#' # extract coefficient "dependency" and "gender" from each model
#' bs$dependency <- map_dbl(bs$models, ~coef(.x)[2])
#' bs$gender <- map_dbl(bs$models, ~coef(.x)[3])
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
#' boot_ci(bs, dependency)
#' boot_ci(bs, dependency, gender)
#' boot_ci(bs, dependency, gender, method = "q")
#'
#'
#' # compare coefficients
#' mean(bs$dependency)
#' boot_est(bs$dependency)
#' coef(fit)[2]
#'
#'
#' # bootstrap() and boot_ci() work fine within pipe-chains
#' efc %>%
#'   bootstrap(100) %>%
#'   mutate(
#'     models = map(strap, ~lm(neg_c_7 ~ e42dep + c161sex, data = .x)),
#'     dependency = map_dbl(models, ~coef(.x)[2])
#'   ) %>%
#'   boot_ci(dependency)
#'
#' # check p-value
#' boot_p(bs$gender)
#' summary(fit)$coefficients[3, ]
#'
#' \dontrun{
#' # 'spread_coef()' from the 'sjmisc'-package makes it easy to generate
#' # bootstrapped statistics like confidence intervals or p-values
#' library(dplyr)
#' library(sjmisc)
#' efc %>%
#'   # generate bootstrap replicates
#'   bootstrap(100) %>%
#'   # apply lm to all bootstrapped data sets
#'   mutate(
#'     models = map(strap, ~lm(neg_c_7 ~ e42dep + c161sex + c172code, data = .x))
#'   ) %>%
#'   # spread model coefficient for all 100 models
#'   spread_coef(models) %>%
#'   # compute the CI for all bootstrapped model coefficients
#'   boot_ci(e42dep, c161sex, c172code)
#'
#' # or...
#' efc %>%
#'   # generate bootstrap replicates
#'   bootstrap(100) %>%
#'   # apply lm to all bootstrapped data sets
#'   mutate(
#'     models = map(strap, ~lm(neg_c_7 ~ e42dep + c161sex + c172code, data = .x))
#'   ) %>%
#'   # spread model coefficient for all 100 models
#'   spread_coef(models, append = FALSE) %>%
#'   # compute the CI for all bootstrapped model coefficients
#'   boot_ci()}
#'
#' @importFrom stats qt quantile
#' @importFrom dplyr quos
#' @importFrom rlang .data
#' @export
boot_ci <- function(data, ..., method = c("dist", "quantile"), ci.lvl = .95) {
  # match arguments
  method <- match.arg(method)

  # evaluate arguments, generate data
  .dat <- get_dot_data(data, dplyr::quos(...))

  # compute confidence intervalls for all values
  transform_boot_result(lapply(.dat, function(x) {
    # check if method should be based on t-distribution of
    # bootstrap values or quantiles
    if (method == "dist") {
      # get bootstrap standard error
      bootse <- stats::qt((1 + ci.lvl) / 2, df = length(x) - 1) * stats::sd(x, na.rm = T)
      # lower and upper confidence interval
      ci <- mean(x, na.rm = T) + c(-bootse, bootse)
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
#' @importFrom stats sd
#' @export
boot_se <- function(data, ...) {
  # evaluate arguments, generate data
  .dat <- get_dot_data(data, dplyr::quos(...))

  # compute confidence intervalls for all values
  transform_boot_result(lapply(.dat, function(x) {
    # get bootstrap standard error
    se <- stats::sd(x, na.rm = T)
    names(se) <- "std.err"
    se
  }))
}


#' @rdname boot_ci
#' @importFrom stats sd pt
#' @export
boot_p <- function(data, ...) {
  # evaluate arguments, generate data
  .dat <- get_dot_data(data, dplyr::quos(...))

  # compute confidence intervalls for all values
  transform_boot_result(lapply(.dat, function(x) {
    # compute t-statistic
    t.stat <- mean(x, na.rm = T) / stats::sd(x, na.rm = T)
    # compute p-value
    p <- 2 * stats::pt(abs(t.stat), df = length(x) - 1, lower.tail = FALSE)
    names(p) <- "p.value"
    p
  }))
}


#' @rdname boot_ci
#' @export
boot_est <- function(data, ...) {
  # evaluate arguments, generate data
  .dat <- get_dot_data(data, dplyr::quos(...))

  # compute mean for all values (= bootstrapped estimate)
  transform_boot_result(lapply(.dat, function(x) {
    estimate <- mean(x, na.rm = T)
    names(estimate) <- "estimate"
    estimate
  }))
}



transform_boot_result <- function(res) {
  # transform a bit, so we have each estimate in a row, and ci's as columns...
  res %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    rownames_as_column(var = "term")
}


#' @importFrom dplyr select
get_dot_data <- function(x, qs) {
  if (sjmisc::is_empty(qs))
    as.data.frame(x)
  else
    suppressWarnings(dplyr::select(x, !!!qs))
}
