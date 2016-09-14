#' @title Standard Error and Confidence Intervals for bootstrapped estimates
#' @name boot_ci
#'
#' @description Compute nonparametric bootstrap standard error, confidence
#'              intervals and p-value for a vector of bootstrap replicate
#'              estimates.
#'
#' @param data A data frame that containts the vector with bootstrapped
#'          estimates, or directly the vector (see 'Examples').
#' @param ... Optional, names of the variables with bootstrapped estimates.
#'          Required, if either \code{data} is a data frame and no vector,
#'          or if only selected variables from \code{data} should be used
#'          in the function.
#'
#' @return A \code{\link[tibble]{data_frame}} with either bootstrap standard error,
#'         the lower and upper confidence intervals or the p-value for all
#'         bootstrapped estimates.
#'
#' @details This method requires one or more vectors of bootstrap replicate estimates
#'          as input. The function then computes the nonparametric bootstrap
#'          standard error by calculating the standard deviation of the input
#'          vector. The mean value of the input vector and its standard error is
#'          used to calculate the lower and upper confidence interval, assuming
#'          a t-distribution of bootstrap estimate replicates. P-values
#'          from \code{boot_p} are also based on t-statistics, assuming normal
#'          distribution.
#'
#' @seealso \code{\link{bootstrap}} to generate nonparametric bootstrap samples.
#'
#' @examples
#' data(efc)
#' bs <- bootstrap(efc, 100)
#'
#' # now run models for each bootstrapped sample
#' bs$models <- lapply(bs$strap, function(x) lm(neg_c_7 ~ e42dep + c161sex, data = x))
#'
#' # extract coefficient "dependency" and "gender" from each model
#' bs$dependency <- unlist(lapply(bs$models, function(x) coef(x)[2]))
#' bs$gender <- unlist(lapply(bs$models, function(x) coef(x)[3]))
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
#'
#'
#' # compare coefficients
#' mean(bs$dependency)
#' coef(fit)[2]
#'
#' # bootstrap() and boot_ci() work fine within pipe-chains
#' library(dplyr)
#' efc %>%
#'   bootstrap(100) %>%
#'   mutate(models = lapply(.$strap, function(x) {
#'     lm(neg_c_7 ~ e42dep + c161sex, data = x)
#'   })) %>%
#'   mutate(dependency = unlist(lapply(.$models, function(x) coef(x)[2]))) %>%
#'   boot_ci(dependency)
#'
#'
#' # check p-value
#' boot_p(bs$gender)
#' summary(fit)$coefficients[3, ]
#'
#'
#' # 'spread_coef()' from the 'sjmisc'-package makes it easy to generate
#' # bootstrapped statistics like confidence intervals or p-values
#' library(dplyr)
#' library(sjmisc)
#' efc %>%
#'   # generate bootstrap replicates
#'   bootstrap(100) %>%
#'   # apply lm to all bootstrapped data sets
#'   mutate(models = lapply(.$strap, function(x) {
#'     lm(neg_c_7 ~ e42dep + c161sex + c172code, data = x)
#'   })) %>%
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
#'   mutate(models = lapply(strap, function(x) {
#'     lm(neg_c_7 ~ e42dep + c161sex + c172code, data = x)
#'   })) %>%
#'   # spread model coefficient for all 100 models
#'   spread_coef(models, append = FALSE) %>%
#'   # compute the CI for all bootstrapped model coefficients
#'   boot_ci()
#'
#' @importFrom stats qt
#' @export
boot_ci <- function(data, ...) {
  # evaluate arguments, generate data
  .dat <- get_boot_data(data, match.call(expand.dots = FALSE)$`...`)
  # compute confidence intervalls for all values
  transform_boot_result(lapply(.dat, function(x) {
    # get bootstrap standard error
    bootse <- stats::qt(.975, df = length(x) - 1) * stats::sd(x, na.rm = T)
    # lower and upper confidence interval
    ci <- mean(x, na.rm = T) + c(-bootse, bootse)
    names(ci) <- c("conf.low", "conf.high")
    ci
  }))
}


#' @rdname boot_ci
#' @importFrom stats sd
#' @export
boot_se <- function(data, ...) {
  # evaluate arguments, generate data
  .dat <- get_boot_data(data, match.call(expand.dots = FALSE)$`...`)
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
  .dat <- get_boot_data(data, match.call(expand.dots = FALSE)$`...`)
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



#' @importFrom tibble rownames_to_column
transform_boot_result <- function(res) {
  # transform a bit, so we have each estimate in a row, and ci's as columns...
  tibble::rownames_to_column(as.data.frame(t(as.data.frame(res))), var = "term")
}


get_boot_data <- function(data, dots) {
  # any dots?
  if (length(dots) > 0) {
    # get variable names
    vars <- dot_names(dots)
  } else {
    vars <- NULL
  }

  # check if data is a data frame
  if (is.data.frame(data)) {
    # do we have any variables specified?
    if (!is.null(vars))
      x <- data[, vars, drop = FALSE]
    else
      x <- data
  } else {
    x <- as.data.frame(data)
  }
  x
}
