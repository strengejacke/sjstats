#' @title Generate bootstrap replications
#' @name bootstrap
#'
#' @description Generates \code{n} bootstrap samples of \code{data} and
#'              returns the bootstrapped data frames as list-variable.
#'
#' @param data A data frame.
#' @param n Number of bootstraps to be generated
#' @param size Optional, size of the bootstrap samples. May either be a number
#'          between 1 and \code{nrow(data)} or a value between 0 and 1 to sample
#'          a proportion of observations from \code{data} (see 'Examples').
#'
#' @return A \code{\link[tibble]{data_frame}} with one column: a list-variable
#'           \code{strap}, which contains the bootstrapped samples from \code{data}.
#'
#' @details By default, each bootstrap sample has the same number of observations
#'            as \code{data}. To generate bootstrap samples without resampling
#'            same observations (i.e. sampling without replacement), use
#'            \code{size} to get bootstrapped data with a specific number
#'            of observations.
#'
#'
#' @examples
#' data(efc)
#' bs <- bootstrap(efc, 5)
#'
#' # now run models for each bootstrapped sample
#' lapply(bs$strap, function(x) lm(neg_c_7 ~ e42dep + c161sex, data = x))
#'
#' # generate bootstrap samples with 600 observations for each sample
#' bs <- bootstrap(efc, 5, 600)
#'
#' # generate bootstrap samples with 70% observations of the original sample size
#' bs <- bootstrap(efc, 5, .7)
#'
#'
#' @importFrom tibble data_frame
#' @export
bootstrap <- function(data, n, size) {
  if (!missing(size) && !is.null(size)) {
    # check for valid range
    if (size < 0 || size > nrow(data))
      stop("`size` must be greater than 0, but not greater than number of rows of `data`.", call. = F)
    # check if we want proportions
    if (size < 1) size <- as.integer(nrow(data) * size)
    # generate bootstraps
    res <- replicate(n, data[sample(nrow(data), size = size, replace = F), , drop = F], simplify = F)
  } else {
    res <- replicate(n, data[sample(nrow(data), replace = TRUE), , drop = F], simplify = F)
  }
  tibble::data_frame(strap = res)
}


#' @title Standard Error and Confidence Intervals for bootstrapped estimates
#' @name boot_ci
#'
#' @description Compute bootstrap standard error and confidence intervals
#'              for a vector of bootstrap replicate estimates.
#'
#' @param x A vector.
#'
#' @return The bootstrap standard error or the lower and upper confidence
#'         intervals of \code{x}.
#'
#' @details This method requires a vector of bootstrap replicate estimates
#'          as input. The function then computes the bootstrap standard error
#'          by calculating the standard deviation of the input vector. The mean
#'          value of the input vector is used to calculate the lower and upper
#'          confidence interval, assuming a t-distribution of bootstrap estimate
#'          replicates.
#'
#' @examples
#' data(efc)
#' bs <- bootstrap(efc, 100)
#'
#' # now run models for each bootstrapped sample
#' bs$models <- lapply(bs$strap, function(x) lm(neg_c_7 ~ e42dep + c161sex, data = x))
#'
#' # extract coefficient "dependency" from each model
#' bs$dependency <- unlist(lapply(bs$models, function(x) coef(x)[2]))
#'
#' # get bootstrapped confidence intervals
#' boot_ci(bs$dependency)
#'
#' # compare with model fit
#' fit <- lm(neg_c_7 ~ e42dep + c161sex, data = efc)
#' confint(fit)[2, ]
#'
#' # compare coefficients
#' mean(bs$dependency)
#' coef(fit)[2]
#'
#' @importFrom stats qt
#' @export
boot_ci <- function(x) {
  # get bootstrap standard error
  bootse <- stats::qt(.975, df = length(x) - 1) * boot_se(x)
  # lower and upper confidence interval
  ci <- mean(x, na.rm = T) + c(-bootse, bootse)
  names(ci) <- c("conf.low", "conf.high")
  ci
}

#' @rdname boot_ci
#' @importFrom stats sd
#' @export
boot_se <- function(x) {
  # compute 1.96 * se for bootstrap replicates
  # see https://www.zoology.ubc.ca/~schluter/R/resample/
  stats::sd(x, na.rm = T)
}
