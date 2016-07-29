#' @title Standard Error and Confidence Intervals for bootstrapped estimates
#' @name boot_ci
#'
#' @description Compute bootstrap standard error, confidence intervals and
#'              p-value for a vector of bootstrap replicate estimates.
#'
#' @param data A data frame that containts the vector with bootstrapped
#'          estimates, or directly the vector (see 'Examples').
#' @param x Name of the variable with bootstrapped estimates. Required, if
#'          \code{data} is a data frame and no vector.
#'
#' @return The bootstrap standard error, the lower and upper confidence
#'         intervals or the p-value of the bootstrapped estimates.
#'
#' @details This method requires a vector of bootstrap replicate estimates
#'          as input. The function then computes the bootstrap standard error
#'          by calculating the standard deviation of the input vector. The mean
#'          value of the input vector is used to calculate the lower and upper
#'          confidence interval, assuming a t-distribution of bootstrap estimate
#'          replicates.
#'
#' @seealso \code{\link{bootstrap}} to genearte bootstrap samples.
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
#' # extract coefficient "gender" from each model
#' bs$gender <- unlist(lapply(bs$models, function(x) coef(x)[3]))
#'
#' # check p-value
#' boot_p(bs$gender)
#' summary(fit)$coefficients[3, ]
#'
#' @importFrom stats qt
#' @export
boot_ci <- function(data, x) {
  # check if data is a data frame
  if (is.data.frame(data)) {
    # evaluate argument
    x <- deparse(substitute(x))
    # get vector
    x <- data[[x]]
  } else {
    x <- data
  }
  # get bootstrap standard error
  bootse <- stats::qt(.975, df = length(x) - 1) * boot_se(data = x)
  # lower and upper confidence interval
  ci <- mean(x, na.rm = T) + c(-bootse, bootse)
  names(ci) <- c("conf.low", "conf.high")
  ci
}


#' @rdname boot_ci
#' @importFrom stats sd
#' @export
boot_se <- function(data, x) {
  # check if data is a data frame
  if (is.data.frame(data)) {
    # evaluate argument
    x <- deparse(substitute(x))
    # get vector
    x <- data[[x]]
  } else {
    x <- data
  }
  # compute SE, see https://www.zoology.ubc.ca/~schluter/R/resample/
  stats::sd(x, na.rm = T)
}


#' @rdname boot_ci
#' @importFrom stats sd pt
#' @export
boot_p <- function(data, x) {
  # check if data is a data frame
  if (is.data.frame(data)) {
    # evaluate argument
    x <- deparse(substitute(x))
    # get vector
    x <- data[[x]]
  } else {
    x <- data
  }
  # compute t-statistic
  t.stat <- mean(x, na.rm = T) / stats::sd(x, na.rm = T)
  # compute p-value
  2 * stats::pt(abs(t.stat), df = length(x) - 1, lower.tail = FALSE)
}
