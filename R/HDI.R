#' @title Compute high density intervals (HDI) for MCMC samples
#' @name hdi
#'
#' @description \code{hdi()} computes the high density interval for values from
#'              MCMC samples.
#'
#' @param x A vector of values from a probability distribution (e.g., posterior
#'        probabilities from MCMC sampling), or a \code{stanreg}-object.
#' @param prob Scalar between 0 and 1, indicating the mass within the credible
#'        interval that is to be estimated.
#' @param trans Name of a function or character vector naming a function, used
#'        to apply transformations on the returned HDI-values.
#'
#'
#' @return If \code{x} is a vector, \code{hdi()} returns a vector of length two
#'         with the lower and upper limit of the HDI; if \code{x} is a
#'         \code{stanreg}-object, returns a tibble with lower and upper HDI-limits
#'         for each predictor.
#'
#' @details Computation is based on the code from Kruschke 2015, pp. 727f.
#'
#' @references Kruschke JK. Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. 2nd edition. Academic Press, 2015
#'
#' @examples
#' library(rstanarm)
#' fit <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#' hdi(fit)
#'
#' # fit logistic regression model
#' fit <- stan_glm(
#'   vs ~ wt + am,
#'   data = mtcars,
#'   family = binomial("logit"),
#'   chains = 1
#' )
#' # compute hdi, transform on "odds ratio scale"
#' hdi(fit, trans = exp)
#'
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom purrr map_dbl map_df
#' @importFrom sjmisc rotate_df
#' @export
hdi <- function(x, prob = .9, trans = NULL) {
  UseMethod("hdi")
}


#' @export
hdi.stanreg <- function(x, prob = .9, trans = NULL) {
  # get posterior data
  dat <- x %>%
    tibble::as_tibble() %>%
    purrr::map_df(~ hdi_helper(.x, prob, trans)) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column()

  colnames(dat) <- c("term", "hdi.low", "hdi.high")

  dat
}


#' @export
hdi.default <- function(x, prob = .9, trans = NULL) {
  hdi_helper(x, prob, trans)
}


# based on Kruschke 2015, pp727f
hdi_helper <- function(x, prob, trans) {
  x <- sort(x)
  ci.index <- ceiling(prob * length(x))
  nCIs <- length(x) - ci.index
  ci.width <- purrr::map_dbl(1:nCIs, ~ x[.x + ci.index] - x[.x])
  HDImin <- x[which.min(ci.width)]
  HDImax <- x[which.min(ci.width) + ci.index]

  # check if we have correct function
  if (!is.null(trans)) {
    trans <- match.fun(trans)
    HDImin <- trans(HDImin)
    HDImax <- trans(HDImax)
  }

  c(HDImin, HDImax)
}
