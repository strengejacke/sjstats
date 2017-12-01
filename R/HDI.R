#' @title Compute statistics for MCMC samples
#' @name hdi
#'
#' @description \code{hdi()} computes the high density interval for values from
#'   MCMC samples. \code{rope()} calculates the proportion of a posterior
#'   distribution that lies within a region of practical equivalence.
#'   \code{n_eff()} calculates the number of effective samples (effective
#'   sample size). \code{mcse()} returns the Monte Carlo standard error.
#'
#' @param x A \code{stanreg}, \code{stanfit}, or \code{brmsfit} object. For
#'   \code{hdi()} and \code{rope()}, may also be a vector of values from a
#'   probability distribution (e.g., posterior probabilities from MCMC sampling).
#' @param prob Scalar between 0 and 1, indicating the mass within the credible
#'   interval that is to be estimated.
#' @param rope Vector of length two, indicating the lower and upper limit of a
#'   range around zero, which indicates the region of practical equivalence.
#'   Values of the posterior distribution within this range are considered as
#'   being "practically equivalent to zero".
#' @param trans Name of a function or character vector naming a function, used
#'   to apply transformations on the returned HDI-values resp.
#'   (for \code{rope()}) on the values of the posterior distribution, before
#'   calculating the rope based on the boundaries given in \code{rope}. Note
#'   that the values in \code{rope} are not transformed.
#' @param type For mixed effects models, specify the type of effects that should
#'   be returned. \code{type = "fixed"} returns fixed effects only,
#'   \code{type = "random"} the random effects and \code{type = "all"} returns
#'   both fixed and random effects.
#'
#'
#' @return For \code{hdi()}, if \code{x} is a vector, returns a vector of length
#'   two with the lower and upper limit of the HDI; if \code{x} is a
#'   \code{stanreg}, \code{stanfit} or \code{brmsfit} object, returns a
#'   tibble with lower and upper HDI-limits for each predictor.
#'   For \code{rope()}, returns the proportion of values from \code{x}
#'   that are within the boundaries of \code{rope}. \code{mcse()} and
#'   \code{n_eff()} return a tibble with two columns: one with the term names
#'   and one with the related statistic.
#'
#' @details Computation for HDI is based on the code from Kruschke 2015, pp. 727f.
#'   For default sampling in Stan (4000 samples), the 90\% intervals for HDI are
#'   more stable than, for instance, 95\% intervals. An effective sample size
#'   (see \code{n_eff()}) of at least 10.000 is recommended if 95\% intervals
#'   should be computed (see Kruschke 2015, p. 183ff). \cr \cr
#'   The formula for \code{mcse()} is from Kruschke 2015, p. 187.
#'
#' @references Kruschke JK. Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. 2nd edition. Academic Press, 2015
#'
#' @examples
#' \dontrun{
#' if (require("rstanarm")) {
#'   fit <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#'   hdi(fit)
#'
#'   # fit logistic regression model
#'   fit <- stan_glm(
#'     vs ~ wt + am,
#'     data = mtcars,
#'     family = binomial("logit"),
#'     chains = 1
#'   )
#'   # compute hdi, transform on "odds ratio scale"
#'   hdi(fit, trans = exp)
#'
#'   # compute rope, on scale of linear predictor. finds proportion
#'   # of posterior distribution values between -1 and 1.
#'   rope(fit, rope = c(-1, 1))
#'
#'   # compute rope, boundaries as "odds ratios". finds proportion of
#'   # posterior distribution values, which - after being exponentiated -
#'   # are between .8 and 1.25 (about -.22 and .22 on linear scale)
#'   rope(fit, rope = c(.8, 1.25), trans = exp)
#' }}
#'
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom purrr map_dbl map_df
#' @importFrom sjmisc rotate_df
#' @export
hdi <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  UseMethod("hdi")
}


#' @export
hdi.stanreg <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  # get posterior data
  dat <- x %>%
    tibble::as_tibble() %>%
    purrr::map_df(~ hdi_helper(.x, prob, trans)) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column()

  colnames(dat) <- c("term", "hdi.low", "hdi.high")

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = FALSE)
}


#' @export
hdi.brmsfit <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  # get posterior data
  dat <- x %>%
    tibble::as_tibble() %>%
    purrr::map_df(~ hdi_helper(.x, prob, trans)) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column()

  colnames(dat) <- c("term", "hdi.low", "hdi.high")

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = TRUE)
}


#' @export
hdi.stanfit <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  # get posterior data
  dat <- x %>%
    as.data.frame() %>%
    purrr::map_df(~ hdi_helper(.x, prob, trans)) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column()

  colnames(dat) <- c("term", "hdi.low", "hdi.high")

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = FALSE)
}


#' @export
hdi.default <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
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
