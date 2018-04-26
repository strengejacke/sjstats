#' @title Compute statistics for MCMC samples
#' @name hdi
#'
#' @description \code{hdi()} computes the highest density interval for values from
#'   MCMC samples. \code{rope()} calculates the proportion of a posterior
#'   distribution that lies within a region of practical equivalence.
#'   \code{n_eff()} calculates the number of effective samples (effective
#'   sample size). \code{mcse()} returns the Monte Carlo standard error.
#'
#' @param x A \code{stanreg}, \code{stanfit}, or \code{brmsfit} object. For
#'   \code{hdi()} and \code{rope()}, may also be a data frame or a vector
#'   of values from a probability distribution (e.g., posterior probabilities
#'   from MCMC sampling).
#' @param prob Vector of scalars between 0 and 1, indicating the mass within
#'   the credible interval that is to be estimated. See \code{\link{hdi}}.
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
#'   tibble with lower and upper HDI-limits for each predictor. To distinguish
#'   multiple HDI values, column names for the HDI get a suffix when \code{prob}
#'   has more than one element.
#'   \cr \cr
#'   For \code{rope()}, returns a tibble with two columns: the proportion of
#'   values from \code{x} that are within and outside the boundaries of
#'   \code{rope}.
#'   \cr \cr
#'   \code{mcse()} and \code{n_eff()} return a tibble with two columns: one
#'   with the term names and one with the related statistic.
#'
#' @details Computation for HDI is based on the code from Kruschke 2015, pp. 727f.
#'   For default sampling in Stan (4000 samples), the 90\% intervals for HDI are
#'   more stable than, for instance, 95\% intervals. An effective sample size
#'   (see \code{n_eff()}) of at least 10.000 is recommended if 95\% intervals
#'   should be computed (see Kruschke 2015, p. 183ff). \cr \cr
#'   The formula for \code{mcse()} is from Kruschke 2015, p. 187.
#'   \cr \cr
#'   For \code{n_eff()}, the ratio of effective number of samples ranges from 0
#'   to 1, and should be close to 1. The closer this ratio comes to zero means
#'   that the chains may be inefficient, but possibly still okay.
#'
#' @references Kruschke JK. Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. 2nd edition. Academic Press, 2015
#'
#' @examples
#' \dontrun{
#' if (require("rstanarm")) {
#'   fit <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#'   hdi(fit)
#'
#'   # return multiple intervals
#'   hdi(fit, prob = c(.5, .7, .9))
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
#' @export
hdi <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  UseMethod("hdi")
}


#' @export
hdi.stanreg <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  dat <- hdi_worker(x = x, prob = prob, trans = trans)

  # check if we need to remove random or fixed effects
  dat <- remove_effects_from_stan(dat, type, is.brms = FALSE)

  class(dat) <- c("sj_hdi", class(dat))
  dat
}


#' @export
hdi.brmsfit <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  dat <- hdi_worker(x = x, prob = prob, trans = trans)

  # check if we need to remove random or fixed effects
  dat <- remove_effects_from_stan(dat, type, is.brms = TRUE)

  class(dat) <- c("sj_hdi", class(dat))
  dat
}


#' @export
hdi.stanfit <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  dat <- hdi_worker(x = x, prob = prob, trans = trans)

  # check if we need to remove random or fixed effects
  dat <- remove_effects_from_stan(dat, type, is.brms = FALSE)

  class(dat) <- c("sj_hdi", class(dat))
  dat
}


#' @export
hdi.data.frame <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  dat <- hdi_worker(x = x, prob = prob, trans = trans)

  class(dat) <- c("sj_hdi", class(dat))
  dat
}


#' @export
hdi.default <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all")) {
  hdi_helper(x, prob, trans)
}


#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom purrr map_df
#' @importFrom sjmisc rotate_df
hdi_worker <- function(x, prob, trans) {
  dat <- purrr::map(
    prob,
    function(i) {

      out <- x %>%
        tibble::as_tibble() %>%
        purrr::map_df(~ hdi_helper(.x, i, trans)) %>%
        sjmisc::rotate_df() %>%
        tibble::rownames_to_column()

      colnames(out) <- c("term", "hdi.low", "hdi.high")
      out

    }) %>%
    dplyr::bind_cols() %>%
    dplyr::select(1, tidyselect::starts_with("hdi."))



  # for multiple HDIs, fix column names

  if (length(prob) > 1) {
    suffix <- prob %>%
      purrr::map(~ rep(.x, 2)) %>%
      purrr::flatten_dbl()

    colnames(dat)[2:ncol(dat)] <-
      sprintf(
        "%s_%s",
        rep(c("hdi.low", "hdi.high"), length(prob)),
        as.character(suffix)
      )
  }

  dat
}


# based on Kruschke 2015, pp727f
#' @importFrom purrr map_dbl map_df
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
