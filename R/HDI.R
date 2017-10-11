#' @title Compute high density intervals (HDI) for MCMC samples
#' @name hdi
#'
#' @description \code{hdi()} computes the high density interval for values from
#'              MCMC samples. \code{rope()} calculates the proportion of a posterior
#'              distribution that lies within a region of practical equivalence.
#'
#' @param x A vector of values from a probability distribution (e.g., posterior
#'        probabilities from MCMC sampling), or a \code{stanreg},
#'        \code{stanfit}, or \code{brmsfit} object.
#' @param prob Scalar between 0 and 1, indicating the mass within the credible
#'        interval that is to be estimated.
#' @param rope Vector of length two, indicating the lower and upper limit of a
#'        range around zero, which indicates the region of practical equivalence.
#'        Values of the posterior distribution within this range are considered as
#'        being "practically equivalent to zero".
#' @param trans Name of a function or character vector naming a function, used
#'        to apply transformations on the returned HDI-values resp.
#'        (for \code{rope()}) on the values of the posterior distribution, before
#'        calculating the rope based on the boundaries given in \code{rope}. Note
#'        that the values in \code{rope} are not transformed.
#' @param type For mixed effects models, specify the type of effects that should
#'        be returned. \code{type = "fixed"} returns fixed effects only,
#'        \code{type = "random"} the random effects and \code{type = "all"} returns
#'        both fixed and random effects.
#'
#'
#' @return For \code{hdi()}, if \code{x} is a vector, returns a vector of length two
#'         with the lower and upper limit of the HDI; if \code{x} is a
#'         \code{stanreg}, \code{stanfit} or \code{brmsfit} object, returns a
#'         tibble with lower and upper HDI-limits for each predictor.
#'         For \code{rope()}, returns the proportion of values from \code{x}
#'         that are within the boundaries of \code{rope}.
#'
#' @details Computation for HDI is based on the code from Kruschke 2015, pp. 727f.
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


#' @rdname hdi
#' @export
rope <- function(x, rope, trans = NULL, type = c("fixed", "random", "all")) {
  UseMethod("rope")
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


#' @export
rope.default <- function(x, rope, trans = NULL, type = c("fixed", "random", "all")) {
  rope_helper(x, rope, trans)
}


#' @export
rope.stanreg <- function(x, rope, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  # get posterior data
  dat <- x %>%
    as.data.frame() %>%
    purrr::map_df(~ rope_helper(.x, rope, trans)) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column()

  colnames(dat) <- c("term", "rope")

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = FALSE)
}


#' @export
rope.brmsfit <- function(x, rope, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  # get posterior data
  dat <- x %>%
    tibble::as_tibble() %>%
    purrr::map_df(~ rope_helper(.x, rope, trans)) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column()

  colnames(dat) <- c("term", "rope")

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = TRUE)
}


#' @export
rope.stanfit <- function(x, rope, trans = NULL, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  # get posterior data
  dat <- x %>%
    as.data.frame() %>%
    purrr::map_df(~ rope_helper(.x, rope, trans)) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column()

  colnames(dat) <- c("term", "rope")

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = FALSE)
}


#' @importFrom dplyr between
rope_helper <- function(x, rope, trans) {
  # stop if argument is not correct
  if (length(rope) != 2)
    stop("Argument `rope` needs to be a vector of length two.", call. = F)

  # switch values, if lower bound is larger than upper bound
  if (rope[1] > rope[2]) {
    tmp <- rope[2]
    rope[2] <- rope[1]
    rope[1] <- tmp
  }

  # check if we have correct function
  if (!is.null(trans)) {
    trans <- match.fun(trans)
    x <- trans(x)
  }

  # sort values, to compute rope
  x <- sort(x)
  r <- dplyr::between(x, rope[1], rope[2])

  # compute proportion of values within boundaries
  sum(r) / length(x)
}
