#' @title Compute statistics for MCMC samples and Stan models
#' @name mcse
#'
#' @description \code{n_eff()} calculates the the number of effective samples
#'   (effective sample size). \code{mcse()} returns the Monte Carlo standard
#'   error.
#'
#' @param x A \code{stanreg}, \code{stanfit}, or \code{brmsfit} object.
#' @param type For mixed effects models, specify the type of effects that should
#'   be returned. \code{type = "fixed"} returns fixed effects only,
#'   \code{type = "random"} the random effects and \code{type = "all"} returns
#'   both fixed and random effects.
#' @param ... Not used
#'
#'
#' @return \code{mcse()} and \code{n_eff()} return a tibble with two columns: one
#'   with the term names and one with the related statistic resp. effective
#'   sample size.
#'
#' @details
#'   \describe{
#'     \item{\strong{MCSE}}{
#'       The Monte Carlo Standard Error is another useful measure of accuracy of
#'       the chains. It is defined as standard deviation of the chains divided by
#'       their effective sample size (the formula for \code{mcse()} is from
#'       Kruschke 2015, p. 187). The MCSE \dQuote{provides a quantitative suggestion
#'       of how big the estimation noise is}.
#'     }
#'     \item{\strong{Number of Effective Samples}}{
#'       The effective sample size divides the actual sample size by the amount
#'       of autocorrelation. The effective sample size is a measure of \dQuote{how
#'       much independent information there is in autocorrelated chains}, or:
#'       \dQuote{What would be the sample size of a completely non-autocorrelated chain
#'       that yielded the same information?} (\emph{Kruschke 2015, p182-3}).
#'       The ratio of effective number of samples and total number of samples
#'       (provided in \code{tidy_stan()}) ranges from 0 to 1, and should be close
#'       to 1. The closer this ratio comes to zero means that the chains may be
#'       inefficient, but possibly still okay.
#'     }
#'   }
#'
#' @references Kruschke JK. Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. 2nd edition. Academic Press, 2015
#'
#' @examples
#' \dontrun{
#' if (require("rstanarm")) {
#'   m <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#'   mcse(m)
#'   n_eff(m)
#' }}
#' @export
mcse <- function(x, ...) {
  UseMethod("mcse")
}


#' @export
mcse.brmsfit <- function(x, type = c("fixed", "random", "all"), ...) {
  # check arguments
  type <- match.arg(type)
  mcse_helper(x, type)
}


#' @rdname mcse
#' @export
mcse.stanmvreg <- function(x, type = c("fixed", "random", "all"), ...) {
  # check arguments
  type <- match.arg(type)

  s <- summary(x)
  dat <- data_frame(
    term = rownames(s),
    mcse = s[, "mcse"]
  )

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = FALSE)
}


#' @rdname mcse
#' @export
mcse.stanreg <- function(x, type = c("fixed", "random", "all"), ...) {
  # check arguments
  type <- match.arg(type)

  if (inherits(x, "stanmvreg"))
    return(mcse.stanmvreg(x = x, type = type, ...))

  mcse_helper(x, type)
}


#' @importFrom purrr map_dbl
#' @importFrom dplyr pull
mcse_helper <- function(x, type) {
  dat <- as.data.frame(x)
  if (inherits(x, "brmsfit")) dat <- brms_clean(dat)

  # get standard deviations from posterior samples
  stddev <- purrr::map_dbl(dat, sd)

  # remove certain terms
  keep <- seq_len(length(names(stddev)))

  for (.x in c("^(?!lp__)", "^(?!log-posterior)", "^(?!mean_PPD)")) {
    keep <- intersect(keep, grep(.x, names(stddev), perl = TRUE))
  }

  stddev <- stddev[keep]


  # get effective sample sizes
  ess <- dplyr::pull(n_eff(x, type = "all"), "n_eff")

  # compute mcse
  dat <- data_frame(
    term = colnames(dat),
    mcse = stddev / sqrt(ess)
  )

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = inherits(x, "brmsfit"))
}
