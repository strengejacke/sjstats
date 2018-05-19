#' @rdname hdi
#' @export
n_eff <- function(x, type = c("fixed", "random", "all")) {
  UseMethod("n_eff")
}


#' @export
n_eff.default <- function(x, type = c("fixed", "random", "all")) {
  NextMethod()
}


#' @importFrom bayesplot neff_ratio
#' @export
n_eff.stanreg <- function(x, type = c("fixed", "random", "all")) {
  type <- match.arg(type)
  n_eff_helper(x, smry = bayesplot::neff_ratio(x), type)
}


#' @importFrom bayesplot neff_ratio
#' @export
n_eff.stanfit <- function(x, type = c("fixed", "random", "all")) {
  type <- match.arg(type)
  n_eff_helper(x, smry = bayesplot::neff_ratio(x), type)
}


#' @export
n_eff.brmsfit <- function(x, type = c("fixed", "random", "all")) {
  type <- match.arg(type)

  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  n_eff_helper(x, smry = brms::neff_ratio(x), type)
}


#' @importFrom tibble as_tibble
#' @importFrom dplyr slice
n_eff_helper <- function(x, smry, type) {
  if (inherits(x, "brmsfit"))
    tn <- colnames(tibble::as_tibble(x))
  else
    tn <- names(smry)

  dat <- tibble::tibble(
    term = tn,
    n_eff = as.vector(smry)
  )

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = inherits(x, "brmsfit"))
}
