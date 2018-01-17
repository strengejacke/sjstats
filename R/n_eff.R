#' @rdname hdi
#' @export
n_eff <- function(x, type = c("fixed", "random", "all")) {
  UseMethod("n_eff")
}


#' @export
n_eff.default <- function(x, type = c("fixed", "random", "all")) {
  NextMethod()
}


#' @export
n_eff.stanreg <- function(x, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)
  n_eff_helper(x, smry = summary(x)[, "n_eff"], type)
}


#' @export
n_eff.brmsfit <- function(x, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)

  if (!requireNamespace("rstan")) stop("Package `rstan` required.", call. = F)
  n_eff_helper(x, smry = rstan::summary(x$fit)$summary[, "n_eff"], type)
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
