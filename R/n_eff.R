#' @rdname hdi
#' @export
n_eff <- function(x, ...) {
  UseMethod("n_eff")
}


#' @export
n_eff.default <- function(x, ...) {
  NextMethod()
}


#' @rdname hdi
#' @export
n_eff.stanreg <- function(x, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  s <- summary(x)
  n_eff_helper(rownames(s), s[, "n_eff"], type, FALSE)
}


#' @export
n_eff.stanfit <- function(x, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  s <- summary(x)
  n_eff_helper(rownames(s), s[, "n_eff"], type, FALSE)
}


#' @rdname hdi
#' @export
n_eff.brmsfit <- function(x, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)

  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Please install and load package `rstan` first.")

  s <- rstan::summary(x$fit)$summary
  n_eff_helper(rownames(s), s[, "n_eff"], type = type, is.brms = TRUE)
}


#' @importFrom dplyr slice
n_eff_helper <- function(tn, effs, type, is.brms) {
  dat <- data_frame(
    term = tn,
    n_eff = effs
  )

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = is.brms)
}
