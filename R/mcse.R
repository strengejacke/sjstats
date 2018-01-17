#' @rdname hdi
#' @export
mcse <- function(x, type = c("fixed", "random", "all")) {
  UseMethod("mcse")
}


#' @export
mcse.brmsfit <- function(x, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)
  mcse_helper(x, type)
}


#' @export
mcse.stanreg <- function(x, type = c("fixed", "random", "all")) {
  # check arguments
  type <- match.arg(type)
  mcse_helper(x, type)
}


#' @importFrom purrr map_dbl
#' @importFrom tibble as_tibble
#' @importFrom dplyr pull
mcse_helper <- function(x, type) {
  dat <- tibble::as_tibble(x)
  if (inherits(x, "brmsfit")) dat <- brms_clean(dat)

  # get standard deviations from posterior samples
  stddev <- purrr::map_dbl(dat, sd)

  # remove certain terms
  keep <-  which(!(names(stddev) %in% c("lp__", "log-posterior", "mean_PPD")))
  stddev <- stddev[keep]

  # get effective sample sizes
  ess <- dplyr::pull(n_eff(x, type = "all"), "n_eff")

  # compute mcse
  dat <- tibble::tibble(
    term = colnames(dat),
    mcse = stddev / sqrt(ess)
  )

  # check if we need to remove random or fixed effects
  remove_effects_from_stan(dat, type, is.brms = inherits(x, "brmsfit"))
}
