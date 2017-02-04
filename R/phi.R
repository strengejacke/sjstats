#' @rdname xtab_statistics
#' @export
phi <- function(tab) {
  # convert to flat table
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)
  tb <- summary(MASS::loglm(~1 + 2, tab))$tests
  phi_val <- sqrt(tb[2, 1] / sum(tab))
  return(phi_val)
}


#' @rdname xtab_statistics
#' @export
cramer <- function(tab) {
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)
  phi_val <- phi(tab)
  cramer <- sqrt(phi_val ^ 2 / min(dim(tab) - 1))
  return(cramer)
}
