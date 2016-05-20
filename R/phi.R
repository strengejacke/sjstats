#' @title Phi value for contingency tables
#' @name phi
#' @description Compute Phi value for a contingency table.
#'
#' @seealso \code{\link{cramer}}

#' @param tab A \code{\link{table}} or \code{\link{ftable}}. Tables of class
#'          \code{\link{xtabs}} and other will be coerced to \code{\link{ftable}}
#'          objects.
#'
#' @return The table's Phi value.
#'
#' @examples
#' tab <- table(sample(1:2, 30, TRUE), sample(1:2, 30, TRUE))
#' phi(tab)
#'
#' @importFrom MASS loglm
#' @importFrom stats ftable
#' @export
phi <- function(tab) {
  # convert to flat table
  if (all(class(tab) != "ftable")) tab <- stats::ftable(tab)
  tb <- summary(MASS::loglm(~1 + 2, tab))$tests
  phi_val <- sqrt(tb[2, 1] / sum(tab))
  return(phi_val)
}
