#' @title Cramer's V for a contingency table
#' @name cramer
#' @description Compute Cramer's V for a table with more than 2x2 fields.
#'
#' @seealso \code{\link{phi}}
#'
#' @inheritParams phi
#'
#' @return The table's Cramer's V.
#'
#' @examples
#' tab <- table(sample(1:2, 30, TRUE), sample(1:3, 30, TRUE))
#' cramer(tab)
#'
#' @importFrom stats ftable
#' @export
cramer <- function(tab) {
  if (all(class(tab) != "ftable")) tab <- stats::ftable(tab)
  phi_val <- phi(tab)
  cramer <- sqrt(phi_val ^ 2 / min(dim(tab) - 1))
  return(cramer)
}
