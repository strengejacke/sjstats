#' @title Measures of associations for contingency tables
#' @name phi
#' @description Compute measures of associtations for a contingency table,
#'              like \emph{Phi coefficient} or \emph{Cramer's V}.
#'
#' @param tab A \code{\link{table}} or \code{\link{ftable}}. Tables of class
#'          \code{\link{xtabs}} and other will be coerced to \code{\link{ftable}}
#'          objects.
#'
#' @return \itemize{
#'           \item For \code{phi()}, the table's Phi value.
#'           \item For \code{cramer()}, the table's Cramer's V.
#'         }
#'
#'
#' @examples
#' # Phi coefficient for 2x2 tables
#' tab <- table(sample(1:2, 30, TRUE), sample(1:2, 30, TRUE))
#' phi(tab)
#'
#' # Cramer's V for nominal variables with more than 2 categories
#' tab <- table(sample(1:2, 30, TRUE), sample(1:3, 30, TRUE))
#' cramer(tab)
#'
#' @importFrom MASS loglm
#' @importFrom stats ftable
#' @export
phi <- function(tab) {
  # convert to flat table
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)
  tb <- summary(MASS::loglm(~1 + 2, tab))$tests
  phi_val <- sqrt(tb[2, 1] / sum(tab))
  return(phi_val)
}


#' @rdname phi
#' @importFrom stats ftable
#' @export
cramer <- function(tab) {
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)
  phi_val <- phi(tab)
  cramer <- sqrt(phi_val ^ 2 / min(dim(tab) - 1))
  return(cramer)
}
