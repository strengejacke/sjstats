#' @rdname xtab_statistics
#' @export
phi <- function(tab) {
  # convert to flat table
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)

  tb <- summary(MASS::loglm(~1 + 2, tab))$tests
  sqrt(tb[2, 1] / sum(tab))
}


#' @rdname xtab_statistics
#' @export
cramer <- function(tab, ...) {
  UseMethod("cramer")
}


#' @export
cramer.table <- function(tab, ...) {
  .cramer(tab)
}

#' @rdname xtab_statistics
#' @export
cramer.formula <- function(formula, data, ...) {
  terms <- all.vars(formula)
  tab <- table(data[[terms[1]]], data[[terms[2]]])
  .cramer(tab)
}

.cramer <- function(tab) {
  # convert to flat table
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)
  sqrt(phi(tab)^2 / min(dim(tab) - 1))
}
