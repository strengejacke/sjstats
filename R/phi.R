#' @rdname crosstable_statistics
#' @export
phi <- function(tab, ...) {
  UseMethod("phi")
}


#' @export
phi.table <- function(tab, ...) {
  .phi(tab)
}


#' @export
phi.ftable <- function(tab, ...) {
  .phi(tab)
}


#' @export
phi.formula <- function(formula, data, ci.lvl = NULL, n = 1000, method = c("dist", "quantile"), ...) {
  formula_terms <- all.vars(formula)
  tab <- table(data[[formula_terms[1]]], data[[formula_terms[2]]])
  method <- match.arg(method)

  if (is.null(ci.lvl) || is.na(ci.lvl)) {
    .cramers_v(tab)
  } else {
    straps <- sjstats::bootstrap(data[formula_terms], n)
    tables <- lapply(straps$strap, function(x) {
      dat <- as.data.frame(x)
      table(dat[[1]], dat[[2]])
    })
    phis <- sapply(tables, .phi)
    ci <- boot_ci(phis, ci.lvl = ci.lvl, method = method)

    data_frame(
      phi = .phi(tab),
      conf.low = ci$conf.low,
      conf.high = ci$conf.high
    )
  }
}


.phi <- function(tab) {
  insight::check_if_installed("MASS")
  # convert to flat table
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)

  tb <- summary(MASS::loglm(~1 + 2, tab))$tests
  sqrt(tb[2, 1] / sum(tab))
}
