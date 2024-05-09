#' @rdname crosstable_statistics
#' @export
cramers_v <- function(tab, ...) {
  UseMethod("cramers_v")
}

#' @rdname crosstable_statistics
#' @export
cramer <- cramers_v

#' @export
cramers_v.table <- function(tab, ...) {
  .cramers_v(tab)
}


#' @export
cramers_v.ftable <- function(tab, ...) {
  .cramers_v(tab)
}


#' @rdname crosstable_statistics
#' @export
cramers_v.formula <- function(formula, data, ci.lvl = NULL, n = 1000, method = c("dist", "quantile"), ...) {
  terms <- all.vars(formula)
  tab <- table(data[[terms[1]]], data[[terms[2]]])
  method <- match.arg(method)

  if (is.null(ci.lvl) || is.na(ci.lvl)) {
    .cramers_v(tab)
  } else {
    straps <- sjstats::bootstrap(data[terms], n)
    tables <- lapply(straps$strap, function(x) {
      dat <- as.data.frame(x)
      table(dat[[1]], dat[[2]])
    })
    cramers <- sapply(tables, function(x) .cramers_v(x))
    ci <- boot_ci(cramers, ci.lvl = ci.lvl, method = method)

    data_frame(
      cramer = .cramers_v(tab),
      conf.low = ci$conf.low,
      conf.high = ci$conf.high
    )
  }
}

.cramers_v <- function(tab) {
  # convert to flat table
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)
  sqrt(phi(tab)^2 / min(dim(tab) - 1))
}
