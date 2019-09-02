#' @rdname xtab_statistics
#' @export
cramer <- function(tab, ...) {
  UseMethod("cramer")
}


#' @export
cramer.table <- function(tab, ...) {
  .cramer(tab)
}


#' @export
cramer.ftable <- function(tab, ...) {
  .cramer(tab)
}


#' @rdname xtab_statistics
#' @export
cramer.formula <- function(formula, data, ci.lvl = NULL, n = 1000, method = c("dist", "quantile"), ...) {
  terms <- all.vars(formula)
  tab <- table(data[[terms[1]]], data[[terms[2]]])
  method <- match.arg(method)

  if (is.null(ci.lvl) || is.na(ci.lvl)) {
    .cramer(tab)
  } else {
    ci <- data[, terms] %>%
      sjstats::bootstrap(n) %>%
      dplyr::mutate(
        tables = lapply(.data$strap, function(x) {
          dat <- as.data.frame(x)
          table(dat[[1]], dat[[2]])
        }),
        cramers = sapply(.data$tables, function(x) .cramer(x))
      ) %>%
      dplyr::pull("cramers") %>%
      boot_ci(ci.lvl = ci.lvl, method = method)

    data_frame(
      cramer = .cramer(tab),
      conf.low = ci$conf.low,
      conf.high = ci$conf.high
    )
  }
}

.cramer <- function(tab) {
  # convert to flat table
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)
  sqrt(phi(tab)^2 / min(dim(tab) - 1))
}
