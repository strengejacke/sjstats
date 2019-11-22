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
        phis = sapply(.data$tables, function(x) .cramer(x))
      ) %>%
      dplyr::pull("phis") %>%
      boot_ci(ci.lvl = ci.lvl, method = method)

    data_frame(
      phi = .phi(tab),
      conf.low = ci$conf.low,
      conf.high = ci$conf.high
    )
  }
}


.phi <- function(tab) {
  # convert to flat table
  if (!inherits(tab, "ftable")) tab <- stats::ftable(tab)

  tb <- summary(MASS::loglm(~1 + 2, tab))$tests
  sqrt(tb[2, 1] / sum(tab))
}
