#' @rdname wtd_sd
#' @export
wtd_chisqtest <- function(data, ...) {
  UseMethod("wtd_chisqtest")
}

#' @importFrom dplyr select
#' @rdname wtd_sd
#' @export
wtd_chisqtest.default <- function(data, x, y, weights, ...) {
  x.name <- deparse(substitute(x))
  y.name <- deparse(substitute(y))
  w.name <- deparse(substitute(weights))

  if (w.name == "NULL") {
    w.name <- "weights"
    data$weights <- 1
  }

  # create string with variable names
  vars <- c(x.name, y.name, w.name)

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  colnames(dat)[3] <- ".weights"
  xtab_statistics(data = dat, statistics = "auto", weights = ".weights", ...)
}


#' @importFrom stats xtabs
#' @rdname wtd_sd
#' @export
wtd_chisqtest.formula <- function(formula, data, ...) {
  vars <- all.vars(formula)

  if (length(vars) < 3) {
    vars <- c(vars, ".weights")
    data$.weights <- 1
  }

  tab <- as.table(round(stats::xtabs(data[[vars[3]]] ~ data[[vars[1]]] + data[[vars[2]]])))
  class(tab) <- "table"
  xtab_statistics(data = tab, statistics = "auto", weights = NULL, ...)
}
