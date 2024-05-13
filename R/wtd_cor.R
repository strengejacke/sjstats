#' @rdname weighted_se
#' @export
weighted_correlation <- function(data, ...) {
  UseMethod("weighted_correlation")
}


#' @rdname weighted_se
#' @export
weighted_correlation.default <- function(data, x, y, weights, ci.lvl = 0.95, ...) {
  if (!missing(ci.lvl) && (length(ci.lvl) != 1 || !is.finite(ci.lvl) || ci.lvl < 0 || ci.lvl > 1))
    insight::format_error("'ci.lvl' must be a single number between 0 and 1.")

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
  dat <- suppressMessages(data[vars])
  dat <- stats::na.omit(dat)

  xv <- dat[[x.name]]
  yv <- dat[[y.name]]
  wv <- dat[[w.name]]

  weighted_correlation_helper(xv, yv, wv, ci.lvl)
}


#' @rdname weighted_se
#' @export
weighted_correlation.formula <- function(formula, data, ci.lvl = 0.95, ...) {

  if (!missing(ci.lvl) && (length(ci.lvl) != 1 || !is.finite(ci.lvl) || ci.lvl < 0 || ci.lvl > 1))
    insight::format_error("'ci.lvl' must be a single number between 0 and 1.")

  vars <- all.vars(formula)

  if (length(vars) < 3) {
    vars <- c(vars, "weights")
    data$weights <- 1
  }

  # get data
  dat <- suppressMessages(data[vars])
  dat <- stats::na.omit(dat)

  xv <- dat[[vars[1]]]
  yv <- dat[[vars[2]]]
  wv <- dat[[vars[3]]]

  weighted_correlation_helper(xv, yv, wv, ci.lvl)
}


weighted_correlation_helper <- function(xv, yv, wv, ci.lvl) {

  x <- xv - weighted_mean(xv, weights = wv)
  y <- yv - weighted_mean(yv, weights = wv)

  x <- x / weighted_sd(x, weights = wv)
  y <- y / weighted_sd(y, weights = wv)

  results <- stats::coef(summary(stats::lm(y ~ x, weights = wv)))[2, ]

  ci <- ci.lvl - ((1 - ci.lvl) / 2)
  ci <- results[1] + (stats::qnorm(ci) * c(-1, 1) * results[2])

  structure(
    class = "sj_wcor",
    list(
      estimate = results[1],
      method = "Pearson's Correlation Coefficient",
      p.value = results[4],
      ci = ci,
      ci.lvl = ci.lvl
    )
  )
}
