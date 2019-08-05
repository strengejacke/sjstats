#' @rdname wtd_sd
#' @export
wtd_cor <- function(data, ...) {
  UseMethod("wtd_cor")
}


#' @rdname wtd_sd
#' @export
wtd_cor.default <- function(data, x, y, weights, ci.lvl = .95, ...) {
  if (!missing(ci.lvl) & (length(ci.lvl) != 1 || !is.finite(ci.lvl) || ci.lvl < 0 || ci.lvl > 1))
    stop("'ci.lvl' must be a single number between 0 and 1")

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

  xv <- dat[[x.name]]
  yv <- dat[[y.name]]
  wv <- dat[[w.name]]

  wtd_cor_helper(xv, yv, wv, ci.lvl)
}


#' @rdname wtd_sd
#' @export
wtd_cor.formula <- function(formula, data, ci.lvl = .95, ...) {

  if (!missing(ci.lvl) & (length(ci.lvl) != 1 || !is.finite(ci.lvl) || ci.lvl < 0 || ci.lvl > 1))
    stop("'ci.lvl' must be a single number between 0 and 1")

  vars <- all.vars(formula)

  if (length(vars) < 3) {
    vars <- c(vars, "weights")
    data$weights <- 1
  }

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  xv <- dat[[vars[1]]]
  yv <- dat[[vars[2]]]
  wv <- dat[[vars[3]]]

  wtd_cor_helper(xv, yv, wv, ci.lvl)
}


#' @importFrom stats cor.test
wtd_cor_helper <- function(xv, yv, wv, ci.lvl) {

  x <- xv - wtd_mean(xv, weights = wv)
  y <- yv - wtd_mean(yv, weights = wv)

  x <- x / wtd_sd(x, weights = wv)
  y <- y / wtd_sd(y, weights = wv)

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
