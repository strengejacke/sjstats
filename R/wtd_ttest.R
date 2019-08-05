#' @rdname wtd_sd
#' @importFrom stats pt qt weighted.mean setNames
#' @importFrom sjmisc is_empty
#' @export
wtd_ttest <- function(data, ...) {
  UseMethod("wtd_ttest")
}

#' @rdname wtd_sd
#' @export
wtd_ttest.default <- function(data, x, y = NULL, weights, mu = 0, paired = FALSE, ci.lvl = 0.95, alternative = c("two.sided", "less", "greater"), ...) {

  if (!missing(ci.lvl) & (length(ci.lvl) != 1 || !is.finite(ci.lvl) || ci.lvl < 0 || ci.lvl > 1))
    stop("'ci.lvl' must be a single number between 0 and 1")

  alternative <- match.arg(alternative)

  x.name <- deparse(substitute(x))
  y.name <- deparse(substitute(y))
  w.name <- deparse(substitute(weights))

  if (y.name == "NULL") y.name <- NULL

  if (w.name == "NULL") {
    w.name <- "weights"
    data$weights <- 1
  }

  # create string with variable names
  vars <- c(x.name, y.name, w.name)

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  if (sjmisc::is_empty(dat) || nrow(dat) == 1) {
    warning("Too less data to compute t-test.")
    return(NULL)
  }

  xv <- dat[[x.name]]
  wx <- wy <- dat[[w.name]]

  if (!is.null(y.name))
    yv <- dat[[y.name]]
  else
    yv <- NULL

  nx <- ny <- nrow(dat)

  wtd_ttest_helper(xv, yv, wx, wy, nx, ny, mu, paired, alternative, ci.lvl, x.name, y.name, NULL)
}


#' @rdname wtd_sd
#' @export
wtd_ttest.formula <- function(formula, data, mu = 0, paired = FALSE, ci.lvl = 0.95, alternative = c("two.sided", "less", "greater"), ...) {

  if (!missing(ci.lvl) & (length(ci.lvl) != 1 || !is.finite(ci.lvl) || ci.lvl < 0 || ci.lvl > 1))
    stop("'ci.lvl' must be a single number between 0 and 1")

  alternative <- match.arg(alternative)

  vars <- all.vars(formula)

  g <- data[[vars[2]]]

  if (is.factor(g))
    grps <- levels(g)
  else
    grps <- na.omit(sort(unique(g)))

  if (length(grps) > 2)
    stop("Grouping factor has more than two levels.")

  if (length(vars) < 3) {
    vars <- c(vars, "weights")
    data$weights <- 1
  }

  x <- data[[vars[1]]]
  y <- data[[vars[2]]]
  w <- data[[vars[3]]]

  xv <- x[y == grps[1]]
  yv <- x[y == grps[2]]
  wx <- w[y == grps[1]]
  wy <- w[y == grps[2]]

  mxv <- is.na(xv)
  xv <- xv[!mxv]
  wx <- wx[!mxv]

  myv <- is.na(yv)
  yv <- yv[!myv]
  wy <- wy[!myv]

  nx <- length(xv)
  ny <- length(yv)

  labs <- sjlabelled::get_labels(
    data[[vars[2]]],
    attr.only = FALSE,
    values = "p",
    drop.na = TRUE,
    drop.unused = TRUE
  )

  wtd_ttest_helper(xv, yv, wx, wy, nx, ny, mu, paired, alternative, ci.lvl, vars[1], vars[2], labs)
}


wtd_ttest_helper <- function(xv, yv, wx, wy, nx, ny, mu, paired, alternative, ci.lvl, x.name, y.name, group.name) {
  if (paired) {
    xv <- xv - yv
    yv <- NULL
  }

  mu.x.w <- stats::weighted.mean(xv, wx)
  var.x.w <- wtd_sd(xv, wx)^2
  se.x <- sqrt(var.x.w / nx)

  if (!is.null(yv)) {
    mu.y.w <- stats::weighted.mean(yv, wy)
    var.y.w <- wtd_sd(yv, wy)^2
    se.y <- sqrt(var.y.w / ny)

    se <- sqrt(se.x^2 + se.y^2)
    df <- se^4 / (se.x^4 / (nx - 1) + se.y^4 / (ny - 1))
    tstat <- (mu.x.w - mu.y.w - mu) / se

    estimate <- c(mu.x.w, mu.y.w)
    names(estimate) <- c("mean of x", "mean of y")
    method <- "Two-Sample t-test"
  } else {
    se <- se.x
    df <- nx - 1
    tstat <- (mu.x.w - mu) / se

    estimate <- stats::setNames(mu.x.w, if (paired)  "mean of the differences" else "mean of x")
    method <- if (paired) "Paired t-test" else "One Sample t-test"
  }



  if (alternative == "less") {
    pval <- stats::pt(tstat, df)
    cint <- c(-Inf, tstat + stats::qt(ci.lvl, df))
  } else if (alternative == "greater") {
    pval <- stats::pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - stats::qt(ci.lvl, df), Inf)
  } else {
    pval <- 2 * stats::pt(-abs(tstat), df)
    alpha <- 1 - ci.lvl
    cint <- stats::qt(1 - alpha / 2, df)
    cint <- tstat + c(-cint, cint)
  }

  cint <- mu + cint * se

  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (paired || !is.null(yv)) "difference in means" else "mean"


  tt <- structure(
    class = "sj_ttest",
    list(
      estimate = estimate,
      statistic = tstat,
      df = df,
      p.value = pval,
      ci = cint,
      alternative = alternative,
      method = method
    )
  )

  attr(tt, "x.name") <- x.name
  attr(tt, "y.name") <- y.name
  attr(tt, "group.name") <- group.name

  tt
}
