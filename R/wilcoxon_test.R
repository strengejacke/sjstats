#' @title Wilcoxon rank sum test
#' @name wilcoxon_test
#' @description This function performs Wilcoxon rank sum tests for one sample
#' or for two _paired_ (dependent) samples. For _unpaired_ (independent)
#' samples, please use the `mann_whitney_test()` function.
#'
#' A Wilcoxon rank sum test is a non-parametric test for the null hypothesis
#' that two samples have identical continuous distributions. The implementation
#' in `wilcoxon_test()` is only used for _paired_, i.e. _dependent_ samples. For
#' independent (unpaired) samples, use `mann_whitney_test()`.
#'
#' `wilcoxon_test()` can be used for ordinal scales or when the continuous
#' variables are not normally distributed. For large samples, or approximately
#' normally distributed variables, the `t_test()` function can be used (with
#' `paired = TRUE`).
#'
#' @inheritParams mann_whitney_test
#' @inherit mann_whitney_test seealso
#'
#' @return A data frame with test results. The function returns p and Z-values
#' as well as effect size r and group-rank-means.
#'
#' @examplesIf requireNamespace("coin")
#' data(mtcars)
#' # one-sample test
#' wilcoxon_test(mtcars, "mpg")
#' # base R equivalent, we set exact = FALSE to avoid a warning
#' wilcox.test(mtcars$mpg ~ 1, exact = FALSE)
#'
#' # paired test
#' wilcoxon_test(mtcars, c("mpg", "hp"))
#' # base R equivalent, we set exact = FALSE to avoid a warning
#' wilcox.test(mtcars$mpg, mtcars$hp, paired = TRUE, exact = FALSE)
#'
#' # when `by` is specified, each group must be of same length
#' data(iris)
#' d <- iris[iris$Species != "setosa", ]
#' wilcoxon_test(d, "Sepal.Width", by = "Species")
#' @export
wilcoxon_test <- function(data,
                          select = NULL,
                          by = NULL,
                          weights = NULL,
                          mu = 0,
                          alternative = "two.sided",
                          ...) {
  insight::check_if_installed("datawizard")
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))

  # sanity checks
  .sanitize_htest_input(data, select, by, weights, test = "wilcoxon_test")

  # alternative only if weights are NULL
  if (!is.null(weights) && alternative != "two.sided") {
    insight::format_error("Argument `alternative` must be `two.sided` if `weights` are specified.")
  }

  # for paired two-sample, do groups all have same length?
  if (!is.null(by)) {
    group_len <- as.numeric(table(as.vector(data[[by]])))
    if (!all(group_len == group_len[1])) {
      insight::format_error("For paired two-sample Wilcoxon test, all groups specified in `by` must have the same length.") # nolint
    }
    # convert to wide format
    out <- split(data[select], as.character(data[[by]]))
    data <- stats::setNames(do.call(cbind, out), names(out))
    select <- colnames(data)
  }

  # value labels
  group_labels <- select

  x <- data[[select[1]]]
  if (length(select) > 1) {
    y <- data[[select[2]]]
  } else {
    y <- NULL
  }

  if (is.null(weights)) {
    .calculate_wilcox(x, y, alternative, mu, group_labels, ...)
  } else {
    .calculate_weighted_mwu(x, y, data[[weights]], group_labels)
  }
}


# Mann-Whitney-Test for two groups --------------------------------------------

.calculate_wilcox <- function(x, y, alternative, mu, group_labels, ...) {
  insight::check_if_installed("coin")
  # for paired Wilcoxon test, we have effect sizes
  if (!is.null(y)) {
    # prepare data
    wcdat <- data.frame(x, y)
    # perfom wilcox test
    wt <- coin::wilcoxsign_test(x ~ y, data = wcdat)
    # compute statistics
    u <- as.numeric(coin::statistic(wt, type = "linear"))
    z <- as.numeric(coin::statistic(wt, type = "standardized"))
    r <- abs(z / sqrt(nrow(wcdat)))
  } else {
    wt <- u <- z <- r <- NULL
  }

  # prepare data
  if (is.null(y)) {
    dv <- x
  } else {
    dv <- x - y
  }
  htest <- suppressWarnings(stats::wilcox.test(
    dv ~ 1,
    alternative = alternative,
    mu = mu,
    ...
  ))
  v <- htest$statistic
  p <- htest$p.value

  out <- data.frame(
    group1 = group_labels[1],
    v = v,
    p = as.numeric(p),
    mu = mu,
    alternative = alternative
  )
  # two groups?
  if (length(group_labels) > 1) {
    out$group2 <- group_labels[2]
  }
  # add effectsizes, when we have
  if (!is.null(wt)) {
    out$u <- u
    out$z <- z
    out$r <- r
  }
  attr(out, "group_labels") <- group_labels
  attr(out, "method") <- "wilcoxon"
  attr(out, "weighted") <- FALSE
  attr(out, "one_sample") <- length(group_labels) == 1
  class(out) <- c("sj_htest_wilcox", "data.frame")

  out
}


# Weighted Mann-Whitney-Test for two groups ----------------------------------

.calculate_weighted_wilcox <- function(x, y, weights, group_labels) {
  # check if pkg survey is available
  insight::check_if_installed("survey")

  # prepare data
  if (is.null(y)) {
    dv <- x
  } else {
    dv <- x - y
  }

  dat <- stats::na.omit(data.frame(dv, weights))
  colnames(dat) <- c("y", "w")

  design <- survey::svydesign(ids = ~0, data = dat, weights = ~w)
  result <- survey::svyranktest(formula = y ~ 1, design, test = "wilcoxon")

  # statistics and effect sizes
  z <- result$statistic
  r <- abs(z / sqrt(nrow(dat)))

  out <- data_frame(
    group1 = group_labels[1],
    estimate = result$estimate,
    z = z,
    r = r,
    p = as.numeric(result$p.value),
    mu = 0,
    alternative = "two.sided"
  )
  # two groups?
  if (length(group_labels) > 1) {
    out$group2 <- group_labels[2]
  }

  attr(out, "group_labels") <- group_labels
  attr(out, "weighted") <- TRUE
  attr(out, "one_sample") <- length(group_labels) == 1
  attr(out, "method") <- "wilcoxon"
  class(out) <- c("sj_htest_wilcox", "data.frame")

  out
}


# methods ---------------------------------------------------------------------

#' @export
print.sj_htest_wilcox <- function(x, ...) {
  # fetch attributes
  group_labels <- attributes(x)$group_labels
  weighted <- attributes(x)$weighted
  one_sample <- attributes(x)$one_sample

  if (weighted) {
    weight_string <- " (weighted)"
  } else {
    weight_string <- ""
  }

  if (one_sample) {
    onesample_string <- "One Sample"
  } else {
    onesample_string <- "Paired"
  }

  # same width
  group_labels <- format(group_labels)

  # header
  insight::print_color(sprintf(
    "# %s Wilcoxon signed rank test%s\n\n",
    onesample_string,
    weight_string
  ), "blue")

  # alternative hypothesis
  if (!is.null(x$alternative) && !is.null(x$mu)) {
    alt_string <- switch(x$alternative,
      two.sided = "not equal to",
      less = "less than",
      greater = "greater than"
    )
    alt_string <- paste("true location shift is", alt_string, x$mu)
    insight::print_color(sprintf("  Alternative hypothesis: %s\n", alt_string), "cyan")
  }

  if (!is.null(x[["v"]])) {
    v_stat <- sprintf("V = %i, ", round(x$v))
  } else {
    v_stat <- ""
  }

  if (!is.null(x[["r"]])) {
    cat(sprintf("\n  %sr = %.2f, Z = %.2f, %s\n\n", v_stat, x$r, x$z, insight::format_p(x$p)))
  } else {
    cat(sprintf("\n  %s%s\n\n", v_stat, insight::format_p(x$p)))
  }
}
