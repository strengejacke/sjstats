#' @title Wilcoxon rank sum test
#' @name wilcoxon_test
#' @description This function performs Wilcoxon rank sum tests for one sample
#' or for two _paired_ (dependent) samples. For _unpaired_ (independent)
#' samples, please use the `mann_whitney_test()` function.
#'
#' @inheritParams mann_whitney_test
#'
#' @return A data frame with test results. The function returns p and Z-values
#' as well as effect size r and group-rank-means.
#'
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
    group_len <- as.numeric(table(data[[by]]))
    if (!all(group_len == group_len[1])) {
      insight::format_error("For paired two-sample Wilcoxon test, all groups specified in `by` must have the same length.") # nolint
    }
    # convert to wide format
    data <- datawizard::data_to_wide(data, values_from = select, names_from = by)
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
    .calculate_weighted_mwu(dv, grp, data[[weights]], group_labels)
  }
}


# Mann-Whitney-Test for two groups --------------------------------------------

.calculate_wilcox <- function(x, y, alternative, mu, group_labels, ...) {
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
  attr(out, "one_sample") <- length(group_labels) > 1
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
    alternative = "two.sided"
  )
  # two groups?
  if (length(group_labels) > 1) {
    out$group2 <- select[2]
  }

  attr(out, "group_labels") <- group_labels
  attr(out, "weighted") <- TRUE
  attr(out, "one_sample") <- length(group_labels) > 1
  attr(out, "method") <- "wilcoxon"
  class(out) <- c("sj_htest_wilcox", "data.frame")

  out
}


# methods ---------------------------------------------------------------------

#' @export
print.sj_htest_wilcox <- function(x, ...) {
  # fetch attributes
  group_labels <- attributes(x)$group_labels
  rank_means <- attributes(x)$rank_means
  n_groups <- attributes(x)$n_groups
  weighted <- attributes(x)$weighted

  if (weighted) {
    weight_string <- " (weighted)"
  } else {
    weight_string <- ""
  }

  # same width
  group_labels <- format(group_labels)

  # header
  insight::print_color(sprintf("# Mann-Whitney test%s\n\n", weight_string), "blue")

  # group-1-info
  insight::print_color(
    sprintf(
      "  Group 1: %s (n = %i, rank mean = %s)\n",
      group_labels[1], n_groups[1], insight::format_value(rank_means[1], protect_integers = TRUE)
    ), "cyan"
  )

  # group-2-info
  insight::print_color(
    sprintf(
      "  Group 2: %s (n = %i, rank mean = %s)\n",
      group_labels[2], n_groups[2], insight::format_value(rank_means[2], protect_integers = TRUE)
    ), "cyan"
  )

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

  cat(sprintf("\n  r = %.3f, Z = %.3f, %s\n\n", x$r, x$z, insight::format_p(x$p)))
}
