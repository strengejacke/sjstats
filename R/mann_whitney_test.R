#' @title Mann-Whitney test
#' @name mann_whitney_test
#' @description This function performs a Mann-Whitney test (or Wilcoxon rank
#' sum test for _unpaired_ samples.
#'
#' A Mann-Whitney test is a non-parametric test for the null hypothesis that two
#' independent samples have identical continuous distributions. It can be used
#' when the two continuous variables are not normally distributed.
#'
#' @param data A data frame.
#' @param select One or more name of the continuous variable (as character
#' vector) to be used as samples for the test. If `select` only specified one
#' variable, a one-sample test is carried out (only applicable for `t_test()`).
#' Else, `by` must be provided to indicate the groups of comparison.
#' @param by Name of the variable indicating the groups. Required if `select`
#' specifies only one variable that contains all samples to be compared in the
#' test. If `by` is not a factor, it will be coerced to a factor. For
#' `chi_squared_test()`, if `probabilities` is provided, `by` must be `NULL`.
#' @param weights Name of an (optional) weighting variable to be used for the test.
#' @param ... Additional arguments passed to `wilcox.test()` (for unweighted
#' tests, i.e. when `weights = NULL`).
#'
#' @return A data frame with test results. The function returns p and Z-values
#' as well as effect size r and group-rank-means.
#'
#' @details This function is based on [`wilcox.test()`] and [`coin::wilcox_test()`]
#' (the latter to extract effect sizes). The weighted version of the test is
#' based on [`survey::svyranktest()`].
#'
#' Interpretation of the effect size **r**, as a rule-of-thumb:
#'
#' - small effect >= 0.1
#' - medium effect >= 0.3
#' - large effect >= 0.5
#'
#' **r** is calcuated as:
#'
#' ```
#' r = |Z| / sqrt(n1 + n2)
#' ```
#'
#' @examplesIf requireNamespace("coin", quietly = TRUE) && requireNamespace("survey", quietly = TRUE)
#' data(efc)
#' # Mann-Whitney-U tests for elder's age by elder's sex.
#' mann_whitney_test(efc, "e17age", by = "e16sex")
#' # base R equivalent
#' wilcox.test(e17age ~ e16sex, data = efc)
#'
#' # when data is in wide-format, specify all relevant continuous
#' # variables in `select` and omit `by`
#' set.seed(123)
#' wide_data <- data.frame(scale1 = runif(20), scale2 = runif(20))
#' mann_whitney_test(wide_data, select = c("scale1", "scale2"))
#' # base R equivalent
#' wilcox.test(wide_data$scale1, wide_data$scale2)
#' # same as if we had data in long format, with grouping variable
#' long_data <- data.frame(
#'   scales = c(wide_data$scale1, wide_data$scale2),
#'   groups = as.factor(rep(c("A", "B"), each = 20))
#' )
#' mann_whitney_test(long_data, select = "scales", by = "groups")
#' # base R equivalent
#' wilcox.test(scales ~ groups, long_data)
#' @export
mann_whitney_test <- function(data,
                              select = NULL,
                              by = NULL,
                              weights = NULL,
                              ...) {
  insight::check_if_installed("datawizard")

  # sanity checks
  .sanitize_htest_input(data, select, by, weights)

  # does select indicate more than one variable?
  if (length(select) > 1) {
    # sanity check - may only specify two variable names
    if (length(select) > 2) {
      insight::format_error("You may only specify two variables for Mann-Whitney test.")
    }
    if (!is.null(by)) {
      insight::format_error("If `select` specifies more than one variable, `by` must be `NULL`.")
    }
    # we convert the data into long format, and create a grouping variable
    data <- datawizard::data_to_long(data[select], names_to = "group", values_to = "scale")
    by <- select[2]
    select <- select[1]
    # after converting to long, we have the "grouping" variable first in the data
    colnames(data) <- c(by, select)
  }

  # get data
  dv <- data[[select]]
  grp <- data[[by]]

  # coerce to factor
  grp <- datawizard::to_factor(grp)

  # only two groups allowed
  if (insight::n_unique(grp) > 2) {
    insight::format_error("Only two groups are allowed for Mann-Whitney test. Please use `kruskal_wallis_test()` for more than two groups.") # nolint
  }

  # value labels
  group_labels <- names(attr(data[[by]], "labels", exact = TRUE))
  if (is.null(group_labels)) {
    group_labels <- levels(droplevels(grp))
  }

  if (is.null(weights)) {
    .calculate_mwu(dv, grp, distribution, group_labels, ...)
  } else {
    .calculate_weighted_mwu(dv, grp, data[[weights]], group_labels)
  }
}


# Mann-Whitney-Test for two groups --------------------------------------------

.calculate_mwu <- function(dv, grp, distribution, group_labels, ...) {
  insight::check_if_installed("coin")
  # prepare data
  wcdat <- data.frame(dv, grp)
  # perfom wilcox test
  wt <- coin::wilcox_test(dv ~ grp, data = wcdat, distribution = distribution)

  # for rank mean
  group_levels <- levels(grp)

  # compute statistics
  u <- as.numeric(coin::statistic(wt, type = "linear"))
  z <- as.numeric(coin::statistic(wt, type = "standardized"))
  r <- abs(z / sqrt(length(dv)))
  htest <- suppressWarnings(stats::wilcox.test(dv ~ grp, data = wcdat, ...))
  w <- htest$statistic
  p <- htest$p.value

  # group means
  dat_gr1 <- stats::na.omit(dv[grp == group_levels[1]])
  dat_gr2 <- stats::na.omit(dv[grp == group_levels[2]])

  rank_mean_1 <- mean(rank(dat_gr1))
  rank_mean_2 <- mean(rank(dat_gr2))

  # compute n for each group
  n_grp1 <- length(dat_gr1)
  n_grp2 <- length(dat_gr2)

  out <- data.frame(
    group1 = group_levels[1],
    group2 = group_levels[2],
    estimate = rank_mean_1 - rank_mean_2,
    u = u,
    w = w,
    z = z,
    r = r,
    p = as.numeric(p)
  )
  attr(out, "rank_means") <- stats::setNames(
    c(rank_mean_1, rank_mean_2),
    c("Mean Group 1", "Mean Group 2")
  )
  attr(out, "n_groups") <- stats::setNames(
    c(n_grp1, n_grp2),
    c("N Group 1", "N Group 2")
  )
  attr(out, "group_labels") <- group_labels
  attr(out, "method") <- "wilcoxon"
  attr(out, "weighted") <- FALSE
  class(out) <- c("sj_htest_mwu", "data.frame")

  out
}


# Weighted Mann-Whitney-Test for two groups ----------------------------------

.calculate_weighted_mwu <- function(dv, grp, weights, group_labels) {
  # check if pkg survey is available
  insight::check_if_installed("survey")

  dat <- stats::na.omit(data.frame(dv, grp, weights))
  colnames(dat) <- c("x", "g", "w")

  design <- survey::svydesign(ids = ~0, data = dat, weights = ~w)
  result <- survey::svyranktest(formula = x ~ g, design, test = "wilcoxon")

  # for rank mean
  group_levels <- levels(droplevels(grp))
  # subgroups
  dat_gr1 <- dat[dat$g == group_levels[1], ]
  dat_gr2 <- dat[dat$g == group_levels[2], ]
  dat_gr1$rank_x <- rank(dat_gr1$x)
  dat_gr2$rank_x <- rank(dat_gr2$x)

  # rank means
  design_mean1 <- survey::svydesign(
    ids = ~0,
    data = dat_gr1,
    weights = ~w
  )
  rank_mean_1 <- survey::svymean(~rank_x, design_mean1)

  design_mean2 <- survey::svydesign(
    ids = ~0,
    data = dat_gr2,
    weights = ~w
  )
  rank_mean_2 <- survey::svymean(~rank_x, design_mean2)

  # group Ns
  n_grp1 <- round(sum(dat_gr1$w))
  n_grp2 <- round(sum(dat_gr2$w))

  # statistics and effect sizes
  z <- result$statistic
  r <- abs(z / sqrt(sum(n_grp1, n_grp2)))

  out <- data.frame(
    group1 = group_levels[1],
    group2 = group_levels[2],
    estimate = result$estimate,
    z = z,
    r = r,
    p = as.numeric(result$p.value)
  )

  attr(out, "rank_means") <- stats::setNames(
    c(rank_mean_1, rank_mean_2),
    c("Mean Group 1", "Mean Group 2")
  )
  attr(out, "n_groups") <- stats::setNames(
    c(n_grp1, n_grp2),
    c("N Group 1", "N Group 2")
  )
  attr(out, "group_labels") <- group_labels
  attr(out, "weighted") <- TRUE
  class(out) <- c("sj_htest_mwu", "data.frame")

  out
}


# helper ----------------------------------------------------------------------

.sanitize_htest_input <- function(data, select, by, weights) {
  # check if arguments are NULL
  if (is.null(select)) {
    insight::format_error("Argument `select` is missing.")
  }

  # check if arguments have correct length or are of correct type
  if (!is.character(select)) {
    insight::format_error("Argument `select` must be a character string with the name(s) of the variable(s).")
  }
  if (!is.null(by) && (length(by) != 1 || !is.character(by))) {
    insight::format_error("Argument `by` must be a character string with the name of a single variable.")
  }
  if (!is.null(weights) && length(weights) != 1) {
    insight::format_error("Argument `weights` must be a character string with the name of a single variable.")
  }

  # check if "select" is in data
  if (!all(select %in% colnames(data))) {
    not_found <- setdiff(select, colnames(data))[1]
    insight::format_error(
      sprintf("Variable '%s' not found in data frame.", not_found),
      .misspelled_string(colnames(data), not_found, "Maybe misspelled?")
    )
  }
  # check if "by" is in data
  if (!is.null(by) && !by %in% colnames(data)) {
    insight::format_error(
      sprintf("Variable '%s' not found in data frame.", by),
      .misspelled_string(colnames(data), by, "Maybe misspelled?")
    )
  }
  # check if "weights" is in data
  if (!is.null(weights) && !weights %in% colnames(data)) {
    insight::format_error(
      sprintf("Weighting variable '%s' not found in data frame.", weights),
      .misspelled_string(colnames(data), weights, "Maybe misspelled?")
    )
  }
}


# methods ---------------------------------------------------------------------

#' @export
print.sj_htest_mwu <- function(x, ...) {
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

  cat(sprintf("\n  r = %.3f, Z = %.3f, %s\n\n", x$r, x$z, insight::format_p(x$p)))
}
