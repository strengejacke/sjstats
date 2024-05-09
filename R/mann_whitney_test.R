#' @title Mann-Whitney-Test
#' @name mann_whitney_test
#' @description This function performs a Mann-Whitney-Test (or Wilcoxon rank
#' sum test for _unpaired_ samples.
#'
#' A Mann-Whitney-Test is a non-parametric test for the null hypothesis that two
#' independent samples have identical continuous distributions. It can be used
#' when the two continuous variables are not normally distributed.
#'
#' @param data A data frame.
#' @param select Name of the dependent variable (as string) to be used for the
#' test. `select` can also be a character vector, specifing the names of
#' multiple continuous variables. In this case, `by` is ignored and variables
#' specified in `select` are used to compute the test. This can be useful if
#' the data is in wide-format and no grouping variable is available.
#' @param by Name of the grouping variable to be used for the test. If `by` is
#' not a factor, it will be coerced to a factor. For `chi_squared_test()`, if
#' `probabilities` is provided, `by` must be `NULL`.
#' @param weights Name of an (optional) weighting variable to be used for the test.
#' @param distribution Indicates how the null distribution of the test statistic
#' should be computed. May be one of `"exact"`, `"approximate"` or `"asymptotic"`
#' (default). See [`coin::wilcox_test()`] for details.
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
#' @examples
#' data(efc)
#' # Mann-Whitney-U-Tests for elder's age by elder's sex.
#' mann_whitney_test(efc, "e17age", by = "e16sex")
#'
#' # when data is in wide-format, specify all relevant continuous
#' # variables in `select` and omit `by`
#' set.seed(123)
#' wide_data <- data.frame(scale1 = runif(20), scale2 = runif(20))
#' mann_whitney_test(wide_data, select = c("scale1", "scale2"))
#'
#' # same as if we had data in long format, with grouping variable
#' long_data <- data.frame(
#'   scales = c(wide_data$scale1, wide_data$scale2),
#'   groups = rep(c("A", "B"), each = 20)
#' )
#' mann_whitney_test(long_data, select = "scales", by = "groups")
#' @export
mann_whitney_test <- function(data,
                              select = NULL,
                              by = NULL,
                              weights = NULL,
                              distribution = "asymptotic") {
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
    .calculate_mwu(dv, grp, distribution, group_labels)
  } else {
    .calculate_weighted_mwu(dv, grp, data[[weights]], group_labels)
  }
}


# Mann-Whitney-Test for two groups --------------------------------------------

.calculate_mwu <- function(dv, grp, distribution, group_labels) {
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
  p <- coin::pvalue(wt)
  r <- abs(z / sqrt(length(dv)))
  w <- suppressWarnings(stats::wilcox.test(dv ~ grp, data = wcdat)$statistic)

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

.misspelled_string <- function(source, searchterm, default_message = NULL) {
  if (is.null(searchterm) || length(searchterm) < 1) {
    return(default_message)
  }
  # used for many matches
  more_found <- ""
  # init default
  msg <- ""
  # remove matching strings
  same <- intersect(source, searchterm)
  searchterm <- setdiff(searchterm, same)
  source <- setdiff(source, same)
  # guess the misspelled string
  possible_strings <- unlist(lapply(searchterm, function(s) {
    source[.fuzzy_grep(source, s)] # nolint
  }), use.names = FALSE)
  if (length(possible_strings)) {
    msg <- "Did you mean "
    if (length(possible_strings) > 1) {
      # make sure we don't print dozens of alternatives for larger data frames
      if (length(possible_strings) > 5) {
        more_found <- sprintf(
          " We even found %i more possible matches, not shown here.",
          length(possible_strings) - 5
        )
        possible_strings <- possible_strings[1:5]
      }
      msg <- paste0(msg, "one of ", toString(paste0("\"", possible_strings, "\"")))
    } else {
      msg <- paste0(msg, "\"", possible_strings, "\"")
    }
    msg <- paste0(msg, "?", more_found)
  } else {
    msg <- default_message
  }
  # no double white space
  insight::trim_ws(msg)
}


.fuzzy_grep <- function(x, pattern, precision = NULL) {
  if (is.null(precision)) {
    precision <- round(nchar(pattern) / 3)
  }
  if (precision > nchar(pattern)) {
    return(NULL)
  }
  p <- sprintf("(%s){~%i}", pattern, precision)
  grep(pattern = p, x = x, ignore.case = FALSE)
}


.sanitize_htest_input <- function(data, select, by, weights) {
  # check if arguments are NULL
  if (is.null(select)) {
    insight::format_error("Argument `select` is missing.")
  }
  # `by` is only allowed to be NULL if `select` specifies more than one variable
  if (is.null(by) && length(select) == 1) {
    insight::format_error("Arguments `by` is missing.")
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
