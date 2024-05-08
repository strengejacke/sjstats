#' @title Chi-Squared Test
#' @name chi_squared_test
#' @description This function performs a Mann-Whitney-Test (or Wilcoxon rank
#' sum test for _unpaired_ samples, see [`wilcox.test()`] and [`coin::wilcox_test()`]).
#'
#' The function reports p and Z-values as well as effect size r and group-rank-means.
#'
#' @param data A data frame.
#' @param select The dependent variable (numeric) to be used for the test.
#' @param by The grouping variable (factor) to be used for the test. If `by` is
#' not a factor, it will be coerced to a factor.
#' @param weights An optional weighting variable (numeric) to be used for the test.
#' @param ... Additional arguments passed down to [`chisq.test()`].
#'
#' @return A data frame.
#'
#' @details This function calls [`coin::wilcox_test()`] to extract effect sizes.
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
#' chi_squared_test(efc, "e17age", "e16sex")
#' @export
chi_squared_test <- function(data,
                             select = NULL,
                             by = NULL,
                             probabilities = NULL,
                             weights = NULL,
                             ...) {
  if (is.null(probabilities)) {
    .calculate_chisq(data, select, by, weights, ...)
  } else {
    .calculate_chisq_gof(data, select, probabilities, weights, ...)
  }
}


# Mann-Whitney-Test for two groups --------------------------------------------

.calculate_chisq <- function(data, select, by, weights, ...) {
  insight::check_if_installed("datawizard")
  # sanity checks
  .sanitize_htest_input(data, select, by, weights)

  # get data
  grp1 <- data[[select]]
  grp2 <- data[[by]]

  # create data frame for table
  x <- data.frame(
    grp1 = datawizard::to_factor(grp1),
    grp2 = datawizard::to_factor(grp2)
  )
  # add weights
  if (!is.null(weights)) {
    x$weights <- data[[weights]]
  }
  # remove missings
  x <- stats::na.omit(x)

  # contingency table
  if (is.null(weights)) {
    tab <- table(x)
  } else {
    tab <- as.table(round(stats::xtabs(data[[3]] ~ data[[1]] + data[[2]])))
    class(tab) <- "table"
  }

  # expected values, to identify whether Fisher's test is needed
  expected_values <- as.table(round(as.array(margin.table(tab, 1)) %*% t(as.array(margin.table(tab, 2))) / margin.table(tab))) # nolint

  # chi-squared test
  htest <- suppressWarnings(stats::chisq.test(tab, ...))
  test_statistic <- htest$statistic

  # need fisher?
  if (min(expected_values$expected) < 5 || (min(expected_values$expected) < 10 && htest$parameter == 1)) {
    htest <- stats::fisher.test(tab, simulate.p.value = TRUE, ...)
  }
  p_value <- htest$p.value

  # effect size
  if (nrow(tab) > 2 || ncol(tab) > 2) {
    effect_size <- stats::setNamed(cramer(tab), "Cramer's V")
  } else {
    effect_size <- stats::setNamed(phi(tab), "Phi")
  }

  # return result
  out <- data.frame(
    statistic_name = "Chi-squared",
    statistic = test_statistic,
    effect_size_name = names(effect_size),
    effect_size = as.numeric(effect_size),
    p = p_value,
    df = (nrow(tab) - 1) * (ncol(tab) - 1),
    n_obs = sum(tab, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  class(out) <- c("sj_htest_chi", "data.frame")
  attr(out, "caption") <- "Contingency Tables"
  out
}


.calculate_chisq_gof <- function(data, select, probabilities, weights, ...) {
  insight::check_if_installed("effectsize")

  # get data
  x <- data.frame(grp = data[[select]])
  # add weights
  if (!is.null(weights)) {
    x$weights <- data[[weights]]
  }
  # remove missings
  x <- stats::na.omit(x)

  # contingency table
  if (is.null(weights)) {
    tab <- table(x)
  } else {
    tab <- as.table(round(stats::xtabs(data[[2]] ~ data[[1]])))
    class(tab) <- "table"
  }

  # sanity check
  if (length(probabilities) != as.numeric(dim(tab))) {
    insight::format_error("Length of probabilities must match number of cells in table (i.e. number of levels of input factor).") # nolint
  }

  # chi-squared test
  htest <- suppressWarnings(stats::chisq.test(tab, p = probabilities, rescale.p = TRUE, ...))
  test_statistic <- htest$statistic
  p_value <- htest$p.value

  effect_size <- effectsize::chisq_to_fei(
    test_statistic,
    n = sum(tab),
    nrow = as.numeric(dim(tab)),
    ncol = 1,
    p = probabilities,
    alternative = "two.sided"
  )$Fei

  # return result
  out <- data.frame(
    statistic_name = "Chi-squared",
    statistic = test_statistic,
    effect_size_name = "Fei",
    effect_size = as.numeric(effect_size),
    p = p_value,
    df = (nrow(tab) - 1) * (ncol(tab) - 1),
    n_obs = sum(tab, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  class(out) <- c("sj_htest_chi", "data.frame")
  attr(out, "caption") <- "given Probabilities"
  out
}


# methods ---------------------------------------------------------------------

#' @export
print.sj_htest_chi <- function(x, ...) {
  # get length of method name, to align output
  l <- max(nchar(c(x$statistic_name, x$effect_size_name, "p-value", "Observations")))
  # headline
  insight::print_color(sprintf("\n# Chi-squared Test for%s\n\n", attributes(x)$caption), "blue")

  # print test statistic
  cat(sprintf("  %*s: %.4f\n", l, x$statistic_name, x$statistic))
  cat(sprintf("  %*s: %.4f\n", l, x$effect_size_name, x$effect_size))
  cat(sprintf("  %*s: %g\n", l, "df", x$df))
  cat(sprintf("  %*s: %s\n", l, "p-value", insight::format_p(x$p.value, stars = TRUE, name = NULL)))
  cat(sprintf("  %*s: %g\n", l, "Observations", x$n_obs))
}
