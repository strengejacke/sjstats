#' @title Chi-Squared Test
#' @name chi_squared_test
#' @description This function performs a Mann-Whitney-Test (or Wilcoxon rank
#' sum test for _unpaired_ samples, see [`wilcox.test()`] and [`coin::wilcox_test()`]).
#'
#' The function reports p and Z-values as well as effect size r and group-rank-means.
#'
#' @param probabilities A numeric vector of probabilities for each cell in the
#' contingency table. The length of the vector must match the number of cells
#' in the table, i.e. the number of unique levels of the variable specified
#' in `select`. If `probabilities` is provided, a chi-squared test for given
#' probabilities is conducted. Furthermore, if `probabilities` is given, `by`
#' must be `NULL`. The probabilities must sum to 1.
#' @param ... Additional arguments passed down to [`chisq.test()`].
#' @inheritParams mann_whitney_test
#'
#' @return A data frame with test results.
#'
#' @examples
#' data(efc)
#' efc$weight <- abs(rnorm(nrow(efc), 1, 0.3))
#' # Chi-squared-test
#' chi_squared_test(efc, "c161sex", by = "e16sex")
#' # weighted Chi-squared-test
#' chi_squared_test(efc, "c161sex", by = "e16sex", weights = "weight")
#' # Chi-squared-test for given probabilities
#' chi_squared_test(efc, "c161sex", probabilities = c(0.3, 0.7))
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

.calculate_chisq <- function(data, select, by, weights, verbose = TRUE, ...) {
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
    tab <- as.table(round(stats::xtabs(x[[3]] ~ x[[1]] + x[[2]])))
    class(tab) <- "table"
  }

  # expected values, to identify whether Fisher's test is needed
  expected_values <- as.table(round(as.array(margin.table(tab, 1)) %*% t(as.array(margin.table(tab, 2))) / margin.table(tab))) # nolint

  # chi-squared test
  htest <- suppressWarnings(stats::chisq.test(tab, ...))
  test_statistic <- htest$statistic

  # need fisher?
  if (min(expected_values) < 5 || (min(expected_values) < 10 && htest$parameter == 1)) {
    htest <- stats::fisher.test(tab, simulate.p.value = TRUE, ...)
  }
  p_value <- htest$p.value

  # effect size
  if (nrow(tab) > 2 || ncol(tab) > 2) {
    effect_size <- stats::setNames(cramer(tab), "Cramer's V")
  } else {
    effect_size <- stats::setNames(phi(tab), "Phi")
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
  attr(out, "weighted") <- !is.null(weights)
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
    tab <- as.table(round(stats::xtabs(x[[2]] ~ x[[1]])))
    class(tab) <- "table"
  }

  # table dimensions
  n_rows <- nlevels(droplevels(as.factor(x$grp)))

  # sanity check
  if (length(probabilities) != n_rows) {
    insight::format_error("Length of probabilities must match number of cells in table (i.e. number of levels of input factor).") # nolint
  }
  if (!isTRUE(all.equal(sum(probabilities), 1))) {
    insight::format_error("Probabilities must sum to 1.")
  }

  # chi-squared test
  htest <- suppressWarnings(stats::chisq.test(tab, p = probabilities, rescale.p = TRUE, ...))
  test_statistic <- htest$statistic
  p_value <- htest$p.value

  effect_size <- effectsize::chisq_to_fei(
    test_statistic,
    n = sum(tab),
    nrow = n_rows,
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
    df = n_rows - 1,
    n_obs = sum(tab, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  class(out) <- c("sj_htest_chi", "data.frame")
  attr(out, "caption") <- "given Probabilities"
  attr(out, "weighted") <- !is.null(weights)
  out
}


# methods ---------------------------------------------------------------------

#' @export
print.sj_htest_chi <- function(x, ...) {
  weighted <- attributes(x)$weighted
  if (weighted) {
    weight_string <- " (weighted)"
  } else {
    weight_string <- ""
  }

  # get length of method name, to align output
  l <- max(nchar(c(x$statistic_name, x$effect_size_name, "p-value", "Observations")))
  # headline
  insight::print_color(sprintf(
    "\n# Chi-Squared Test for %s%s\n\n",
    attributes(x)$caption,
    weight_string
  ), "blue")

  stat_symbol <- .format_symbols(effect_size_name)

  # print test statistic
  cat(sprintf("  %*s: %.4f\n", l, x$statistic_name, x$statistic))
  cat(sprintf("  %*s: %.4f\n", l, stat_symbol, x$effect_size))
  cat(sprintf("  %*s: %g\n", l, "df", x$df))
  cat(sprintf("  %*s: %s\n", l, "p-value", insight::format_p(x$p, stars = TRUE, name = NULL)))
  cat(sprintf("  %*s: %g\n", l, "Observations", x$n_obs))
}
