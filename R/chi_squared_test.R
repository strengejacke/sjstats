#' @title Chi-Squared test
#' @name chi_squared_test
#' @description This function performs a \eqn{chi}^2 test for contingency
#' tables or tests for given probabilities. The returned effects sizes are
#' Cramer's V for tables with more than two rows and columns, Phi (\eqn{\phi})
#' for 2x2 tables, and \ifelse{latex}{\eqn{Fei}}{פ (Fei)} for tests against
#' given probabilities (see _Ben-Shachar et al. 2023_).
#'
#' @param probabilities A numeric vector of probabilities for each cell in the
#' contingency table. The length of the vector must match the number of cells
#' in the table, i.e. the number of unique levels of the variable specified
#' in `select`. If `probabilities` is provided, a chi-squared test for given
#' probabilities is conducted. Furthermore, if `probabilities` is given, `by`
#' must be `NULL`. The probabilities must sum to 1.
#' @param paired Logical, if `TRUE`, a McNemar test is conducted for 2x2 tables.
#' Note that `paired` only works for 2x2 tables.
#' @param ... Additional arguments passed down to [`chisq.test()`].
#' @inheritParams mann_whitney_test
#'
#' @return A data frame with test results. The returned effects sizes are
#' Cramer's V for tables with more than two rows and columns, Phi (\eqn{\phi})
#' for 2x2 tables, and \ifelse{latex}{\eqn{Fei}}{פ (Fei)} for tests against
#' given probabilities.
#'
#' @details The function is a wrapper around [`chisq.test()`] and
#' [`fisher.test()`] (for small expected values) for contingency tables, and
#' `chisq.test()` for given probabilities. When `probabilities` are provided,
#' these are rescaled to sum to 1 (i.e. `rescale.p = TRUE`). When `fisher.test()`
#' is called, simulated p-values are returned (i.e. `simulate.p.value = TRUE`,
#' see `?fisher.test`). If `paired = TRUE` and a 2x2 table is provided,
#' a McNemar test (see [`mcnemar.test()`]) is conducted.
#'
#' The weighted version of the chi-squared test is based on the a weighted
#' table, using [`xtabs()`] as input for `chisq.test()`.
#'
#' @references Ben-Shachar, M.S., Patil, I., Thériault, R., Wiernik, B.M.,
#' Lüdecke, D. (2023). Phi, Fei, Fo, Fum: Effect Sizes for Categorical Data
#' That Use the Chi‑Squared Statistic. Mathematics, 11, 1982.
#' \doi{10.3390/math11091982}
#'
#' @examplesIf requireNamespace("effectsize")
#' data(efc)
#' efc$weight <- abs(rnorm(nrow(efc), 1, 0.3))
#'
#' # Chi-squared test
#' chi_squared_test(efc, "c161sex", by = "e16sex")
#'
#' # weighted Chi-squared test
#' chi_squared_test(efc, "c161sex", by = "e16sex", weights = "weight")
#'
#' # Chi-squared test for given probabilities
#' chi_squared_test(efc, "c161sex", probabilities = c(0.3, 0.7))
#' @export
chi_squared_test <- function(data,
                             select = NULL,
                             by = NULL,
                             probabilities = NULL,
                             weights = NULL,
                             paired = FALSE,
                             ...) {
  if (is.null(probabilities)) {
    .calculate_chisq(data, select, by, weights, paired, ...)
  } else {
    # sanity check - `paired = TRUE` is not available for given probabilities
    if (paired) {
      insight::format_error("When `probabilities` are provided, `paired = TRUE` is not available.") # nolint
    }
    .calculate_chisq_gof(data, select, probabilities, weights, ...)
  }
}


# Mann-Whitney-Test for two groups --------------------------------------------

.calculate_chisq <- function(data, select, by, weights, paired = FALSE, ...) {
  insight::check_if_installed("datawizard")
  # sanity checks
  .sanitize_htest_input(data, select, by, weights)

  # get data
  grp1 <- data[[select]]
  grp2 <- data[[by]]

  # if paired = TRUE, we only allow a 2x2 table
  if (paired && (length(stats::na.omit(unique(grp1))) != 2 || length(stats::na.omit(unique(grp2))) != 2)) {
    insight::format_error("When `paired = TRUE`, only 2x2 tables are allowed (i.e. both variables must have exactly two levels).") # nolint
  }

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

  # paired? mc-nemar test
  if (paired) {
    htest <- suppressWarnings(stats::mcnemar.test(tab, ...))
    test_statistic <- htest$statistic
  } else {
    # chi-squared test
    htest <- suppressWarnings(stats::chisq.test(tab, ...))
    test_statistic <- htest$statistic
    # need fisher?
    if (min(expected_values) < 5 || (min(expected_values) < 10 && htest$parameter == 1)) {
      htest <- stats::fisher.test(tab, simulate.p.value = TRUE, ...)
    }
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
    data = paste(select, "by", by),
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
  attr(out, "fisher") <- isTRUE(startsWith(htest$method, "Fisher"))
  attr(out, "mcnemar") <- isTRUE(paired)
  attr(out, "caption") <- "contingency tables"
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
    data = paste(
      select,
      "against probabilities",
      datawizard::text_concatenate(sprintf("%i%%", round(100 * probabilities)))
    ),
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
  attr(out, "caption") <- "given probabilities"
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

  fisher <- attributes(x)$fisher
  mcnemar <- attributes(x)$mcnemar

  # headline
  insight::print_color(sprintf(
    "\n# Chi-squared test for %s%s\n",
    attributes(x)$caption,
    weight_string
  ), "blue")

  # Fisher's exact test?
  if (isTRUE(fisher)) {
    insight::print_color("  (using Fisher's exact test due to small expected values)\n", "blue") # nolint
  } else if (isTRUE(mcnemar)) {
    insight::print_color("  (using McNemar's test for paired data)\n", "blue") # nolint
  }

  cat("\n")

  # data info
  insight::print_color(
    sprintf("  Data: %s (n = %i)\n", x$data, round(x$n_obs)),
    "cyan"
  )

  # prepare and align strings
  eff_symbol <- .format_symbols(x$effect_size_name)
  stat_symbol <- .format_symbols(x$statistic_name)

  cat(sprintf(
    "\n  %s = %.4f, %s = %.4f, df = %i, %s\n\n",
    stat_symbol, x$statistic, eff_symbol, x$effect_size, round(x$df), insight::format_p(x$p)
  ))
}
