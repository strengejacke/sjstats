#' @title Kruskal-Wallis test
#' @name kruskal_wallis_test
#' @description This function performs a Kruskal-Wallis rank sum test, to test
#' the null hypothesis that the population median of all of the groups are
#' equal. The alternative is that they differ in at least one. If `paired = TRUE`,
#' a paired Friedman test is conducted.
#'
#' @inheritParams mann_whitney_test
#' @param paired Logical, if `TRUE`, a paired Friedman test is conducted (see
#' [`friedman.test()`]).
#'
#' @return A data frame with test results.
#'
#' @details The function simply is a wrapper around [`kruskal.test()`]. The
#' weighted version of the Kruskal-Wallis test is based on the `survey` package,
#' using [`survey::svyranktest()`]. When `paired = TRUE`, a paired Friedman test
#' is conducted (see [`friedman.test()`]).
#'
#' @examples
#' data(efc)
#' # Kruskal-Wallis test for elder's age by education
#' kruskal_wallis_test(efc, "e17age", by = "c172code")
#' @export
kruskal_wallis_test <- function(data,
                                select = NULL,
                                by = NULL,
                                weights = NULL,
                                paired = FALSE) {
  insight::check_if_installed("datawizard")

  # sanity checks
  .sanitize_htest_input(data, select, by, weights)

  # get data
  dv <- data[[select]]
  grp <- data[[by]]

  # coerce to factor
  grp <- datawizard::to_factor(grp)

  # only two groups allowed
  if (insight::n_unique(grp) < 2) {
    insight::format_error("At least two groups are required, i.e. data must have at least two unique levels in `by` for `kruskal_wallis_test()`.") # nolint
  }
  if (is.null(weights)) {
    .calculate_kw(dv, grp, paired)
  } else {
    .calculate_weighted_kw(dv, grp, data[[weights]], paired = TRUE)
  }
}


# Kruskal-Wallis-Test --------------------------------------------

.calculate_kw <- function(dv, grp, paired = FALSE) {
  # prepare data
  wcdat <- data.frame(dv, grp)
  if (paired) {
    # perfom friedman test for paired data
    wt <- stats::friedman.test(table(wcdat))
  } else {
    # perfom kruskal wallis test
    wt <- stats::kruskal.test(dv ~ grp, data = wcdat)
  }
  # number of groups
  n_groups <- vapply(
    stats::na.omit(unique(grp)),
    function(g) sum(grp == g, na.rm = TRUE),
    numeric(1)
  )

  out <- data.frame(
    data = wt$data.name,
    Chi2 = wt$statistic,
    df = wt$parameter,
    p = as.numeric(wt$p.value),
    stringsAsFactors = FALSE
  )

  attr(out, "n_groups") <- n_groups
  attr(out, "method") <- ifelse(paired, "friedman", "kruskal")
  attr(out, "weighted") <- FALSE
  class(out) <- c("sj_htest_kw", "data.frame")

  out
}


# Weighted Mann-Whitney-Test for two groups ----------------------------------

.calculate_weighted_kw <- function(dv, grp, weights, paired = FALSE) {
  # check if pkg survey is available
  insight::check_if_installed("survey")

  dat <- stats::na.omit(data.frame(dv, grp, weights))
  colnames(dat) <- c("x", "g", "w")

  # number of groups
  n_groups <- vapply(stats::na.omit(unique(grp)), function(g) {
    sum(dat$w[dat$grp == g], na.rm = TRUE)
  }, numeric(1))

  if (paired) {
    tab <- as.table(round(stats::xtabs(x[[3]] ~ x[[1]] + x[[2]])))
    class(tab) <- "table"
    # perfom friedman test for paired data
    result <- stats::friedman.test(tab)
  } else {
    design <- survey::svydesign(ids = ~0, data = dat, weights = ~w)
    result <- survey::svyranktest(formula = x ~ g, design, test = "KruskalWallis")
  }

  out <- data.frame(
    data = paste(dv, "by", grp),
    Chi2 = result$statistic,
    df = result$parameter,
    p = as.numeric(result$p.value),
    stringsAsFactors = FALSE
  )

  attr(out, "n_groups") <- n_groups
  attr(out, "method") <- ifelse(paired, "friedman", "kruskal")
  attr(out, "weighted") <- TRUE
  class(out) <- c("sj_htest_kw", "data.frame")

  out
}


# methods ---------------------------------------------------------------------

#' @export
print.sj_htest_kw <- function(x, ...) {
  insight::check_if_installed("datawizard")
  # fetch attributes
  n_groups <- attributes(x)$n_groups
  weighted <- attributes(x)$weighted
  method <- attributes(x)$method

  if (weighted) {
    weight_string <- " (weighted)"
  } else {
    weight_string <- ""
  }

  # header
  if (identical(method, "kruskal")) {
    insight::print_color(sprintf("# Kruskal-Wallis test%s\n\n", weight_string), "blue")
  } else {
    insight::print_color(sprintf("# Friedman test%s\n\n", weight_string), "blue")
  }

  # data info
  insight::print_color(
    sprintf(
      "  Data: %s (%i groups, n = %s)\n",
      x$data, length(n_groups), datawizard::text_concatenate(n_groups)
    ), "cyan"
  )

  stat_symbol <- .format_symbols("Chi2")
  cat(sprintf(
    "\n  %s = %.3f, df = %i, %s\n\n",
    stat_symbol, x$Chi2, round(x$df), insight::format_p(x$p)
  ))
}
