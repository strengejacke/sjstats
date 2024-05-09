#' @title Kruskal-Wallis-Test
#' @name kruskal_wallis_test
#' @description This function performs a Kruskal-Wallis rank sum test, to test
#' the null hypothesis that the population median of all of the groups are
#' equal. The alternative is that they differ in at least one.
#'
#' @inheritParams mann_whitney_test
#'
#' @return A data frame with test results.
#'
#' @details The function simply is a wrapper around [`kruskal.test()`]. The
#' weighted version of the Kruskal-Wallis test is based on the `survey` package,
#' using [`survey::svyranktest()`].
#'
#' @examples
#' data(efc)
#' # Kruskal-Wallis-Test for elder's age by education
#' kruskal_wallis_test(efc, "e17age", by = "c172code")
#' @export
kruskal_wallis_test <- function(data,
                                select = NULL,
                                by = NULL,
                                weights = NULL) {
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
    .calculate_kw(dv, grp)
  } else {
    .calculate_weighted_kw(dv, grp, data[[weights]])
  }
}


# Kruskal-Wallis-Test --------------------------------------------

.calculate_kw <- function(dv, grp) {
  # prepare data
  wcdat <- data.frame(dv, grp)
  # perfom wilcox test
  wt <- stats::kruskal.test(dv ~ grp, data = wcdat)
  # number of groups
  n_groups <- vapply(unique(grp), function(g) sum(grp == g, na.rm = TRUE), numeric(1))

  out <- data.frame(
    data = wt$data.name,
    Chi2 = wt$statistic,
    df = wt$parameter,
    p = as.numeric(wt$p.value),
    stringsAsFactors = FALSE
  )

  attr(out, "n_groups") <- n_groups
  attr(out, "method") <- "kruskal"
  attr(out, "weighted") <- FALSE
  class(out) <- c("sj_htest_kw", "data.frame")

  out
}


# Weighted Mann-Whitney-Test for two groups ----------------------------------

.calculate_weighted_kw <- function(dv, grp, weights) {
  # check if pkg survey is available
  insight::check_if_installed("survey")

  dat <- stats::na.omit(data.frame(dv, grp, weights))
  colnames(dat) <- c("x", "g", "w")

  design <- survey::svydesign(ids = ~0, data = dat, weights = ~w)
  result <- survey::svyranktest(formula = x ~ g, design, test = "KruskalWallis")

  # number of groups
  n_groups <- vapply(unique(grp), function(g) {
    sum(dat$w[dat$grp == g], na.rm = TRUE)
  }, numeric(1))

  out <- data.frame(
    data = paste(dv, "by", grp),
    Chi2 = result$statistic,
    df = result$parameter,
    p = as.numeric(result$p.value),
    stringsAsFactors = FALSE
  )

  attr(out, "n_groups") <- n_groups
  attr(out, "method") <- "kruskal"
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

  if (weighted) {
    weight_string <- " (weighted)"
  } else {
    weight_string <- ""
  }

  # header
  insight::print_color(sprintf("# Kruskal-Wallis-Test%s\n\n", weight_string), "blue")

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
