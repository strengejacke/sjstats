#' @title Weighted statistics for tests
#' @name weighted_mannwhitney
#' @description \code{weighted_ttest()} computes a weighted t-test, while \code{weighted_mannwhitney()}
#'   computes a weighted Mann-Whitney-U test or a Kruskal-Wallis test
#'   (for more than two groups). `weighted_ranktest()` is an alias for `weighted_mannwhitney()`.
#'   \code{weighted_chisqtest()} computes a weighted Chi-squared test for contingency tables.
#'
#' @param x (Numeric) vector or a data frame. For \code{survey_median()}, \code{weighted_ttest()},
#'    \code{weighted_mannwhitney()} and \code{weighted_chisqtest()} the bare (unquoted) variable
#'    name, or a character vector with the variable name.
#' @param weights Bare (unquoted) variable name, or a character vector with
#'    the variable name of the numeric vector of weights. If \code{weights = NULL},
#'    unweighted statistic is reported.
#' @param data A data frame.
#' @param formula A formula of the form \code{lhs ~ rhs1 + rhs2} where \code{lhs} is a
#'    numeric variable giving the data values and \code{rhs1} a factor with two
#'    levels giving the corresponding groups and \code{rhs2} a variable with weights.
#' @param y Optional, bare (unquoted) variable name, or a character vector with
#'    the variable name.
#' @param grp Bare (unquoted) name of the cross-classifying variable, where
#'    \code{x} is grouped into the categories represented by \code{grp},
#'    or a character vector with the variable name.
#' @param mu A number indicating the true value of the mean (or difference in
#'   means if you are performing a two sample test).
#' @param ci.lvl Confidence level of the interval.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}. You can specify just the initial letter.
#' @param paired Logical, whether to compute a paired t-test.
#' @param ... For \code{weighted_ttest()} and \code{weighted_mannwhitney()}, currently not used.
#'   For \code{weighted_chisqtest()}, further arguments passed down to
#'   \code{\link[stats]{chisq.test}}.
#'
#' @return The weighted (test) statistic.
#'
#' @note \code{weighted_chisq()} is a convenient wrapper for \code{\link{crosstable_statistics}}.
#'    For a weighted one-way Anova, use \code{means_by_group()} with
#'    \code{weights}-argument.
#'    \cr \cr
#'    \code{weighted_ttest()} assumes unequal variance between the two groups.
#'
#' @examples
#' # weighted sd and se ----
#' weighted_sd(rnorm(n = 100, mean = 3), runif(n = 100))
#'
#' data(efc)
#' weighted_sd(efc[, 1:3], runif(n = nrow(efc)))
#' weighted_se(efc[, 1:3], runif(n = nrow(efc)))
#'
#' # survey_median ----
#' # median for variables from weighted survey designs
#' if (require("survey")) {
#'   data(nhanes_sample)
#'
#'   des <- svydesign(
#'     id = ~SDMVPSU,
#'     strat = ~SDMVSTRA,
#'     weights = ~WTINT2YR,
#'     nest = TRUE,
#'     data = nhanes_sample
#'   )
#'
#'   survey_median(total, des)
#'   survey_median("total", des)
#' }
#'
#' # weighted t-test ----
#' efc$weight <- abs(rnorm(nrow(efc), 1, .3))
#' weighted_ttest(efc, e17age, weights = weight)
#' weighted_ttest(efc, e17age, c160age, weights = weight)
#' weighted_ttest(e17age ~ e16sex + weight, efc)
#'
#' # weighted Mann-Whitney-U-test ----
#' weighted_mannwhitney(c12hour ~ c161sex + weight, efc)
#'
#' # weighted Chi-squared-test ----
#' weighted_chisqtest(efc, c161sex, e16sex, weights = weight, correct = FALSE)
#' weighted_chisqtest(c172code ~ c161sex + weight, efc)
#'
#' # weighted Chi-squared-test for given probabilities ----
#' weighted_chisqtest(c172code ~ weight, efc, p = c(.33, .33, .34))
#' @export
weighted_mannwhitney <- function(data, ...) {
  UseMethod("weighted_mannwhitney")
}


#' @importFrom dplyr select
#' @rdname weighted_mannwhitney
#' @export
weighted_mannwhitney.default <- function(data, x, grp, weights, ...) {
  x.name <- deparse(substitute(x))
  g.name <- deparse(substitute(grp))
  w.name <- deparse(substitute(weights))

  # create string with variable names
  vars <- c(x.name, g.name, w.name)

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  weighted_mannwhitney_helper(dat)
}


#' @importFrom dplyr select
#' @rdname weighted_mannwhitney
#' @export
weighted_mannwhitney.formula <- function(formula, data, ...) {
  vars <- all.vars(formula)

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  weighted_mannwhitney_helper(dat)
}

weighted_mannwhitney_helper <- function(dat, vars) {
  # check if pkg survey is available
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package `survey` needed to for this function to work. Please install it.", call. = FALSE)
  }

  x.name <- colnames(dat)[1]
  group.name <- colnames(dat)[2]

  colnames(dat) <- c("x", "g", "w")

  if (dplyr::n_distinct(dat$g, na.rm = TRUE) > 2) {
    m <- "Weighted Kruskal-Wallis test"
    method <- "KruskalWallis"
  } else {
    m <- "Weighted Mann-Whitney-U test"
    method <- "wilcoxon"
  }

  design <- survey::svydesign(ids = ~0, data = dat, weights = ~w)
  mw <- survey::svyranktest(formula = x ~ g, design, test = method)

  attr(mw, "x.name") <- x.name
  attr(mw, "group.name") <- group.name
  class(mw) <- c("sj_wmwu", "list")

  mw$method <- m

  mw
}


#' @rdname weighted_mannwhitney
#' @export
weighted_ranktest <- weighted_mannwhitney.formula
