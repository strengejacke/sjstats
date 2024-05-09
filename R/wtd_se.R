#' @title Weighted statistics for tests and variables
#' @name weighted_se
#' @description \strong{Weighted statistics for variables}
#'   \cr \cr
#'   \code{weighted_sd()}, \code{weighted_se()}, \code{weighted_mean()} and \code{weighted_median()}
#'   compute weighted standard deviation, standard error, mean or median for a
#'   variable or for all variables of a data frame. \code{survey_median()} computes the
#'   median for a variable in a survey-design (see \code{\link[survey]{svydesign}}).
#'   \code{weighted_correlation()} computes a weighted correlation for a two-sided alternative
#'   hypothesis.
#'   \cr \cr
#'   \strong{Weighted tests}
#'   \cr \cr
#'   \code{weighted_ttest()} computes a weighted t-test, while \code{weighted_mannwhitney()}
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
#' @inheritParams svyglm.nb
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
#' @export
weighted_se <- function(x, weights = NULL) {
  UseMethod("weighted_se")
}


#' @export
weighted_se.data.frame <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ weighted_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
weighted_se.matrix <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ weighted_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
weighted_se.default <- function(x, weights = NULL) {
  weighted_se_helper(x, weights)
}

weighted_se_helper <- function(x, weights) {
  if (is.null(weights)) weights <- rep(1, length(x))
  sqrt(weighted_variance(x, weights) / length(stats::na.omit(x)))
}
