#' @title Weighted statistics for tests and variables
#' @name weighted_se
#' @description
#' `weighted_se()` computes weighted standard errors of a variable or for
#' all variables of a data frame. `survey_median()` computes the median for
#' a variable in a survey-design (see [`survey::svydesign()]`).
#' `weighted_correlation()` computes a weighted correlation for a two-sided
#' alternative hypothesis.
#'
#' @param x (Numeric) vector or a data frame. For `survey_median()` or `weighted_ttest()`,
#' the bare (unquoted) variable name, or a character vector with the variable name.
#' @param weights Bare (unquoted) variable name, or a character vector with
#' the variable name of the numeric vector of weights. If `weights = NULL`,
#' unweighted statistic is reported.
#' @param data A data frame.
#' @param formula A formula of the form `lhs ~ rhs1 + rhs2` where `lhs` is a
#' numeric variable giving the data values and `rhs1` a factor with two
#' levels giving the corresponding groups and `rhs2` a variable with weights.
#' @param y Optional, bare (unquoted) variable name, or a character vector with
#' the variable name.
#' @param ci.lvl Confidence level of the interval.
#' @param ... Currently not used.
#'
#' @inheritParams svyglm.nb
#'
#' @return The weighted (test) statistic.
#'
#' @examplesIf requireNamespace("survey")
#' data(efc)
#' weighted_se(efc$c12hour, abs(runif(n = nrow(efc))))
#'
#' # survey_median ----
#' # median for variables from weighted survey designs
#' data(nhanes_sample)
#'
#' des <- survey::svydesign(
#'   id = ~SDMVPSU,
#'   strat = ~SDMVSTRA,
#'   weights = ~WTINT2YR,
#'   nest = TRUE,
#'   data = nhanes_sample
#' )
#' survey_median(total, des)
#' survey_median("total", des)
#' @export
weighted_se <- function(x, weights = NULL) {
  UseMethod("weighted_se")
}


#' @export
weighted_se.data.frame <- function(x, weights = NULL) {
  se_result <- vapply(x, weighted_se_helper, numeric(1), weights = weights)
  names(se_result) <- colnames(x)

  se_result
}

#' @export
weighted_se.matrix <- function(x, weights = NULL) {
  se_result <- vapply(x, weighted_se_helper, numeric(1), weights = weights)
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
