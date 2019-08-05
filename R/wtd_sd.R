#' @title Weighted statistics for tests and variables
#' @name wtd_sd
#' @description \strong{Weighted statistics for variables}
#'   \cr \cr
#'   \code{wtd_sd()}, \code{wtd_se()}, \code{wtd_mean()} and \code{wtd_median()}
#'   compute weighted standard deviation, standard error, mean or median for a
#'   variable or for all variables of a data frame. \code{svy_md()} computes the
#'   median for a variable in a survey-design (see \code{\link[survey]{svydesign}}).
#'   \code{wtd_cor()} computes a weighted correlation for a two-sided alternative
#'   hypothesis.
#'   \cr \cr
#'   \strong{Weighted tests}
#'   \cr \cr
#'   \code{wtd_ttest()} computes a weighted t-test, while \code{wtd_mwu()}
#'   computes a weighted Mann-Whitney-U test or a Kruskal-Wallis test
#'   (for more than two groups). \code{wtd_chisqtest()} computes a weighted
#'   Chi-squared test for contigency tables.
#'
#' @param x (Numeric) vector or a data frame. For \code{svy_md()}, \code{wtd_ttest()},
#'    \code{wtd_mwu()} and \code{wtd_chisqtest()} the bare (unquoted) variable
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
#' @param ... For \code{wtd_ttest()} and \code{wtd_mwu()}, currently not used.
#'   For \code{wtd_chisqtest()}, further arguments passed down to
#'   \code{\link[stats]{chisq.test}}.
#'
#' @inheritParams svyglm.nb
#' @inheritParams grpmean
#'
#' @return The weighted (test) statistic.
#'
#' @note \code{wtd_chisq()} is a convenient wrapper for \code{\link{xtab_statistics}}.
#'    For a weighted one-way Anova, use \code{grpmean()} with
#'    \code{weights}-argument.
#'    \cr \cr
#'    \code{wtd_ttest()} assumes unequal variance between the two groups.
#'
#' @examples
#' # weighted sd and se ----
#'
#' wtd_sd(rnorm(n = 100, mean = 3), runif(n = 100))
#'
#' data(efc)
#' wtd_sd(efc[, 1:3], runif(n = nrow(efc)))
#' wtd_se(efc[, 1:3], runif(n = nrow(efc)))
#'
#' # svy_md ----
#'
#' # median for variables from weighted survey designs
#' library(survey)
#' data(nhanes_sample)
#'
#' des <- svydesign(
#'   id = ~SDMVPSU,
#'   strat = ~SDMVSTRA,
#'   weights = ~WTINT2YR,
#'   nest = TRUE,
#'   data = nhanes_sample
#' )
#'
#' svy_md(total, des)
#' svy_md("total", des)
#'
#' # weighted t-test ----
#'
#' efc$weight <- abs(rnorm(nrow(efc), 1, .3))
#' wtd_ttest(efc, e17age, weights = weight)
#' wtd_ttest(efc, e17age, c160age, weights = weight)
#' wtd_ttest(e17age ~ e16sex + weight, efc)
#'
#' # weighted Mann-Whitney-U-test ----
#'
#' wtd_mwu(c12hour ~ c161sex + weight, efc)
#'
#' # weighted Chi-squared-test ----
#'
#' wtd_chisqtest(efc, c161sex, e16sex, weights = weight, correct = FALSE)
#' wtd_chisqtest(c172code ~ c161sex + weight, efc)
#'
#' @export
wtd_sd <- function(x, weights = NULL) {
  UseMethod("wtd_sd")
}


#' @export
wtd_sd.data.frame <- function(x, weights = NULL) {
  sd_result <- purrr::map_dbl(x, ~ sqrt(wtd_var(.x, weights)))
  names(sd_result) <- colnames(x)

  sd_result
}

#' @export
wtd_sd.matrix <- function(x, weights = NULL) {
  sd_result <- purrr::map_dbl(x, ~ sqrt(wtd_var(.x, weights)))
  names(sd_result) <- colnames(x)

  sd_result
}

#' @export
wtd_sd.default <- function(x, weights = NULL) {
  sqrt(wtd_var(x, weights))
}
