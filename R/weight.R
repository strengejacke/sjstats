#' @title Weight a variable
#' @name weight
#' @description These functions weight the variable \code{x} by
#'                a specific vector of \code{weights}.
#'
#' @param x (Unweighted) variable.
#' @param weights Vector with same length as \code{x}, which
#'          contains weight factors. Each value of \code{x} has a
#'          specific assigned weight in \code{weights}.
#' @param digits Numeric value indicating the number of decimal places to be
#'          used for rounding the weighted values. By default, this value is
#'          \code{0}, i.e. the returned values are integer values.
#'
#' @return The weighted \code{x}.
#'
#' @details \code{weight2()} sums up all \code{weights} values of the associated
#'            categories of \code{x}, whereas \code{weight()} uses a
#'            \code{\link[stats]{xtabs}} formula to weight cases. Thus, \code{weight()}
#'            may return a vector of different length than \code{x}.
#'
#' @note The values of the returned vector are in sorted order, whereas the values'
#'        order of the original \code{x} may be spread randomly. Hence, \code{x} can't be
#'        used, for instance, for further cross tabulation. In case you want to have
#'        weighted contingency tables or (grouped) box plots etc., use the \code{weightBy}
#'        argument of most functions.
#'
#' @examples
#' v <- sample(1:4, 20, TRUE)
#' table(v)
#' w <- abs(rnorm(20))
#' table(weight(v, w))
#' table(weight2(v, w))
#'
#' set.seed(1)
#' x <- sample(letters[1:5], size = 20, replace = TRUE)
#' w <- runif(n = 20)
#'
#' table(x)
#' table(weight(x, w))
#'
#' @importFrom stats na.pass xtabs
#' @importFrom sjlabelled as_numeric
#' @export
weight <- function(x, weights, digits = 0) {
  # remember if x is numeric
  x.is.num <- is.numeric(x)

  # init values
  weightedvar <- c()

  wtab <- round(stats::xtabs(weights ~ x,
                             data = data.frame(weights = weights, x = x),
                             na.action = stats::na.pass,
                             exclude = NULL),
                digits = digits)

  # iterate all table values
  for (w in seq_len(length(wtab))) {
    # retrieve count of each table cell
    w_count <- wtab[[w]]

    # retrieve "cell name" which is identical to the variable value
    # first check whether values are numeric or not
    nval_ <- suppressWarnings(as.numeric(names(wtab[w])))

    # if value is not numeric, use as is
    if (is.na(nval_))
      w_value <- names(wtab[w])
    else
      # else, use numeric value
      w_value <- nval_

    # append variable value, repeating it "w_count" times.
    weightedvar <- c(weightedvar, rep(w_value, w_count))
  }

  # if we have NA values, weighted var is coerced to character.
  # coerce back to numeric then here
  if (!is.numeric(weightedvar) && x.is.num) weightedvar <- sjlabelled::as_numeric(weightedvar)

  # return result
  weightedvar
}

#' @rdname weight
#' @export
weight2 <- function(x, weights) {
  items <- unique(x)
  newvar <- c()

  for (i in seq_len(length(items))) {
    newcount <- round(sum(weights[which(x == items[i])]))
    newvar <- c(newvar, rep(items[i], newcount))
  }

  newvar
}


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


#' @rdname wtd_sd
#' @export
wtd_mean <- function(x, weights = NULL) {
  UseMethod("wtd_mean")
}

#' @importFrom stats weighted.mean
#' @export
wtd_mean.default <- function(x, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  stats::weighted.mean(x, w = weights, na.rm = TRUE)
}

#' @importFrom stats weighted.mean
#' @importFrom purrr map_dbl
#' @importFrom dplyr select_if
#' @export
wtd_mean.data.frame <- function(x, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  dplyr::select_if(x, is.numeric) %>%
    purrr::map_dbl(~ weighted.mean(.x, w = weights))
}


#' @rdname wtd_sd
#' @export
wtd_se <- function(x, weights = NULL) {
  UseMethod("wtd_se")
}


#' @export
wtd_se.data.frame <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ wtd_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
wtd_se.matrix <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ wtd_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
wtd_se.default <- function(x, weights = NULL) {
  wtd_se_helper(x, weights)
}

wtd_se_helper <- function(x, weights) {
  if (is.null(weights)) weights <- rep(1, length(x))
  sqrt(wtd_var(x, weights) / length(stats::na.omit(x)))
}


#' @rdname wtd_sd
#' @export
wtd_median <- function(x, weights = NULL) {
  UseMethod("wtd_median")
}

#' @export
wtd_median.default <- function(x, weights = NULL) {
  wtd_md_helper(x, w = weights, p = 0.5)
}

#' @importFrom purrr map_dbl
#' @importFrom dplyr select_if
#' @export
wtd_median.data.frame <- function(x, weights = NULL) {
  dplyr::select_if(x, is.numeric) %>%
    purrr::map_dbl(~ wtd_md_helper(.x, w = weights, p = 0.5))
}

wtd_md_helper <- function(x, w, p = .5) {
  if (is.null(w)) w <- rep(1, length(x))

  x[is.na(w)] <- NA
  w[is.na(x)] <- NA

  w <- na.omit(w)
  x <- na.omit(x)

  order <- order(x)
  x <- x[order]
  w <- w[order]

  rw <- cumsum(w) / sum(w)
  md.values <- min(which(rw >= p))

  if (rw[md.values] == p)
    q <- mean(x[md.values:(md.values + 1)])
  else
    q <- x[md.values]

  q
}


#' @rdname wtd_sd
#' @importFrom stats as.formula
#' @export
svy_md <- function(x, design) {
  # check if pkg survey is available
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package `survey` needed to for this function to work. Please install it.", call. = FALSE)
  }

  # deparse
  v <- stats::as.formula(paste("~", as.character(substitute(x))))

  as.vector(
    survey::svyquantile(
      v,
      design = design,
      quantiles = 0.5,
      ci = FALSE,
      na.rm = TRUE
    )
  )
}


wtd_var <- function(x, w) {
  if (is.null(w)) w <- rep(1, length(x))

  x[is.na(w)] <- NA
  w[is.na(x)] <- NA

  w <- na.omit(w)
  x <- na.omit(x)

  xbar <- sum(w * x) / sum(w)
  sum(w * ((x - xbar)^2)) / (sum(w) - 1)
}
