#' @title Measures of association for contingency tables
#' @name crosstable_statistics
#'
#' @description This function calculates various measure of association for
#'              contingency tables and returns the statistic and p-value.
#'              Supported measures are Cramer's V, Phi, Spearman's rho,
#'              Kendall's tau and Pearson's r.
#'
#' @param data A data frame or a table object. If a table object, `x1` and
#' `x2` will be ignored. For Kendall's _tau_, Spearman's _rho_ or Pearson's
#' product moment correlation coefficient, `data` needs to be a data frame.
#' If `x1` and `x2` are not specified, the first two columns of the data
#' frames are used as variables to compute the crosstab.
#' @param formula A formula of the form `lhs ~ rhs` where `lhs` is a
#'  numeric variable giving the data values and `rhs` a factor giving the
#'  corresponding groups.
#' @param tab A [`table()`] or [`ftable()`]. Tables of class [`xtabs()`] and
#' other will be coerced to `ftable` objects.
#' @param x1 Name of first variable that should be used to compute the
#' contingency table. If `data` is a table object, this argument will be
#' irgnored.
#' @param x2 Name of second variable that should be used to compute the
#' contingency table. If `data` is a table object, this argument will be
#' irgnored.
#' @param statistics Name of measure of association that should be computed. May
#' be one of `"auto"`, `"cramer"`, `"phi"`, `"spearman"`, `"kendall"`,
#' `"pearson"` or `"fisher"`. See 'Details'.
#' @param ci.lvl Scalar between 0 and 1. If not `NULL`, returns a data
#' frame including lower and upper confidence intervals.
#' @param weights Name of variable in `x` that indicated the vector of weights
#' that will be applied to weight all observations. Default is `NULL`, so no
#' weights are used.
#' @param ... Other arguments, passed down to the statistic functions
#' [`chisq.test()`], [`fisher.test()`] or [`cor.test()`].
#'
#' @inheritParams bootstrap
#' @inheritParams boot_ci
#'
#' @return For [`phi()`], the table's Phi value. For [`cramers_v()]`, the
#' table's Cramer's V.
#'
#' For `crosstable_statistics()`, a list with following components:
#'
#' - `estimate`: the value of the estimated measure of association.
#' - `p.value`: the p-value for the test.
#' - `statistic`: the value of the test statistic.
#' - `stat.name`: the name of the test statistic.
#' - `stat.html`: if applicable, the name of the test statistic, in HTML-format.
#' - `df`: the degrees of freedom for the contingency table.
#' - `method`: character string indicating the name of the measure of association.
#' - `method.html`: if applicable, the name of the measure of association, in HTML-format.
#' - `method.short`: the short form of association measure, equals the `statistics`-argument.
#' - `fisher`: logical, if Fisher's exact test was used to calculate the p-value.
#'
#' @details The p-value for Cramer's V and the Phi coefficient are based
#' on `chisq.test()`. If any expected value of a table cell is smaller than 5,
#' or smaller than 10 and the df is 1, then `fisher.test()` is used to compute
#' the p-value, unless `statistics = "fisher"`; in this case, the use of
#' `fisher.test()` is forced to compute the p-value. The test statistic is
#' calculated with `cramers_v()` resp. `phi()`.
#'
#' Both test statistic and p-value for Spearman's rho, Kendall's tau and
#' Pearson's r are calculated with `cor.test()`.
#'
#' When `statistics = "auto"`, only Cramer's V or Phi are calculated, based on
#' the dimension of the table (i.e. if the table has more than two rows or
#' columns, Cramer's V is calculated, else Phi).
#'
#' @references Ben-Shachar, M.S., Patil, I., Thériault, R., Wiernik, B.M.,
#' Lüdecke, D. (2023). Phi, Fei, Fo, Fum: Effect Sizes for Categorical Data
#' That Use the Chi‑Squared Statistic. Mathematics, 11, 1982.
#' \doi{10.3390/math11091982}
#'
#' @examples
#' # Phi coefficient for 2x2 tables
#' tab <- table(sample(1:2, 30, TRUE), sample(1:2, 30, TRUE))
#' phi(tab)
#'
#' # Cramer's V for nominal variables with more than 2 categories
#' tab <- table(sample(1:2, 30, TRUE), sample(1:3, 30, TRUE))
#' cramer(tab)
#'
#' # formula notation
#' data(efc)
#' cramer(e16sex ~ c161sex, data = efc)
#'
#' # bootstrapped confidence intervals
#' cramer(e16sex ~ c161sex, data = efc, ci.lvl = .95, n = 100)
#'
#' # 2x2 table, compute Phi automatically
#' crosstable_statistics(efc, e16sex, c161sex)
#'
#' # more dimensions than 2x2, compute Cramer's V automatically
#' crosstable_statistics(efc, c172code, c161sex)
#'
#' # ordinal data, use Kendall's tau
#' crosstable_statistics(efc, e42dep, quol_5, statistics = "kendall")
#'
#' # calcilate Spearman's rho, with continuity correction
#' crosstable_statistics(efc,
#'   e42dep,
#'   quol_5,
#'   statistics = "spearman",
#'   exact = FALSE,
#'   continuity = TRUE
#' )
#' @export
crosstable_statistics <- function(data, x1 = NULL, x2 = NULL, statistics = c("auto", "cramer", "phi", "spearman", "kendall", "pearson", "fisher"), weights = NULL, ...) {
  # match arguments
  statistics <- match.arg(statistics)

  # name for test statistics in HTML
  stat.html <- NULL

  # check if data is a table
  if (is.table(data)) {
    # 'data' is a table - copy to table object
    tab <- data
    # check if statistics are possible to compute
    if (statistics %in% c("spearman", "kendall", "pearson")) {
      stop(
        sprintf(
          "Need arguments `data`, `x1` and `x2` to compute %s-statistics.",
          statistics
        ),
        call. = FALSE
      )
    }
  } else {
    # evaluate unquoted names
    x1 <- deparse(substitute(x1))
    x2 <- deparse(substitute(x2))
    weights <- deparse(substitute(weights))

    # if names were quotes, remove quotes
    x1 <- gsub("\"", "", x1, fixed = TRUE)
    x2 <- gsub("\"", "", x2, fixed = TRUE)
    weights <- gsub("\"", "", weights, fixed = TRUE)

    if (insight::is_empty_object(weights) || weights == "NULL")
      weights <- NULL
    else
      weights <- data[[weights]]


    # check for "NULL" and get data
    if (x1 != "NULL" && x2 != "NULL")
      data <- data[, c(x1, x2)]
    else
      data <- data[, 1:2]

    if (!is.null(weights)) data <- cbind(data, weights)

    # make table
    if (!is.null(weights)) {
      tab <- as.table(round(stats::xtabs(data[[3]] ~ data[[1]] + data[[2]])))
      class(tab) <- "table"
    } else {
      tab <- table(data)
    }
  }

  # get expected values
  tab.val <- table_values(tab)

  # remember whether fisher's exact test was used or not
  use.fisher <- FALSE

  # select statistics automatically, based on number of rows/columns
  if (statistics %in% c("auto", "cramer", "phi", "fisher")) {
    # get chisq-statistics, for df and p-value
    chsq <- suppressWarnings(stats::chisq.test(tab, ...))
    pv <- chsq$p.value
    test <- chsq$statistic

    # set statistics name
    names(test) <- "Chi-squared"
    stat.html <- "&chi;<sup>2</sup>"

    # check row/columns
    if ((nrow(tab) > 2 || ncol(tab) > 2 || statistics %in% c("cramer", "fisher")) && statistics != "phi") {
      # get cramer's V
      s <- cramer(tab)

      # if minimum expected values below 5, compute fisher's exact test
      if (statistics == "fisher" ||
          min(tab.val$expected) < 5 ||
          (min(tab.val$expected) < 10 && chsq$parameter == 1)) {
        pv <- stats::fisher.test(tab, simulate.p.value = TRUE, ...)$p.value
        use.fisher <- TRUE
      }

      # set statistics
      statistics <- "cramer"
    } else {
      # get Phi
      s <- phi(tab)

      # if minimum expected values below 5 and df=1, compute fisher's exact test
      if (min(tab.val$expected) < 5 ||
          (min(tab.val$expected) < 10 && chsq$parameter == 1)) {
        pv <- stats::fisher.test(tab, ...)$p.value
        use.fisher <- TRUE
      }

      # set statistics
      statistics <- "phi"
    }
  } else {
    # compute correlation coefficient
    cv <- stats::cor.test(x = data[[1]], y = data[[2]], method = statistics, ...)
    # get statistics and p-value
    s <- cv$estimate
    pv <- cv$p.value
    test <- cv$statistic
    stat.html <- names(test)
  }

  # compute method string
  method <- ifelse(statistics == "kendall", "Kendall's tau",
    ifelse(statistics == "spearman", "Spearman's rho", # nolint
      ifelse(statistics == "pearson", "Pearson's r", # nolint
        ifelse(statistics == "cramer", "Cramer's V", "Phi") # nolint
      )
    )
  )

  # compute method string
  method.html <- ifelse(statistics == "kendall", "Kendall's &tau;",
    ifelse(statistics == "spearman", "Spearman's &rho;", # nolint
      ifelse(statistics == "pearson", "Pearson's r", # nolint
        ifelse(statistics == "cramer", "Cramer's V", "&phi") # nolint
      )
    )
  )

  # return result
  structure(class = "sj_xtab_stat", list(
    estimate = s,
    p.value = pv,
    statistic = test,
    stat.name = names(test),
    stat.html = stat.html,
    df = (nrow(tab) - 1) * (ncol(tab) - 1),
    n_obs = sum(tab, na.rm = TRUE),
    method = method,
    method.html = method.html,
    method.short = statistics,
    fisher = use.fisher
  ))
}

#' @rdname crosstable_statistics
#' @export
xtab_statistics <- crosstable_statistics
