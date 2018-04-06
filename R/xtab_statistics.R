#' @title Measures of association for contingency tables
#' @name xtab_statistics
#'
#' @description This function calculates various measure of association for
#'              contingency tables and returns the statistic and p-value.
#'              Supported measures are Cramer's V, Phi, Spearman's rho,
#'              Kendall's tau and Pearson's r.
#'
#' @param data A data frame or a table object. If a table object, \code{x1} and
#'          \code{x2} will be ignored. For Kendall's \emph{tau}, Spearman's \emph{rho}
#'          or Pearson's product moment correlation coefficient, \code{data} needs
#'          to be a data frame. If \code{x1} and \code{x2} are not specified,
#'          the first two columns of the data frames are used as variables
#'          to compute the crosstab.
#' @param tab A \code{\link{table}} or \code{\link{ftable}}. Tables of class
#'          \code{\link{xtabs}} and other will be coerced to \code{\link{ftable}}
#'          objects.
#' @param x1 Name of first variable that should be used to compute the
#'           contingency table. If \code{data} is a table object, this argument
#'           will be irgnored.
#' @param x2 Name of second variable that should be used to compute the
#'           contingency table. If \code{data} is a table object, this argument
#'           will be irgnored.
#' @param statistics Name of measure of association that should be computed. May
#'          be one of \code{"auto"}, \code{"cramer"}, \code{"phi"}, \code{"spearman"},
#'          \code{"kendall"} or \code{"pearson"}. See 'Details'.
#' @param ... Other arguments, passed down to the statistic functions
#'          \code{\link[stats]{chisq.test}}, \code{\link[stats]{fisher.test}} or
#'          \code{\link[stats]{cor.test}}.
#'
#' @return For \code{phi()}, the table's Phi value. For \code{cramer()}, the
#'         table's Cramer's V.
#'         \cr \cr
#'         For \code{xtab_statistics()}, a list with following components:
#'         \describe{
#'           \item{\code{estimate}}{the value of the estimated measure of association.}
#'           \item{\code{p.value}}{the p-value for the test.}
#'           \item{\code{statistic}}{the value of the test statistic.}
#'           \item{\code{stat.name}}{the name of the test statistic.}
#'           \item{\code{stat.html}}{if applicable, the name of the test statistic, in HTML-format.}
#'           \item{\code{df}}{the degrees of freedom for the contingency table.}
#'           \item{\code{method}}{character string indicating the name of the measure of association.}
#'           \item{\code{method.html}}{if applicable, the name of the measure of association, in HTML-format.}
#'           \item{\code{method.short}}{the short form of association measure, equals the \code{statistics}-aergument.}
#'           \item{\code{fisher}}{logical, if Fisher's exact test was used to calculate the p-value.}
#'         }
#'
#' @details The p-value for Cramer's V and the Phi coefficient are based
#'          on \code{chisq.test()}. If any expected value of a table cell is
#'          smaller than 5, or smaller than 10 and the df is 1, then \code{fisher.test()}
#'          is used to compute the p-value. The test statistic is calculated
#'          with \code{cramer()} resp. \code{phi()}.
#'          \cr \cr
#'          Both test statistic and p-value for Spearman's rho, Kendall's tau
#'          and Pearson's r are calculated with \code{cor.test()}.
#'          \cr \cr
#'          When \code{statistics = "auto"}, only Cramer's V or Phi are calculated,
#'          based on the dimension of the table (i.e. if the table has more than
#'          two rows or columns, Cramer's V is calculated, else Phi).
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
#' data(efc)
#' # 2x2 table, compute Phi automatically
#' xtab_statistics(efc, e16sex, c161sex)
#'
#' # more dimensions than 2x2, compute Cramer's V automatically
#' xtab_statistics(efc, c172code, c161sex)
#'
#' # ordinal data, use Kendall's tau
#' xtab_statistics(efc, e42dep, quol_5, statistics = "kendall")
#'
#' # calcilate Spearman's rho, with continuity correction
#' xtab_statistics(efc,
#'   e42dep,
#'   quol_5,
#'   statistics = "spearman",
#'   exact = FALSE,
#'   continuity = TRUE
#' )
#'
#' @importFrom stats fisher.test chisq.test cor.test ftable
#' @importFrom dplyr case_when
#' @importFrom MASS loglm
#' @export
xtab_statistics <- function(data, x1 = NULL, x2 = NULL, statistics = c("auto", "cramer", "phi", "spearman", "kendall", "pearson"), ...) {
  # match arguments
  statistics <- match.arg(statistics)

  # name for test statistics in HTML
  stat.html <- NULL

  # check if data is a table
  if (!is.table(data)) {
    # evaluate unquoted names
    x1 <- deparse(substitute(x1))
    x2 <- deparse(substitute(x2))

    # if names were quotes, remove quotes
    x1 <- gsub("\"", "", x1, fixed = T)
    x2 <- gsub("\"", "", x2, fixed = T)

    # check for "NULL" and get data
    if (x1 != "NULL" && x2 != "NULL")
      data <- data[, c(x1, x2)]
    else
      data <- data[, 1:2]

    # make simple table
    tab <- table(data)
  } else {
    # 'data' is a table - copy to table object
    tab <- data
    # check if statistics are possible to compute
    if (statistics %in% c("spearman", "kendall", "pearson")) {
      stop(
        sprintf(
          "Need arguments `data`, `x1` and `x2` to compute %s-statistics.",
          statistics
        ),
        call. = F
      )
    }
  }

  # get expected values
  tab.val <- table_values(tab)

  # remember whether fisher's exact test was used or not
  use.fisher <- FALSE

  # select statistics automatically, based on number of rows/columns
  if (statistics %in% c("auto", "cramer", "phi")) {
    # get chisq-statistics, for df and p-value
    chsq <- suppressWarnings(stats::chisq.test(tab, ...))
    pv <- chsq$p.value
    test <- chsq$statistic

    # set statistics name
    names(test) <- "Chi-squared"
    stat.html <- "&chi;<sup>2</sup>"

    # check row/column
    if ((nrow(tab) > 2 || ncol(tab) > 2 || statistics == "cramer") && statistics != "phi") {
      # get cramer's V
      s <- cramer(tab)

      # if minimum expected values below 5, compute fisher's exact test
      if (min(tab.val$expected) < 5 ||
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
  method <- dplyr::case_when(
    statistics == "kendall" ~ "Kendall's tau",
    statistics == "spearman" ~ "Spearman's rho",
    statistics == "pearson" ~ "Pearson's r",
    statistics == "cramer" ~ "Cramer's V",
    statistics == "phi" ~ "Phi"
  )

  # compute method string
  method.html <- dplyr::case_when(
    statistics == "kendall" ~ "Kendall's &tau;",
    statistics == "spearman" ~ "Spearman's &rho;",
    statistics == "pearson" ~ "Pearson's r",
    statistics == "cramer" ~ "Cramer's V",
    statistics == "phi" ~ "&phi;"
  )

  # return result
  structure(class = "sj_xtab_stat", list(
    estimate = s,
    p.value = pv,
    statistic = test,
    stat.name = names(test),
    stat.html = stat.html,
    df = (nrow(tab) - 1) * (ncol(tab) - 1),
    method = method,
    method.html = method.html,
    method.short = statistics,
    fisher = use.fisher
  ))
}
