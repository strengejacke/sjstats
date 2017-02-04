#' @title Check internal consistency of a test or questionnaire
#' @name reliab_test
#'
#' @description These function compute various measures of internal consistencies
#'                for tests or item-scales of questionnaires.
#'
#' @param x Depending on the function, \code{x} may be a \code{matrix} as
#'          returned by the \code{\link{cor}}-function, or a data frame
#'          with items (e.g. from a test or questionnaire).
#' @param scale.items Logical, if \code{TRUE}, the data frame's vectors will be scaled. Recommended,
#'          when the variables have different measures / scales.
#' @param digits Amount of digits for returned values.
#' @param cor.method Correlation computation method. May be one of
#'          \code{"spearman"} (default), \code{"pearson"} or \code{"kendall"}.
#'          You may use initial letter only.
#'
#' @return \describe{
#'            \item{\code{reliab_test()}}{
#'              A data frame with the corrected item-total correlations (item discrimination,
#'              column \code{item.discr}) and Cronbach's alpha (if item deleted, column
#'              \code{alpha.if.deleted}) for each item of the scale, or \code{NULL}
#'              if data frame had too less columns.
#'            }
#'            \item{\code{split_half()}}{
#'              A list with two values: the split-half reliability \code{splithalf} and
#'              the Spearman-Brown corrected split-half reliability \code{spearmanbrown}.
#'            }
#'            \item{\code{cronb()}}{
#'              The Cronbach's Alpha value for \code{x}.
#'            }
#'            \item{\code{mic()}}{
#'              The mean inter-item-correlation value for \code{x}.
#'            }
#'          }
#'
#' @details \describe{
#'            \item{\code{reliab_test()}}{
#'              This function calculates the item discriminations (corrected item-total
#'              correlations for each item of \code{x} with the remaining items) and
#'              the Cronbach's alpha for each item, if it was deleted from the scale.
#'            }
#'            \item{\code{split_half()}}{
#'              This function calculates the split-half reliability for items in
#'              the data frame \code{x}, including the Spearman-Brown adjustment.
#'              Splitting is done by selecting odd versus even columns in \code{x}.
#'            }
#'            \item{\code{cronb()}}{
#'              The Cronbach's Alpha value for \code{x}.
#'            }
#'            \item{\code{mic()}}{
#'              This function calculates a mean inter-item-correlation, i.e.
#'              a correlation matrix of \code{x} will be computed (unless
#'              \code{x} is already a matrix as returned by the
#'              \code{\link{cor}}-function) and the mean
#'              of the sum of all item's correlation values is returned.
#'              Requires either a data frame or a computed \code{\link{cor}}-object.
#'            }
#'          }
#'
#' @note \code{reliab_test()} is similar to a basic reliability test
#'         in SPSS. The correlations in the Item-Total-Statistic are a computed
#'         correlation of each item against the sum of the remaining items
#'         (which are thus treated as one item).
#'         \cr \cr
#'         For \code{split_half()} and \code{cronb()}, a value closer  to 1 indicates
#'         greater internal consistency.
#'         \cr \cr
#'         For the mean inter-item-correlation:
#'         \dQuote{Ideally, the average inter-item correlation for a set of
#'          items should be between .20 and .40, suggesting that while the
#'          items are reasonably homogenous, they do contain sufficiently
#'          unique variance so as to not be isomorphic with each other.
#'          When values are lower than .20, then the items may not be
#'          representative of the same content domain. If values are higher than
#'          .40, the items may be only capturing a small bandwidth of the construct.}
#'          \cite{(Piedmont 2014)}
#'
#' @references Spearman C. 1910. Correlation calculated from faulty data. British Journal of Psychology (3): 271–295. \doi{10.1111/j.2044-8295.1910.tb00206.x}
#'             \cr \cr
#'             Brown W. 1910. Some experimental results in the correlation of mental abilities. British Journal of Psychology (3): 296–322. \doi{10.1111/j.2044-8295.1910.tb00207.x}
#'             \cr \cr
#'             Piedmont RL. 2014. Inter-item Correlations. In: Michalos AC (eds) Encyclopedia of Quality of Life and Well-Being Research. Dordrecht: Springer, 3303-3304. \doi{10.1007/978-94-007-0753-5_1493}
#'
#' @examples
#' library(sjmisc)
#' # Data from the EUROFAMCARE sample dataset
#' data(efc)
#'
#' # retrieve variable and value labels
#' varlabs <- get_label(efc)
#'
#' # recveive first item of COPE-index scale
#' start <- which(colnames(efc) == "c82cop1")
#' # recveive last item of COPE-index scale
#' end <- which(colnames(efc) == "c90cop9")
#'
#' # create data frame with COPE-index scale
#' x <- efc[, c(start:end)]
#' colnames(x) <- varlabs[c(start:end)]
#'
#' # reliability tests
#' reliab_test(x)
#'
#' # split-half-reliability
#' split_half(x)
#'
#' # cronbach's alpha
#' cronb(x)
#'
#' # mean inter-item-correlation
#' mic(x)
#'
#' \dontrun{
#' library(sjPlot)
#' sjt.df(reliab_test(x), describe = FALSE, show.cmmn.row = TRUE,
#'        string.cmmn = sprintf("Cronbach's &alpha;=%.2f", cronb(x)))
#'
#' # Compute PCA on Cope-Index, and perform a
#' # reliability check on each extracted factor.
#' factors <- sjt.pca(x)$factor.index
#' findex <- sort(unique(factors))
#' library(sjPlot)
#' for (i in seq_len(length(findex))) {
#'  rel.df <- subset(x, select = which(factors == findex[i]))
#'  if (ncol(rel.df) >= 3) {
#'    sjt.df(reliab_test(rel.df), describe = FALSE, show.cmmn.row = TRUE,
#'           use.viewer = FALSE, title = "Item-Total-Statistic",
#'           string.cmmn = sprintf("Scale's overall Cronbach's &alpha;=%.2f",
#'                                 cronb(rel.df)))
#'    }
#'  }}
#'
#' @importFrom stats cor
#' @importFrom tibble tibble
#' @importFrom sjmisc std
#' @export
reliab_test <- function(x, scale.items = FALSE, digits = 3) {
  # check param
  if (!is.matrix(x) && !is.data.frame(x)) {
    warning("`x` needs to be a data frame or matrix.", call. = F)
    return(NULL)
  }

  # remove missings, so correlation works
  x <- stats::na.omit(x)

  # remember item (column) names for return value
  # return value gets column names of initial data frame
  df.names <- colnames(x)

  # check for minimum amount of columns
  # can't be less than 3, because the reliability
  # test checks for Cronbach's alpha if a specific
  # item is deleted. If data frame has only two columns
  # and one is deleted, Cronbach's alpha cannot be calculated.
  if (ncol(x) > 2) {
    # Check whether items should be scaled. Needed,
    # when items have different measures / scales
    if (scale.items) x <- sjmisc::std(x)

    # init vars
    totalCorr <- c()
    cronbachDeleted <- c()

    # iterate all items
    for (i in seq_len(ncol(x))) {
      # create subset with all items except current one
      # (current item "deleted")
      sub.df <- subset(x, select = c(-i))

      # calculate cronbach-if-deleted
      cronbachDeleted <- c(cronbachDeleted, cronb(sub.df))

      # calculate corrected total-item correlation
      totalCorr <- c(totalCorr, stats::cor(x[, i],
                                           apply(sub.df, 1, sum),
                                           use = "pairwise.complete.obs"))
    }

    # create return value
    ret.df <- tibble::tibble(term = df.names,
                             alpha.if.deleted = round(cronbachDeleted, digits),
                             item.discr = round(totalCorr, digits))
  } else {
    warning("Data frame needs at least three columns for reliability-test.", call. = F)
    ret.df <- NULL
  }
  return(ret.df)
}



#' @rdname reliab_test
#' @importFrom stats cor
#' @export
split_half <- function(x, digits = 3) {
  # Calculating total score for even items
  score_e <- rowMeans(x[, c(TRUE, FALSE)], na.rm = TRUE)
  # Calculating total score for odd items
  score_o <- rowMeans(x[, c(FALSE, TRUE)], na.rm = TRUE)

  # Correlating scores from even and odd items
  shr <- stats::cor(score_e, score_o, use = "complete.obs")

  # Adjusting with the Spearman-Brown prophecy formula
  sb.shr <- (2 * shr) / (1 + shr)

  structure(class = "sj_splithalf",
            list(splithalf = shr, spearmanbrown = sb.shr))
}



#' @rdname reliab_test
#' @importFrom stats na.omit var
#' @export
cronb <- function(x) {
  # remove missings
  .data <- stats::na.omit(x)

  # we need at least two columns for Cronach's Alpha
  if (is.null(ncol(.data)) || ncol(.data) < 2) {
    warning("Too less columns in `x` to compute Cronbach's Alpha.", call. = F)
    return(NULL)
  }

  # Compute Cronb. Alpha
  return(dim(.data)[2] / (dim(.data)[2] - 1) * (1 - sum(apply(.data, 2, var)) / stats::var(rowSums(.data))))
}



#' @rdname reliab_test
#' @importFrom stats cor na.omit
#' @export
mic <- function(x, cor.method = c("pearson", "spearman", "kendall")) {
  # Check parameter
  cor.method <- match.arg(cor.method)

  # Mean-interitem-corelation
  if (class(x) == "matrix") {
    corr <- x
  } else {
    x <- stats::na.omit(x)
    corr <- stats::cor(x, method = cor.method)
  }

  # Sum up all correlation values
  meanic <- c()
  for (j in seq_len((ncol(corr) - 1))) {
    # first correlation is always "1" (self-correlation)
    for (i in (j + 1):nrow(corr)) {
      # check for valid bound
      if (i <= nrow(corr) && j <= ncol(corr)) {
        # add up all subsequent values
        meanic <- c(meanic, corr[i, j])
      } else {
        meanic <- c(meanic, "NA")
      }
    }
  }
  return(mean(meanic))
}
