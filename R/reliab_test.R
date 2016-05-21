#' @title Performs a reliability test on an item scale
#' @name reliab_test
#' @description This function calculates the item discriminations (corrected item-total
#'                correlations for each item of \code{x} with the remaining items) and
#'                the Cronbach's alpha for each item, if it was deleted from the
#'                scale.
#'
#' @seealso \code{\link{cronb}}
#'
#' @param x \code{data.frame} with items (from a scale).
#' @param scale.items Logical, if \code{TRUE}, the data frame's vectors will be scaled. Recommended,
#'          when the variables have different measures / scales.
#' @param digits Amount of digits for Cronbach's Alpha and correlation values in
#'          returned data frame.
#' @return A data frame with the corrected item-total correlations (item discrimination)
#'           and Cronbach's alpha (if item deleted) for each item of the scale, or
#'           \code{NULL} if data frame had too less columns.
#'
#' @note This function is similar to a basic reliability test in SPSS. The correlations in
#'         the Item-Total-Statistic are a computed correlation of each item against the sum
#'         of the remaining items (which are thus treated as one item).
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
#' x <- data.frame(efc[, c(start:end)])
#' colnames(x) <- varlabs[c(start:end)]
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
#' for (i in 1:length(findex)) {
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
#' @export
reliab_test <- function(x,
                        scale.items = FALSE,
                        digits = 3) {
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
    if (scale.items) x <- data.frame(scale(x, center = TRUE, scale = TRUE))

    # init vars
    totalCorr <- c()
    cronbachDeleted <- c()

    # iterate all items
    for (i in 1:ncol(x)) {
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
    ret.df <- data.frame(cbind(round(cronbachDeleted, digits),
                               round(totalCorr, digits)))

    # set names of data frame
    colnames(ret.df) <- c("Cronbach's &alpha; if item deleted", "Item discrimination")
    rownames(ret.df) <- df.names
  } else {
    warning("Data frame needs at least three columns for reliability-test.", call. = F)
    ret.df <- NULL
  }
  return(ret.df)
}
