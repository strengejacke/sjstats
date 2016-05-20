#' @title Mean Inter-Item-Correlation
#' @name mic
#' @description This function calculates a mean inter-item-correlation, i.e.
#'                a correlation matrix of \code{data} will be computed (unless
#'                \code{data} is already a matrix as returned by the
#'                \code{\link{cor}}-function) and the mean
#'                of the sum of all item's correlation values is returned.
#'                Requires either a data frame or a computed \code{\link{cor}}-object.
#'
#' @param data A \code{matrix} as returned by the \code{\link{cor}}-function, or
#'          a \code{data.frame}, where correlations will be calculated.
#' @param cor.method Indicates the correlation computation method. May be one of
#'          \code{"spearman"} (default), \code{"pearson"} or \code{"kendall"}.
#'          You may use initial letter only.
#' @return The value of the computed mean inter-item-correlation.
#'
#' @note \dQuote{Ideally, the average inter-item correlation for a set of
#'          items should be between .20 and .40, suggesting that while the
#'          items are reasonably homogenous, they do contain sufficiently
#'          unique variance so as to not be isomorphic with each other.
#'
#'          When values are lower than .20, then the items may not be
#'          representative of the same content domain. If values are higher than
#'          .40, the items may be only capturing a small bandwidth of the construct.}
#'          \emph{(Piedmont 2014)}
#'
#' @references Piedmont RL (2014) Inter-item Correlations. In: Michalos AC (eds)
#'             Encyclopedia of Quality of Life and Well-Being Research.
#'             Dordrecht: Springer, 3303-3304
#'             \doi{10.1007/978-94-007-0753-5_1493}
#'
#' @examples
#' # Data from the EUROFAMCARE sample dataset
#' data(efc)
#' # recveive first item of COPE-index scale
#' start <- which(colnames(efc) == "c82cop1")
#' # recveive last item of COPE-index scale
#' end <- which(colnames(efc) == "c90cop9")
#' # create data frame with COPE-index scale
#' mydat <- data.frame(efc[, c(start:end)])
#'
#' mic(mydat)
#'
#' @importFrom stats cor na.omit
#' @export
mic <- function(data, cor.method = c("pearson", "spearman", "kendall")) {
  # Check parameter
  cor.method <- match.arg(cor.method)

  # Mean-interitem-corelation
  if (class(data) == "matrix") {
    corr <- data
  } else {
    data <- stats::na.omit(data)
    corr <- stats::cor(data, method = cor.method)
  }

  # Sum up all correlation values
  meanic <- c()
  for (j in 1:(ncol(corr) - 1)) {
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
