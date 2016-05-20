#' @title Cronbach's Alpha for a matrix or data frame
#' @name cronb
#' @description This function calculates the Cronbach's alpha value
#'                of a data frame or matrix.
#'
#' @seealso \code{\link{reliab_test}}
#'
#' @param data \code{data.frame} or matrix with more than 2 columns.
#' @return The Cronbach's alpha value for \code{data}.
#'
#' @note See 'Examples' from \code{\link[sjPlot]{sjp.pca}} and \code{\link[sjPlot]{sjt.pca}}.
#'
#' @importFrom stats na.omit var
#' @export
cronb <- function(data) {
  .data <- stats::na.omit(data)
  if (is.null(ncol(.data)) || ncol(.data) < 2) {
    warning("Too less columns in `data` to compute Cronbach's Alpha.", call. = F)
    return(NULL)
  }
  return(dim(.data)[2] / (dim(.data)[2] - 1) * (1 - sum(apply(.data, 2, var)) / stats::var(rowSums(.data))))
}
