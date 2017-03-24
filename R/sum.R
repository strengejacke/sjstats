#' @title Sum, mean and median for vectors
#' @name mn
#'
#' @description \code{mn()}, \code{md()} and \code{sm()} calculate the mean,
#'        median or sum of values in a vector, but have set argument \code{na.rm}
#'        to \code{TRUE} by default.
#'
#' @param x A vector.
#' @param na.rm Logical, whether to remove NA values from \code{x} before computing
#'          mean, median or sum.
#'
#' @return The mean, median or sum of \code{x}.
#'
#' @examples
#' data(efc)
#' md(efc$neg_c_7)
#' mn(efc$neg_c_7)
#' mean(efc$neg_c_7)
#' sm(efc$neg_c_7 > 15)
#'
#' @importFrom stats median
#' @export
mn <- function(x, na.rm = TRUE) mean(x = x, na.rm = na.rm)

#' @rdname mn
#' @export
sm <- function(x, na.rm = TRUE) sum(x = x, na.rm = na.rm)


#' @rdname mn
#' @export
md <- function(x, na.rm = TRUE) stats::median(x = x, na.rm = na.rm)
