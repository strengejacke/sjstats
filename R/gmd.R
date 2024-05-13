#' @title Gini's Mean Difference
#' @name gmd
#' @description `gmd()` computes Gini's mean difference for a numeric vector
#'   or for all numeric vectors in a data frame.
#'
#' @param x A vector or data frame.
#' @param select Optional, names of variables as character vector that should be
#' selected for further processing. Required, if `x` is a data frame (and no vector)
#' and only selected variables from `x` should be processed.
#'
#' @return For numeric vectors, Gini's mean difference. For non-numeric vectors
#' or vectors of length < 2, returns `NA`.
#'
#' @note Gini's mean difference is defined as the mean absolute difference between
#' any two distinct elements of a vector. Missing values from `x` are silently
#' removed.
#'
#' @references David HA. Gini's mean difference rediscovered. Biometrika 1968(55): 573-575
#'
#' @examples
#' data(efc)
#' gmd(efc$e17age)
#' gmd(efc, c("e17age", "c160age", "c12hour"))
#'
#' @export
gmd <- function(x, select = NULL) {
  if (is.data.frame(x)) {
    do.call(rbind, lapply(select, function(i) {
      data.frame(
        variable = i,
        gmd = gmd_helper(x[[i]])
      )
    }))
  } else {
    gmd_helper(x)
  }
}


gmd_helper <- function(x) {
  if (!is.numeric(x)) return(NA)

  x <- stats::na.omit(x)
  n <- length(x)

  if (n < 2) return(NA)

  w <- 4 * ((1:n) - (n - 1) / 2) / n / (n - 1)
  sum(w * sort(x - mean(x)))
}
