#' @title Gini's Mean Difference
#' @name gmd
#' @description \code{gmd()} computes Gini's mean difference for a numeric vector
#'   or for all numeric vectors in a data frame.
#'
#' @param x A vector or data frame.
#' @param ... Optional, unquoted names of variables that should be selected for
#'   further processing. Required, if \code{x} is a data frame (and no vector)
#'   and only selected variables from \code{x} should be processed. You may also
#'   use functions like \code{:} or tidyselect's \code{\link[tidyselect]{select_helpers}}.
#'
#' @return For numeric vectors, Gini's mean difference. For non-numeric vectors
#'   or vectors of length < 2, returns \code{NA}.
#'
#' @note Gini's mean difference is defined as the mean absolute difference between
#'   any two distinct elements of a vector. Missing values from \code{x} are
#'   silently removed.
#'
#' @references David HA. Gini's mean difference rediscovered. Biometrika 1968(55): 573–575
#'
#' @examples
#' data(efc)
#' gmd(efc$e17age)
#' gmd(efc, e17age, c160age, c12hour)
#'
#' @importFrom dplyr quos select
#' @importFrom purrr map_df
#' @importFrom sjmisc is_empty
#' @export
gmd <- function(x, ...) {
  # evaluate dots
  qs <- dplyr::quos(...)
  if (!sjmisc::is_empty(qs)) x <- suppressMessages(dplyr::select(x, !!!qs))

  if (is.data.frame(x))
    purrr::map_df(x, gmd_helper)
  else
    gmd_helper(x)
}


#' @importFrom stats na.omit
gmd_helper <- function(x) {
  if (!is.numeric(x)) return(NA)

  x <- stats::na.omit(x)
  n <- length(x)

  if (n < 2) return(NA)

  w <- 4 * ((1:n) - (n - 1) / 2) / n / (n - 1)
  sum(w * sort(x - mean(x)))
}
