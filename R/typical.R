#' @title Return the typical value of a vector
#' @name typical_value
#'
#' @description This function returns the "typical" value of a variable.
#'
#'
#' @param x A variable.
#' @param fun Character vector, naming the function to be applied to numeric
#'        \code{x}. Currently, \code{"mean"}, \code{"weighted.mean"},
#'        \code{"median"} and \code{"mode"} are supported, which call the
#'        corresponding R functions (except \code{"mode"}, which calls an
#'        internal function to compute the most common value). \code{"zero"}
#'        simply returns 0.
#' @param ... Further arguments, passed down to \code{fun}.
#'
#' @return The "typical" value of \code{x}.
#'
#' @details By default, for numeric variables, \code{typical_value()} returns the
#'          mean value of \code{x} (unless changed with the \code{fun}-argument).
#'          For factors, the reference level and for character vectors, to the
#'          most common value (mode) is returned.

#'
#' @examples
#' data(iris)
#' typical_value(iris$Sepal.Length)
#'
#' library(purrr)
#' map(iris, ~ typical_value(.x))
#'
#' # example from ?stats::weighted.mean
#' wt <- c(5,  5,  4,  1) / 15
#' x <- c(3.7, 3.3, 3.5, 2.8)
#'
#' typical_value(x, "weighted.mean")
#' typical_value(x, "weighted.mean", w = wt)
#'
#' @export
typical_value <- function(x, fun = c("mean", "median", "mode", "weighted.mean", "zero"), ...) {
  fun <- match.arg(fun)

  if (fun == "median")
    myfun <- get("median", asNamespace("stats"))
  else if (fun == "weighted.mean")
    myfun <- get("weighted.mean", asNamespace("stats"))
  else if (fun == "mode")
    myfun <- get("mode_value", asNamespace("sjstats"))
  else if (fun == "zero")
    return(0)
  else
    myfun <- get("mean", asNamespace("base"))

  if (is.numeric(x))
    do.call(myfun, args = list(x = x, na.rm = TRUE, ...))
  else if (is.factor(x))
    levels(x)[1]
  else
    mode_value(x)
}


mode_value <- function(x, ...) {
  # create frequency table, to find most common value
  counts <- table(x)
  modus <- names(counts)[max(counts) == counts]

  # check if it's numeric
  if (!is.na(as.numeric(modus)))
    as.numeric(modus)
  else
    modus
}
