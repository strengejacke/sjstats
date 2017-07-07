#' @title Return the typical value of a vector
#' @name typical_value
#'
#' @description This function returns the "typical" value of a variable.
#'
#'
#' @param x A variable.
#' @param fun Character vector, naming the function to be applied to numeric
#'        \code{x}. Currently, only \code{"mean"} and \code{"median"} are
#'        supported.
#'
#' @return The "typical" value of \code{x}.
#'
#' @details By default, for numeric variables, \code{typical_value()} returns the
#'          mean value of \code{x} (unless changed with the \code{fun}-argument).
#'          For factors, the reference level and for character vectors, to the
#'          most common value is returned.

#'
#' @examples
#' data(iris)
#' typical_value(iris$Sepal.Length)
#'
#' library(purrr)
#' map(iris, ~ typical_value(.x))
#'
#' @export
typical_value <- function(x, fun = c("mean", "median")) {
  fun <- match.arg(fun)

  if (fun == "median")
    myfun <- get("median", asNamespace("stats"))
  else
    myfun <- get("mean", asNamespace("base"))

  if (is.numeric(x))
    do.call(myfun, args = list(x = x, na.rm = TRUE))
  else if (is.factor(x))
    levels(x)[1]
  else {
    counts <- table(x)
    names(counts)[max(counts) == counts]
  }
}
