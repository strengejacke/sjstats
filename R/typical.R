#' @title Return the typical value of a vector
#' @name typical_value
#'
#' @description This function returns the "typical" value of a variable.
#'
#'
#' @param x A variable.
#' @param fun Character vector, naming the function to be applied to
#'        \code{x}. Currently, \code{"mean"}, \code{"weighted.mean"},
#'        \code{"median"} and \code{"mode"} are supported, which call the
#'        corresponding R functions (except \code{"mode"}, which calls an
#'        internal function to compute the most common value). \code{"zero"}
#'        simply returns 0. \strong{Note:} By default, if \code{x} is a factor,
#'        only \code{fun = "mode"} is applicable; for all other functions (including
#'        the default, \code{"mean"}) the reference level of \code{x} is returned.
#'        For character vectors, only the mode is returned. You can use a named
#'        vector to apply other different functions to numeric and categorical
#'        \code{x}, where factors are first converted to numeric vectors, e.g.
#'        \code{fun = c(numeric = "median", factor = "mean")}. See 'Examples'.
#' @param ... Further arguments, passed down to \code{fun}.
#'
#' @inheritParams grpmean
#'
#' @return The "typical" value of \code{x}.
#'
#' @details By default, for numeric variables, \code{typical_value()} returns the
#'          mean value of \code{x} (unless changed with the \code{fun}-argument).
#'          \cr \cr
#'          For factors, the reference level is returned or the most common value
#'          (if \code{fun = "mode"}), unless \code{fun} is a named vector. If
#'          \code{fun} is a named vector, specify the function for numeric
#'          and categorical variables as element names, e.g.
#'          \code{fun = c(numeric = "median", factor = "mean")}. In this case,
#'          factors are converted to numeric values (using \code{\link[sjmisc]{to_value}})
#'          and the related function is applied. You may abbreviate the names
#'          \code{fun = c(n = "median", f = "mean")}. See also 'Examples'.
#'          \cr \cr
#'          For character vectors the most common value (mode) is returned.
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
#' typical_value(x, "weighted.mean", weight.by = wt)
#'
#' # for factors, return either reference level or mode value
#' set.seed(123)
#' x <- sample(iris$Species, size = 30, replace = TRUE)
#' typical_value(x)
#' typical_value(x, fun = "mode")
#'
#' # for factors, use a named vector to apply other functions than "mode"
#' map(iris, ~ typical_value(.x, fun = c(n = "median", f = "mean")))
#'
#'
#' @export
typical_value <- function(x, fun = "mean", weight.by = NULL, ...) {

  # check if we have named vectors and find the requested function
  # for special functions for factors, convert to numeric first

  fnames <- names(fun)

  if (!is.null(fnames)) {
    if (is.numeric(x)) {
      fun <- fun[which(fnames %in% c("numeric", "n"))]
    } else if (is.factor(x)) {
      fun <- fun[which(fnames %in% c("factor", "f"))]
      x <- sjmisc::to_value(x, keep.labels = FALSE)
    }
  }


  if (!fun %in% c("mean", "median", "mode", "weighted.mean", "zero"))
    stop("`fun` must be one of \"mean\", \"median\", \"mode\", \"weighted.mean\" or \"zero\".", call. = F)


  # for weighted mean, check that weights are of same length as x

  if (fun == "weighted.mean" && !is.null(weight.by)) {

    # make sure weights and x have same length

    if (length(weight.by) != length(x)) {
      # if not, tell user and change function to mean
      warning("Vector of weights is of different length than `x`. Using `mean` as function for typical value.", call. = F)
      fun <- "mean"
    }


    # make sure weights are differen from 1

    if (all(weight.by == 1)) {
      # if not, tell user and change function to mean
      warning("All weight values are `1`. Using `mean` as function for typical value.", call. = F)
      fun <- "mean"
    }
  }


  # no weights, than use normal mean function

  if (fun == "weighted.mean" && is.null(weight.by)) fun <- "mean"


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

  if (is.numeric(x)) {
    if (fun == "weighted.mean")
      do.call(myfun, args = list(x = x, na.rm = TRUE, w = weight.by, ...))
    else
      do.call(myfun, args = list(x = x, na.rm = TRUE, ...))
  } else if (is.factor(x)) {
    if (fun != "mode")
      levels(x)[1]
    else
      mode_value(x)
  } else {
    mode_value(x)
  }
}


mode_value <- function(x, ...) {
  # create frequency table, to find most common value
  counts <- table(x)
  modus <- names(counts)[max(counts) == counts]

  # in case all values appear equally often, use first value
  if (length(modus) > 1) modus <- modus[1]

  # check if it's numeric
  if (!is.na(suppressWarnings(as.numeric(modus))))
    as.numeric(modus)
  else
    modus
}
