#' @title Weight a variable
#' @name weight
#' @description These functions weight the variable \code{x} by
#'                a specific vector of \code{weights}.
#'
#' @param x (Unweighted) variable.
#' @param weights Vector with same length as \code{x}, which
#'          contains weight factors. Each value of \code{x} has a
#'          specific assigned weight in \code{weights}.
#' @param digits Numeric value indicating the number of decimal places to be
#'          used for rounding the weighted values. By default, this value is
#'          \code{0}, i.e. the returned values are integer values.
#'
#' @return The weighted \code{x}.
#'
#' @details \code{weight2()} sums up all \code{weights} values of the associated
#'            categories of \code{x}, whereas \code{weight()} uses a
#'            \code{\link[stats]{xtabs}} formula to weight cases. Thus, \code{weight()}
#'            may return a vector of different length than \code{x}.
#'
#' @note The values of the returned vector are in sorted order, whereas the values'
#'        order of the original \code{x} may be spread randomly. Hence, \code{x} can't be
#'        used, for instance, for further cross tabulation. In case you want to have
#'        weighted contingency tables or (grouped) box plots etc., use the \code{weightBy}
#'        argument of most functions.
#'
#' @examples
#' v <- sample(1:4, 20, TRUE)
#' table(v)
#' w <- abs(rnorm(20))
#' table(weight(v, w))
#' table(weight2(v, w))
#'
#' set.seed(1)
#' x <- sample(letters[1:5], size = 20, replace = TRUE)
#' w <- runif(n = 20)
#'
#' table(x)
#' table(weight(x, w))
#'
#' @importFrom stats na.pass xtabs
#' @importFrom sjlabelled as_numeric
#' @export
weight <- function(x, weights, digits = 0) {
  # remember if x is numeric
  x.is.num <- is.numeric(x)

  # init values
  weightedvar <- c()

  wtab <- round(stats::xtabs(weights ~ x,
                             data = data.frame(weights = weights, x = x),
                             na.action = stats::na.pass,
                             exclude = NULL),
                digits = digits)

  # iterate all table values
  for (w in seq_len(length(wtab))) {
    # retrieve count of each table cell
    w_count <- wtab[[w]]

    # retrieve "cell name" which is identical to the variable value
    # first check whether values are numeric or not
    nval_ <- suppressWarnings(as.numeric(names(wtab[w])))

    # if value is not numeric, use as is
    if (is.na(nval_))
      w_value <- names(wtab[w])
    else
      # else, use numeric value
      w_value <- nval_

    # append variable value, repeating it "w_count" times.
    weightedvar <- c(weightedvar, rep(w_value, w_count))
  }

  # if we have NA values, weighted var is coerced to character.
  # coerce back to numeric then here
  if (!is.numeric(weightedvar) && x.is.num) weightedvar <- sjlabelled::as_numeric(weightedvar)

  # return result
  weightedvar
}

#' @rdname weight
#' @export
weight2 <- function(x, weights) {
  items <- unique(x)
  newvar <- c()

  for (i in seq_len(length(items))) {
    newcount <- round(sum(weights[which(x == items[i])]))
    newvar <- c(newvar, rep(items[i], newcount))
  }

  newvar
}
