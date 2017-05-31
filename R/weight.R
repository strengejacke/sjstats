#' @title Weight a variable
#' @name weight
#' @description These functions weight the variable \code{x} by
#'                a specific vector of \code{weights}.
#'
#' @param x (Unweighted) variable
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
#'            \code{\link{xtabs}} formula to weight cases. Thus, \code{weight()}
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
  return(weightedvar)
}

#' @rdname weight
#' @export
weight2 <- function(x, weights) {
  items <- unique(x)
  newvar <- c()
  for (i in seq_len(length(items))) {
    newcount = round(sum(weights[which(x == items[i])]))
    newvar <- c(newvar, rep(items[i], newcount))
  }
  return(newvar)
}


#' @title Weighted statistics for variables
#' @name wtd_sd
#' @description Compute weighted standard deviation or standard error for a
#'                variable or for all variables of a data frame.
#'
#' @param x (Numeric) vector or a data frame.
#' @param weights Numeric vector of weights.
#' @return The weighted standard deviation or standard error of \code{x},
#'           or for each variable if \code{x} is a data frame.
#'
#' @examples
#' wtd_sd(rnorm(n = 100, mean = 3),
#'        runif(n = 100))
#'
#' data(efc)
#' wtd_sd(efc[, 1:3], runif(n = nrow(efc)))
#' wtd_se(efc[, 1:3], runif(n = nrow(efc)))
#'
#' @export
wtd_sd <- function(x, weights = NULL) {
  # check if suggested packages are available
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop("Package `Hmisc` needed for this function to work. Please install it.", call. = FALSE)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    # init return variables
    stdd <- c()
    stdd_names <- c()
    # iterate all columns
    for (i in seq_len(ncol(x))) {
      # get and save standard error for each variable
      # of the data frame
      stdd <- c(stdd,
                sqrt(Hmisc::wtd.var(x[[i]], weights = weights, na.rm = TRUE)))
      # save column name as variable name
      stdd_names <- c(stdd_names, colnames(x)[i])
    }
    # set names to return vector
    names(stdd) <- stdd_names
    # return results
    return(stdd)
  } else {
    return(sqrt(Hmisc::wtd.var(x, weights = weights, na.rm = TRUE)))
  }
}


#' @rdname wtd_sd
#' @export
wtd_se <- function(x, weights = NULL) {
  # check if suggested packages are available
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop("Package `Hmisc` needed for this function to work. Please install it.", call. = FALSE)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    # init return variables
    stde <- c()
    stde_names <- c()
    # iterate all columns
    for (i in seq_len(ncol(x))) {
      # get and save standard error for each variable
      # of the data frame
      stde <- c(stde,
                sqrt(Hmisc::wtd.var(x[[i]], weights = weights, na.rm = TRUE) / length(stats::na.omit(x[[i]]))))
      # save column name as variable name
      stde_names <- c(stde_names, colnames(x)[i])
    }
    # set names to return vector
    names(stde) <- stde_names
    # return results
    return(stde)
  } else {
    return(sqrt(Hmisc::wtd.var(x, weights = weights, na.rm = TRUE) / length(stats::na.omit(x))))
  }
}
