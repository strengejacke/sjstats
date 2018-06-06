#' @title Row means with min amount of valid values
#' @name mean_n
#' @description This function is similar to the SPSS \code{MEAN.n} function and computes
#'                row means from a \code{data.frame} or \code{matrix} if at least \code{n}
#'                values of a row are valid (and not \code{NA}).
#'
#' @param dat A data frame with at least two columns, where row means are applied.
#' @param n May either be
#'          \itemize{
#'            \item a numeric value that indicates the amount of valid values per row to calculate the row mean;
#'            \item or a value between 0 and 1, indicating a proportion of valid values per row to calculate the row mean (see 'Details').
#'          }
#'          If a row's sum of valid values is less than \code{n}, \code{NA} will be returned as row mean value.
#' @param digits Numeric value indicating the number of decimal places to be used for rounding mean
#'          value. Negative values are allowed (see ‘Details’).
#'
#' @return A vector with row mean values of \code{df} for those rows with at least \code{n}
#'           valid values. Else, \code{NA} is returned.
#'
#' @details Rounding to a negative number of \code{digits} means rounding to a power of
#'            ten, so for example mean_n(df, 3, digits = -2) rounds to the
#'            nearest hundred. \cr \cr
#'          For \code{n}, must be a numeric value from \code{0} to \code{ncol(dat)}. If
#'            a \emph{row} in \code{dat} has at least \code{n} non-missing values, the
#'            row mean is returned. If \code{n} is a non-integer value from 0 to 1,
#'            \code{n} is considered to indicate the proportion of necessary non-missing
#'            values per row. E.g., if \code{n = .75}, a row must have at least \code{ncol(dat) * n}
#'            non-missing values for the row mean to be calculated. See 'Examples'.
#'
#' @references \href{http://r4stats.com/2014/09/03/adding-the-spss-mean-n-function-to-r/}{r4stats.com}
#'
#' @examples
#' dat <- data.frame(c1 = c(1,2,NA,4),
#'                   c2 = c(NA,2,NA,5),
#'                   c3 = c(NA,4,NA,NA),
#'                   c4 = c(2,3,7,8))
#'
#' # needs at least 4 non-missing values per row
#' mean_n(dat, 4) # 1 valid return value
#'
#' # needs at least 3 non-missing values per row
#' mean_n(dat, 3) # 2 valid return values
#'
#' # needs at least 2 non-missing values per row
#' mean_n(dat, 2)
#'
#' # needs at least 1 non-missing value per row
#' mean_n(dat, 1) # all means are shown
#'
#' # needs at least 50% of non-missing values per row
#' mean_n(dat, .5) # 3 valid return values
#'
#' # needs at least 75% of non-missing values per row
#' mean_n(dat, .75) # 2 valid return values
#'
#' @export
mean_n <- function(dat, n, digits = 2) {
  # is 'n' indicating a proportion?
  digs <- n %% 1
  if (digs != 0) n <- round(ncol(dat) * digs)

  # coerce matrix to data frame
  if (is.matrix(dat)) dat <- as.data.frame(dat)

  # check if we have a data framme with at least two columns
  if (!is.data.frame(dat) || ncol(dat) < 2) {
    warning("`dat` must be a data frame with at least two columns.", call. = TRUE)
    return(NA)
  }

  # n may not be larger as df's amount of columns
  if (ncol(dat) < n) {
    warning("`n` must be smaller or equal to number of columns in data frame.", call. = TRUE)
    return(NA)
  }

  round(apply(dat, 1, function(x) ifelse(sum(!is.na(x)) >= n, mean(x, na.rm = TRUE), NA)), digits)
}
