#' @title Calculate population variance and standard deviation
#' @name var_pop
#' @description Calculate the population variance or standard deviation of a vector.
#'
#' @param x (Numeric) vector.
#'
#' @return The population variance or standard deviation of \code{x}.
#'
#' @details Unlike \code{\link[stats]{var}}, which returns the sample variance,
#'            \code{var_pop()} returns the population variance. \code{sd_pop()}
#'            returns the standard deviation based on the population variance.
#'
#' @examples
#' data(efc)
#'
#' # sampling variance
#' var(efc$c12hour, na.rm = TRUE)
#' # population variance
#' var_pop(efc$c12hour)
#'
#' # sampling sd
#' sd(efc$c12hour, na.rm = TRUE)
#' # population sd
#' sd_pop(efc$c12hour)
#'
#' @importFrom stats na.omit var
#' @importFrom sjmisc is_num_fac to_value
#' @export
var_pop <- function(x) {
  # check for categorical
  if (is.factor(x)) {
    # only allow numeric factors
    if (!sjmisc::is_num_fac(x)) {
      warning("`x` must be numeric vector or a factor with numeric levels.", call. = F)
      return(NA)
    }
    # convert factor to numeric
    x <- sjmisc::to_value(x)
  }
  # remove NA
  x <- stats::na.omit(x)
  n <- length(x)
  # population variance
  stats::var(x) * ((n - 1) / n)
}


#' @rdname var_pop
#' @importFrom stats na.omit var
#' @export
sd_pop <- function(x) {
  # get population variance
  pv <- var_pop(x)
  # factors with non-numeric level return NULL
  if (!is.null(pv) && !is.na(pv))
    sqrt(pv)
  else
    NA
}
