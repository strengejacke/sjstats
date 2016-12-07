#' @title Standardize and center variables
#' @name std
#' @description \code{std()} computes a z-transformation (standardized and centered)
#'              on the input. \code{center()} centers the input.
#'
#' @param data A variable that should be standardized or centered, or a
#'          data frame with such variables.
#' @param ... Optional, names of the variables with that should be standardized
#'          or centered. Required, if either \code{data} is a data frame and no vector,
#'          and only selected variables from \code{data} should be used
#'          in the function.
#' @param include.fac Logical, if \code{TRUE}, factors will be converted to numeric
#'          vectors and also standardized or centered.
#' @param append Logical, if \code{TRUE} and \code{data} is a data frame, the
#'          standardized or centered variables will be appended as new columns
#'          to \code{data}. Variable (column) names will be the same, however,
#'          standardized variables get the suffix \code{"_z"}, while centered
#'          variables get the suffix \code{"_c"}.
#'
#' @return Either a vector with standardized or centered variables, if \code{data}
#'         was a vector; or a \code{\link[tibble]{tibble}} with standardized or
#'         centered variables, if \code{data} was a data frame. If \code{append = TRUE},
#'         standardized and centered variables are added as new columns to \code{data}.
#'
#' @note \code{std()} and \code{center()} only return a vector, if \code{data} is
#'         a vector. If \code{data} is a data frame and only one variable is specified
#'         in the \code{...}-ellipses argument, both functions do return a
#'         data frame (see 'Examples').
#'
#' @examples
#' library(dplyr)
#' data(efc)
#' std(efc$c160age) %>% head()
#' std(efc, e17age, c160age) %>% head()
#' std(efc, e17age, c160age, append = TRUE) %>% head()
#'
#' center(efc$c160age) %>% head()
#' center(efc, e17age, c160age) %>% head()
#' center(efc, e17age, c160age, append = TRUE) %>% head()
#'
#' # NOTE!
#' std(efc$e17age) # returns a vector
#' std(efc, e17age) # returns a tibble
#'
#' @export
std <- function(data, ..., include.fac = TRUE, append = FALSE) {
  if (is.data.frame(data)) {
    .dat <- get_boot_data(data, match.call(expand.dots = FALSE)$`...`)
    tmp <- tibble::as_tibble(lapply(.dat, FUN = std_helper,
                                    include.fac = include.fac,
                                    standardize = TRUE))
    # check whether standardized variables should be appended to data frame
    if (append) {
      # change variable names, add suffix "_z"
      colnames(tmp) <- sprintf("%s_z", colnames(tmp))
      # append data
      tmp <- dplyr::bind_cols(data, tmp)
    }
    return(tmp)
  } else {
    return(std_helper(data, include.fac = include.fac, standardize = TRUE))
  }
}


#' @rdname std
#' @export
center <- function(data, ..., include.fac = TRUE, append = FALSE) {
  if (is.data.frame(data)) {
    .dat <- get_boot_data(data, match.call(expand.dots = FALSE)$`...`)
    tmp <- tibble::as_tibble(lapply(.dat, FUN = std_helper,
                                    include.fac = include.fac,
                                    standardize = FALSE))
    # check whether standardized variables should be appended to data frame
    if (append) {
      # change variable names, add suffix "_z"
      colnames(tmp) <- sprintf("%s_c", colnames(tmp))
      # append data
      tmp <- dplyr::bind_cols(data, tmp)
    }
    return(tmp)
  } else {
    return(std_helper(data, include.fac = include.fac, standardize = FALSE))
  }
}


#' @importFrom stats sd na.omit
#' @importFrom sjmisc to_value get_label set_label
std_helper <- function(data, include.fac, standardize) {
  # check whether factors should also be standardized
  if (is.factor(data)) {
    if (include.fac)
      data <- sjmisc::to_value(data, keep.labels = FALSE)
    else
      return(data)
  }
  # remove missings
  tmp <- stats::na.omit(data)
  # save value label, if any
  lab <- sjmisc::get_label(data)
  # now center and standardize
  tmp <- tmp - mean(tmp)
  if (standardize) tmp <- tmp / stats::sd(tmp)
  # and fill in values in original vector
  data[!is.na(data)] <- tmp
  # add back label
  sjmisc::set_label(data, lab)
}