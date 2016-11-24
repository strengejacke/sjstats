#' @title Proportion of values in a vector
#' @name prop
#'
#' @description This function calculates the proportion of a value or category
#'                in a variable.
#'
#' @param data A data frame.
#' @param ... One or more value pairs. Put variable names the left-hand-side and
#' @param na.rm Logical, whether to remove NA values from the vector when the
#'          proportion is calculated. \code{na.rm = FALSE} gives you the raw
#'          percentage of a value in a vector, \code{na.rm = TRUE} the valid
#'          percentage.
#'
#' @return Numeric, the proportion of the values inside a vector.
#'
#' @examples
#' data(efc)
#'
#' # proportion of value 1 in e42dep
#' prop(efc, e42dep == 1)
#'
#' # proportion of value 1 in e42dep, and all values greater
#' # than 2 in e42dep, excluding missing values
#' prop(efc, e42dep == 1, e42dep > 2, na.rm = T)
#'
#'
#' # for factors or character vectors, use quoted or unquoted values
#' library(sjmisc)
#' # convert numeric to factor, using labels as factor levels
#' efc$e16sex <- to_label(efc$e16sex)
#'
#' # get proportion of female older persons
#' prop(efc, e16sex == female)
#'
#' # get proportion of male older persons
#' prop(efc, e16sex == "male")
#'
#' @export
prop <- function(data, ..., na.rm = FALSE) {
  # check argument
  if (!is.data.frame(data)) stop("`data` needs to be a data frame.", call. = F)

  # get dots
  dots <- match.call(expand.dots = FALSE)$`...`

  # iterate dots
  result <- lapply(dots, function(x) {
    # to character, and remove spaces and quotes
    x <- gsub(" ", "", deparse(x), fixed = T)
    x <- gsub("\"", "", x, fixed = TRUE)
    # split expression at ==, < or >
    x.parts <- unlist(regmatches(x, gregexpr("[!=]=|[<>]|(?:(?![=!]=)[^<>])+", x, perl = TRUE)))
    # shorter version, however, does not split variable names with dots
    # x.parts <- unlist(regmatches(x, regexec("(\\w+)(\\W+)(\\w+)", x)))[-1]

    # correct == assignment?
    if (length(x.parts) < 3) {
      message("?Syntax error in argument. You possibly used `=` instead of `==`.")
      return(NULL)
    }

    # get variable from data and value from equation
    f <- data[[x.parts[1]]]
    v <- suppressWarnings(as.numeric(x.parts[3]))
    # if we have factor, values maybe non-numeric
    if (is.na(v)) v <- x.parts[3]

    # get proportions
    if (x.parts[2] == "==")
      dummy <- f == v
    else if (x.parts[2] == "!=")
      dummy <- f != v
    else if (x.parts[2] == "<")
      dummy <- f < v
    else if (x.parts[2] == ">")
      dummy <- f > v
    else
      dummy <- f == v

    # remove missings?
    if (na.rm) dummy <- na.omit(dummy)

    # get proportion
    res <- sum(dummy, na.rm = T) / length(dummy)
  })

  unlist(result)
}