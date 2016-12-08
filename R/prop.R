#' @title Proportion of values in a vector
#' @name prop
#'
#' @description This function calculates the proportion of a value or category
#'                in a variable.
#'
#' @param data A data frame.
#' @param ... One or more value pairs of comparisons (logical predicates). Put
#'              variable names the left-hand-side and values to match on the
#'              right hand side. Expressions may be quoted or unquoted. See
#'              'Examples'.
#' @param weight.by Vector of weights that will be applied to weight all observations.
#'          Must be a vector of same length as the input vector. Default is
#'          \code{NULL}, so no weights are used.
#' @param na.rm Logical, whether to remove NA values from the vector when the
#'          proportion is calculated. \code{na.rm = FALSE} gives you the raw
#'          percentage of a value in a vector, \code{na.rm = TRUE} the valid
#'          percentage.
#'
#' @return For one condition, a numeric value with the proportion of the values
#'         inside a vector. For more than one condition, a tibble with one column
#'         of conditions and one column with proportions.
#'
#' @examples
#' data(efc)
#'
#' # proportion of value 1 in e42dep
#' prop(efc, e42dep == 1)
#'
#' # expression may also be completely quotes
#' prop(efc, "e42dep == 1")
#'
#' # proportion of value 1 in e42dep, and all values greater
#' # than 2 in e42dep, excluding missing values. will return a tibble
#' prop(efc, e42dep == 1, e42dep > 2, na.rm = TRUE)
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
#' # also works with pipe-chains
#' library(dplyr)
#' efc %>% prop(e17age > 70)
#' efc %>% summarise(age70 = prop(., e17age > 70))
#'
#' # and with group_by
#' efc %>%
#'   group_by(e16sex) %>%
#'   summarise(hi.dependency = prop(., e42dep > 2))
#'
#' @importFrom tibble tibble
#' @export
prop <- function(data, ..., weight.by = NULL, na.rm = FALSE) {
  # check argument
  if (!is.data.frame(data)) stop("`data` needs to be a data frame.", call. = F)

  # get dots
  dots <- match.call(expand.dots = FALSE)$`...`

  # remember comparisons
  comparisons <- lapply(dots, function(x) {
    # to character, and remove spaces and quotes
    x <- gsub(" ", "", deparse(x), fixed = T)
    x <- gsub("\"", "", x, fixed = TRUE)
    x
  })

  # iterate dots
  result <- lapply(dots, function(x) {
    # check if we have a structured pairlist (e.g. from 'group_by()' of dplyr)
    if (startsWith(deparse(x)[1], "structure(")) {
      # in this case, we just need to evaluate the expression, because the
      # data values are already given as structure
      dummy <- eval(x)
    } else {
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

      # weight vector?
      if (!is.null(weight.by)) f <- weight(f, weights = weight.by)

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
    }

    # remove missings?
    if (na.rm) dummy <- na.omit(dummy)

    # get proportion
    sum(dummy, na.rm = T) / length(dummy)
  })

  # if we have more than one proportion, return a tibble. this allows us
  # to save more information, the condition and the proportion value
  if (length(comparisons) > 1) {
    return(tibble::tibble(condition = as.character(unlist(comparisons)), prop = unlist(result)))
  }

  unlist(result)
}