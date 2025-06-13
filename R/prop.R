#' @title Proportions of values in a vector
#' @name prop
#'
#' @description `prop()` calculates the proportion of a value or category
#'                in a variable. `props()` does the same, but allows for
#'                multiple logical conditions in one statement. It is similar
#'                to `mean()` with logical predicates, however, both
#'                `prop()` and `props()` work with grouped data frames.
#'
#' @param data A data frame. May also be a grouped data frame (see 'Examples').
#' @param ... One or more value pairs of comparisons (logical predicates). Put
#'              variable names the left-hand-side and values to match on the
#'              right hand side. Expressions may be quoted or unquoted. See
#'              'Examples'.
#' @param weights Vector of weights that will be applied to weight all observations.
#'          Must be a vector of same length as the input vector. Default is
#'          `NULL`, so no weights are used.
#' @param na.rm Logical, whether to remove NA values from the vector when the
#'          proportion is calculated. `na.rm = FALSE` gives you the raw
#'          percentage of a value in a vector, `na.rm = TRUE` the valid
#'          percentage.
#' @param digits Amount of digits for returned values.
#'
#' @details `prop()` only allows one logical statement per comparison,
#'          while `props()` allows multiple logical statements per comparison.
#'          However, `prop()` supports weighting of variables before calculating
#'          proportions, and comparisons may also be quoted. Hence, `prop()`
#'          also processes comparisons, which are passed as character vector
#'          (see 'Examples').
#'
#'
#' @return For one condition, a numeric value with the proportion of the values
#'         inside a vector. For more than one condition, a data frame with one column
#'         of conditions and one column with proportions. For grouped data frames,
#'         returns a data frame with one column per group with grouping categories,
#'         followed by one column with proportions per condition.
#'
#' @examplesIf getRversion() >= "4.2.0" && requireNamespace("datawizard", quietly = TRUE)
#' data(efc)
#'
#' # proportion of value 1 in e42dep
#' prop(efc, e42dep == 1)
#'
#' # expression may also be completely quoted
#' prop(efc, "e42dep == 1")
#'
#' # use "props()" for multiple logical statements
#' props(efc, e17age > 70 & e17age < 80)
#'
#' # proportion of value 1 in e42dep, and all values greater
#' # than 2 in e42dep, including missing values. will return a data frame
#' prop(efc, e42dep == 1, e42dep > 2, na.rm = FALSE)
#'
#' # for factors or character vectors, use quoted or unquoted values
#' library(datawizard)
#' # convert numeric to factor, using labels as factor levels
#' efc$e16sex <- to_factor(efc$e16sex)
#' efc$n4pstu <- to_factor(efc$n4pstu)
#'
#' # get proportion of female older persons
#' prop(efc, e16sex == female)
#'
#' # get proportion of male older persons
#' prop(efc, e16sex == "male")
#'
#' # "props()" needs quotes around non-numeric factor levels
#' props(efc,
#'   e17age > 70 & e17age < 80,
#'   n4pstu == 'Care Level 1' | n4pstu == 'Care Level 3'
#' )
#' @export
prop <- function(data, ..., weights = NULL, na.rm = TRUE, digits = 4) {
  # check argument
  if (!is.data.frame(data)) {
    insight::format_error("`data` needs to be a data frame.")
  }
  dots <- match.call(expand.dots = FALSE)[["..."]]
  .proportions(data, dots = dots, weight.by = weights, na.rm, digits, multi_logical = FALSE)
}


#' @rdname prop
#' @export
props <- function(data, ..., na.rm = TRUE, digits = 4) {
  # check argument
  if (!is.data.frame(data)) {
    insight::format_error("`data` needs to be a data frame.")
  }
  dots <- match.call(expand.dots = FALSE)[["..."]]
  .proportions(data, dots = dots, NULL, na.rm, digits, multi_logical = TRUE)
}


.proportions <- function(data, dots, weight.by, na.rm, digits, multi_logical) {
  # remember comparisons
  comparisons <- lapply(dots, function(x) {
    # to character, and remove spaces and quotes
    x <- gsub(" ", "", deparse(x), fixed = TRUE)
    x <- gsub("\"", "", x, fixed = TRUE)
    x
  })

  if (inherits(data, "grouped_df")) {
    grps <- attributes(data)$groups
    result <- lapply(grps[[".rows"]], function(x) {
      .process_prop(data[x, , drop = FALSE], comparisons, dots, multi_logical, na.rm, digits, weight.by)
    })
  } else {
    result <- .process_prop(data, comparisons, dots, multi_logical, na.rm, digits, weight.by)
  }
  result
}


.process_prop <- function(data, comparisons, dots, multi_logical, na.rm, digits, weight.by) {
  # iterate dots (comparing conditions)
  if (multi_logical)
    result <- lapply(dots, get_multiple_proportion, data, na.rm, digits)
  else
    result <- lapply(dots, get_proportion, data, weight.by, na.rm, digits)

  # if we have more than one proportion, return a data frame. this allows us
  # to save more information, the condition and the proportion value
  if (length(comparisons) > 1) {
    return(data_frame(
      condition = as.character(unlist(comparisons)),
      prop = unlist(result)
    ))
  }

  unlist(result)
}


get_proportion <- function(x, data, weight.by, na.rm, digits) {
  # to character, and remove spaces and quotes
  x <- gsub(" ", "", deparse(x), fixed = TRUE)
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
  dummy <- switch(x.parts[2],
    "==" = f == v,
    "!=" = f != v,
    "<" = f < v,
    ">" = f > v,
    f == v
  )

  # remove missings?
  if (na.rm) dummy <- stats::na.omit(dummy)

  # get proportion
  round(sum(dummy, na.rm = TRUE) / length(dummy), digits = digits)
}


get_multiple_proportion <- function(x, data, na.rm, digits) {
  # evaluate argument
  dummy <- with(data, eval(parse(text = deparse(x))))

  # remove missings?
  if (na.rm) dummy <- stats::na.omit(dummy)

  # get proportion
  round(sum(dummy, na.rm = TRUE) / length(dummy), digits = digits)
}
