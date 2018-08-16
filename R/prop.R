#' @title Proportions of values in a vector
#' @name prop
#'
#' @description \code{prop()} calculates the proportion of a value or category
#'                in a variable. \code{props()} does the same, but allows for
#'                multiple logical conditions in one statement. It is similar
#'                to \code{mean()} with logical predicates, however, both
#'                \code{prop()} and \code{props()} work with grouped data frames.
#'
#' @param data A data frame. May also be a grouped data frame (see 'Examples').
#' @param ... One or more value pairs of comparisons (logical predicates). Put
#'              variable names the left-hand-side and values to match on the
#'              right hand side. Expressions may be quoted or unquoted. See
#'              'Examples'.
#' @param weights Vector of weights that will be applied to weight all observations.
#'          Must be a vector of same length as the input vector. Default is
#'          \code{NULL}, so no weights are used.
#' @param na.rm Logical, whether to remove NA values from the vector when the
#'          proportion is calculated. \code{na.rm = FALSE} gives you the raw
#'          percentage of a value in a vector, \code{na.rm = TRUE} the valid
#'          percentage.
#' @inheritParams reliab_test
#'
#' @details \code{prop()} only allows one logical statement per comparison,
#'          while \code{props()} allows multiple logical statements per comparison.
#'          However, \code{prop()} supports weighting of variables before calculating
#'          proportions, and comparisons may also be quoted. Hence, \code{prop()}
#'          also processes comparisons, which are passed as character vector
#'          (see 'Examples').
#'
#'
#' @return For one condition, a numeric value with the proportion of the values
#'         inside a vector. For more than one condition, a tibble with one column
#'         of conditions and one column with proportions. For grouped data frames,
#'         returns a tibble with one column per group with grouping categories,
#'         followed by one column with proportions per condition.
#'
#' @examples
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
#' # than 2 in e42dep, including missing values. will return a tibble
#' prop(efc, e42dep == 1, e42dep > 2, na.rm = FALSE)
#'
#' # for factors or character vectors, use quoted or unquoted values
#' library(sjmisc)
#' # convert numeric to factor, using labels as factor levels
#' efc$e16sex <- to_label(efc$e16sex)
#' efc$n4pstu <- to_label(efc$n4pstu)
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
#'
#' # also works with pipe-chains
#' library(dplyr)
#' efc %>% prop(e17age > 70)
#' efc %>% prop(e17age > 70, e16sex == 1)
#'
#' # and with group_by
#' efc %>%
#'   group_by(e16sex) %>%
#'   prop(e42dep > 2)
#'
#' efc %>%
#'   select(e42dep, c161sex, c172code, e16sex) %>%
#'   group_by(c161sex, c172code) %>%
#'   prop(e42dep > 2, e16sex == 1)
#'
#' # same for "props()"
#' efc %>%
#'   select(e42dep, c161sex, c172code, c12hour, n4pstu) %>%
#'   group_by(c161sex, c172code) %>%
#'   props(
#'     e42dep > 2,
#'     c12hour > 20 & c12hour < 40,
#'     n4pstu == 'Care Level 1' | n4pstu == 'Care Level 3'
#'   )
#'
#' @importFrom dplyr bind_cols bind_rows
#' @importFrom sjlabelled get_label get_labels as_numeric
#' @export
prop <- function(data, ..., weights = NULL, na.rm = TRUE, digits = 4) {
  # check argument
  if (!is.data.frame(data)) stop("`data` needs to be a data frame.", call. = F)

  # get dots
  dots <- match.call(expand.dots = FALSE)$`...`

  proportions(data, dots, weight.by = weights, na.rm, digits, multi_logical = FALSE)
}


#' @rdname prop
#' @export
props <- function(data, ..., na.rm = TRUE, digits = 4) {
  # check argument
  if (!is.data.frame(data)) stop("`data` needs to be a data frame.", call. = F)

  # get dots
  dots <- match.call(expand.dots = FALSE)$`...`

  proportions(data, dots, NULL, na.rm, digits, multi_logical = TRUE)
}


#' @importFrom purrr map_df
proportions <- function(data, dots, weight.by, na.rm, digits, multi_logical) {
  # remember comparisons
  comparisons <- lapply(dots, function(x) {
    # to character, and remove spaces and quotes
    x <- gsub(" ", "", deparse(x), fixed = T)
    x <- gsub("\"", "", x, fixed = TRUE)
    x
  })

  # do we have a grouped data frame?
  if (inherits(data, "grouped_df")) {

    # remember order of values
    reihenfolge <- NULL

    # get grouped data
    grps <- get_grouped_data(data)

    # now get proportions for each subset
    fr <- purrr::map_df(
      seq_len(nrow(grps)),
      function(i) {
        # get data from grouped data frame
        .d <- grps$data[[i]]

        # iterate dots (comparing conditions)
        if (multi_logical)
          result <- lapply(dots, get_multiple_proportion, .d, na.rm, digits)
        else
          result <- lapply(dots, get_proportion, .d, weight.by, na.rm, digits)

        as.data.frame(t(unlist(result)))
      }
    )


    # now we need the values from the groups of the grouped data frame
    for (i in (ncol(grps) - 1):1) {
      # get value label
      var.name <- colnames(grps)[i]
      val.labels <- suppressWarnings(
        rep(sjlabelled::get_labels(data[[var.name]]), length.out = nrow(fr))
      )

      # if we have no value labels, use values instead
      if (is.null(val.labels)) {
        val.labels <-
          rep(unique(sort(data[[var.name]])), length.out = nrow(fr))
      }

      # add row order, based on values of grouping variables
      reihenfolge <- rep(sort(unique(sjlabelled::as_numeric(data[[var.name]]))), length.out = nrow(fr)) %>%
        as.data.frame() %>%
        dplyr::bind_cols(reihenfolge)

      # bind values as column
      fr <- dplyr::bind_cols(data.frame(val.labels, stringsAsFactors = FALSE), fr)
    }

    # get column names. we need variable labels as column names
    var.names <- colnames(grps)[seq_len(ncol(grps) - 1)]
    var.labels <- sjlabelled::get_label(data[, var.names], def.value = var.names)

    # set variable labels and comparisons as colum names
    colnames(fr) <- c(var.labels, comparisons)

    # order rows by values of grouping variables
    fr <- fr[do.call(order, reihenfolge), ]

    return(fr)

  } else {
    # iterate dots (comparing conditions)
    if (multi_logical)
      result <- lapply(dots, get_multiple_proportion, data, na.rm, digits)
    else
      result <- lapply(dots, get_proportion, data, weight.by, na.rm, digits)

    # if we have more than one proportion, return a tibble. this allows us
    # to save more information, the condition and the proportion value
    if (length(comparisons) > 1) {
      return(data_frame(
        condition = as.character(unlist(comparisons)),
        prop = unlist(result)
      ))
    }

    return(unlist(result))
  }
}


get_proportion <- function(x, data, weight.by, na.rm, digits) {
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

  # remove missings?
  if (na.rm) dummy <- na.omit(dummy)

  # get proportion
  round(sum(dummy, na.rm = T) / length(dummy), digits = digits)
}


get_multiple_proportion <- function(x, data, na.rm, digits) {
  # evaluate argument
  dummy <- with(data, eval(parse(text = deparse(x))))

  # remove missings?
  if (na.rm) dummy <- na.omit(dummy)

  # get proportion
  round(sum(dummy, na.rm = T) / length(dummy), digits = digits)
}
