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


#' @title Weighted statistics for variables
#' @name wtd_sd
#' @description \code{wtd_sd()} and \code{wtd_se()} compute weighted standard
#'   deviation or standard error for a variable or for all variables of a data
#'   frame. \code{svy_md()} computes the median for a variable in a survey-design
#'   (see \code{\link[survey]{svydesign}}).
#'
#' @param x (Numeric) vector or a data frame. For \code{svy_md()}, the bare
#'          (unquoted) variable name, or a character vector with the variable
#'          name.
#' @param weights Numeric vector of weights.
#'
#' @inheritParams svyglm.nb
#'
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
#'
#' # median for variables from weighted survey designs
#' library(survey)
#' data(nhanes_sample)
#'
#' des <- svydesign(
#'   id = ~SDMVPSU,
#'   strat = ~SDMVSTRA,
#'   weights = ~WTINT2YR,
#'   nest = TRUE,
#'   data = nhanes_sample
#' )
#'
#' svy_md(total, des)
#' svy_md("total", des)
#'
#' @export
wtd_sd <- function(x, weights = NULL) {
  # check if suggested packages are available
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop("Package `Hmisc` needed for this function to work. Please install it.", call. = FALSE)
  }

  UseMethod("wtd_sd")
}


#' @export
wtd_sd.data.frame <- function(x, weights = NULL) {
  sd_result <- purrr::map_dbl(x, ~ sqrt(Hmisc::wtd.var(.x, weights = weights, na.rm = TRUE)))
  names(sd_result) <- colnames(x)

  sd_result
}

#' @export
wtd_sd.matrix <- function(x, weights = NULL) {
  sd_result <- purrr::map_dbl(x, ~ sqrt(Hmisc::wtd.var(.x, weights = weights, na.rm = TRUE)))
  names(sd_result) <- colnames(x)

  sd_result
}

#' @export
wtd_sd.default <- function(x, weights = NULL) {
  sqrt(Hmisc::wtd.var(x, weights = weights, na.rm = TRUE))
}



#' @rdname wtd_sd
#' @export
wtd_se <- function(x, weights = NULL) {
  # check if suggested packages are available
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop("Package `Hmisc` needed for this function to work. Please install it.", call. = FALSE)
  }

  UseMethod("wtd_se")
}


#' @export
wtd_se.data.frame <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ wtd_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
wtd_se.matrix <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ wtd_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
wtd_se.default <- function(x, weights = NULL) {
  wtd_se_helper(x, weights)
}

wtd_se_helper <- function(x, weights) {
  sqrt(Hmisc::wtd.var(x, weights = weights, na.rm = TRUE) / length(stats::na.omit(x)))
}


#' @rdname wtd_sd
#' @importFrom stats as.formula
#' @export
svy_md <- function(x, design) {
  # check if pkg survey is available
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package `survey` needed to for this function to work. Please install it.", call. = FALSE)
  }

  # deparse
  v <- stats::as.formula(paste("~", as.character(substitute(x))))

  as.vector(
    survey::svyquantile(
      v,
      design = design,
      quantiles = 0.5,
      ci = FALSE,
      na.rm = TRUE
    )
  )
}
