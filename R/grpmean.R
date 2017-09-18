#' @title Summary of mean values by group
#' @name grpmean
#'
#' @description Computes mean, sd and se for each sub-group (indicated by \code{grp})
#'                of \code{dv}.
#'
#' @param x A (grouped) data frame.
#' @param dv Name of the dependent variable, for which the mean value, grouped
#'        by \code{grp}, is computed.
#' @param grp Factor with the cross-classifying variable, where \code{dv} is
#'        grouped into the categories represented by \code{grp}. Numeric vectors
#'        are coerced to factors.
#' @param weight.by Vector of weights that will be applied to weight all cases.
#'        Must be a vector of same length as the input vector. Default is
#'        \code{NULL}, so no weights are used.
#' @param digits Numeric, amount of digits after decimal point when rounding estimates and values.
#'
#' @return For non-grouped data frames, \code{grpmean()} returns a data frame with
#'         following columns: \code{term}, \code{mean}, \code{N}, \code{std.dev},
#'         \code{std.error} and \code{p.value}. For grouped data frames, returns
#'         a list of such data frames.
#'
#' @details This function performs a One-Way-Anova with \code{dv} as dependent
#'            and \code{grp} as independent variable, by calling
#'            \code{lm(count ~ as.factor(grp))}, to get p-values for each
#'            sub-group and the complete "model". P-values indicate whether
#'            each group-mean is significantly different from the reference
#'            group (reference level of \code{grp}).
#'
#' @examples
#' data(efc)
#' grpmean(efc, c12hour, e42dep)
#'
#' data(iris)
#' grpmean(iris, Sepal.Width, Species)
#'
#' # also works for grouped data frames
#' library(dplyr)
#' efc %>%
#'   group_by(c172code) %>%
#'   grpmean(c12hour, e42dep)
#'
#' @importFrom sjlabelled get_label drop_labels get_labels
#' @importFrom stats lm na.omit sd weighted.mean
#' @importFrom purrr map_chr map_df
#' @importFrom tibble tibble add_row add_column
#' @importFrom sjmisc to_value
#' @importFrom rlang enquo .data quo_name
#' @export
grpmean <- function(x, dv, grp, weight.by = NULL, digits = 2) {
  # create quosures
  grp.name <- rlang::quo_name(rlang::enquo(grp))
  dv.name <- rlang::quo_name(rlang::enquo(dv))

  # weights need extra checking, might be NULL
  if (!is.null(weight.by))
    weights <- rlang::quo_name(rlang::enquo(weight.by))
  else
    weights <- NULL

  # create string with variable names
  vars <- c(grp.name, dv.name, weights)

  # get data
  x <- suppressMessages(dplyr::select(x, !! vars))

  # set value and row labels
  varGrpLabel <- sjlabelled::get_label(x[[grp.name]], def.value = grp.name)
  varCountLabel <- sjlabelled::get_label(x[[dv.name]], def.value = dv.name)

  # first, drop unused labels
  x[[grp.name]] <- sjlabelled::drop_labels(x[[grp.name]], drop.na = TRUE)

  # now get valid value labels
  value.labels <- sjlabelled::get_labels(
    x[[grp.name]], attr.only = F, include.values = NULL, include.non.labelled = TRUE
  )

  # return values
  dataframes <- list()

  # do we have a grouped data frame?
  if (inherits(x, "grouped_df")) {
    # get grouped data
    grps <- get_grouped_data(x)

    # now plot everything
    for (i in seq_len(nrow(grps))) {
      # copy back labels to grouped data frame
      tmp <- sjlabelled::copy_labels(grps$data[[i]], x)

      # get grouped means table
      dummy <- grpmean_helper(
        x = tmp,
        dv = dv.name,
        grp = grp.name,
        weight.by = weights,
        digits = digits,
        value.labels = value.labels,
        varCountLabel = varCountLabel,
        varGrpLabel = varGrpLabel
      )

      attr(dummy, "group") <- get_grouped_title(x, grps, i, sep = "\n")

      # save data frame for return value
      dataframes[[length(dataframes) + 1]] <- dummy

      # add class-attr for print-method()
      class(dataframes) <- c("sj_grpmeans", "list")
    }
  } else {
    dataframes <- grpmean_helper(
      x = x,
      dv = dv.name,
      grp = grp.name,
      weight.by = weights,
      digits = digits,
      value.labels = value.labels,
      varCountLabel = varCountLabel,
      varGrpLabel = varGrpLabel
    )

    # add class-attr for print-method()
    class(dataframes) <- c("sj_grpmean", class(dataframes))
  }

  dataframes
}


grpmean_helper <- function(x, dv, grp, weight.by, digits, value.labels, varCountLabel, varGrpLabel) {
  # copy vectors from data frame
  dv <- x[[dv]]
  grp <- x[[grp]]
  if (!is.null(weight.by)) weight.by <- x[[weight.by]]

  # convert values to numeric
  dv <- sjmisc::to_value(dv)

  # compute anova statistics for mean table
  if (!is.null(weight.by)) {
    fit <- stats::lm(dv ~ as.factor(grp), weights = weight.by)
  } else {
    fit <- stats::lm(dv ~ as.factor(grp))
  }

  # p-values of means
  means.p <- p_value(fit)[["p.value"]]

  # create string with p-values
  pval <- purrr::map_chr(seq_len(length(means.p)), function(i) {
    if (means.p[i] < 0.001) {
      "<0.001"
    } else {
      sprintf("%.*f", digits, means.p[i])
    }
  })

  # retrieve group indices
  indices <- sort(unique(stats::na.omit(grp)))

  dat <- purrr::map_df(seq_len(length(indices)), function(i) {
    # do we have weighted means?
    if (!is.null(weight.by)) {
      mw <- stats::weighted.mean(
        dv[grp == indices[i]],
        w = weight.by[grp == indices[i]],
        na.rm = TRUE
      )
    } else {
      mw <- mean(dv[grp == indices[i]], na.rm = TRUE)
    }

    # add new row to data frame with mean, N, sd and se of dv
    # for each sub-group (indicated by indices)
    tibble::tibble(
      mean = sprintf("%.*f", digits, mw),
      N = length(stats::na.omit(dv[grp == indices[i]])),
      std.dev = sprintf("%.*f", digits, stats::sd(dv[grp == indices[i]], na.rm = TRUE)),
      std.error = sprintf("%.*f", digits, se(dv[grp == indices[i]])),
      p.value = pval[i]
    )
  })

  # total mean
  if (!is.null(weight.by)) {
    mw <- weighted.mean(dv, w = weight.by, na.rm = TRUE)
  } else {
    mw <- mean(dv, na.rm = TRUE)
  }

  # finally, add total-row
  dat <- tibble::add_row(
    dat,
    mean = sprintf("%.*f", digits, mw),
    N = length(stats::na.omit(dv)),
    std.dev = sprintf("%.*f", digits, stats::sd(dv, na.rm = TRUE)),
    std.error = sprintf("%.*f", digits, se(dv)),
    p.value = ""
  )

  # add row labels
  dat <- tibble::add_column(
    dat,
    term = c(value.labels, "Total"),
    .before = 1
  )

  # get anova statistics for mean table
  sum.fit <- summary(fit)
  # r-squared values
  r2 <- sum.fit$r.squared
  r2.adj <- sum.fit$adj.r.squared
  # F-statistics
  fstat <- sum.fit$fstatistic[1]

  # copy as attributes
  attr(dat, "r2") <- r2
  attr(dat, "adj.r2") <- r2.adj
  attr(dat, "fstat") <- fstat
  attr(dat, "dv.label") <- varCountLabel
  attr(dat, "grp.label") <- varGrpLabel

  dat
}


get_grouped_title <- function(x, grps, i, sep = "\n") {
  # create title for first grouping level
  tp <- get_title_part(x, grps, 1, i)
  title <- sprintf("%s: %s", tp[1], tp[2])

  # do we have another groupng variable?
  if (length(attr(x, "vars", exact = T)) > 1) {
    tp <- get_title_part(x, grps, 2, i)
    title <- sprintf("%s%s%s: %s", title, sep, tp[1], tp[2])
  }

  # return title
  title
}


get_title_part <- function(x, grps, level, i) {
  # prepare title for group
  var.name <- colnames(grps)[level]

  # get values from value labels
  vals <- sjlabelled::get_values(x[[var.name]])
  # if we have no value labels, get values directly
  if (is.null(vals)) vals <- unique(x[[var.name]])
  # find position of value labels for current group
  lab.pos <- which(vals == grps[[var.name]][i])

  # get variable and value labels
  t1 <- sjlabelled::get_label(x[[var.name]], def.value = var.name)
  t2 <- sjlabelled::get_labels(x[[var.name]])[lab.pos]

  # if we have no value label, use value instead
  if (is.null(t2)) t2 <- vals[lab.pos]

  # generate title
  c(t1, t2)
}
