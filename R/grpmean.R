#' @title Summary of mean values by group
#' @name grpmean
#'
#' @description Computes mean, sd and se for each sub-group (indicated by \code{grp})
#'                of \code{count}.
#'
#' @param x A data frame.
#' @param dv Name of the dependent variable, for which the mean value, grouped
#'        by \code{grp}, is computed.
#' @param grp Factor with the cross-classifying variable, where \code{dv} is
#'        grouped into the categories represented by \code{grp}.
#' @param weight.by Vector of weights that will be applied to weight all cases.
#'        Must be a vector of same length as the input vector. Default is
#'        \code{NULL}, so no weights are used.
#' @param digits Numeric, amount of digits after decimal point when rounding estimates and values.
#'
#' @return ...
#'
#' @details This function performs a One-Way-Anova with \code{count} as dependent
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
#' @importFrom sjlabelled get_label drop_labels get_labels
#' @importFrom stats lm na.omit sd weighted.mean
#' @importFrom purrr map_chr map_df
#' @importFrom tibble tibble add_row add_column
#' @importFrom sjmisc to_value
#' @export
grpmean <- function(x, dv, grp, weight.by = NULL, digits = 2) {

  # get names of variables. used, if these are not labelled
  grp.name <- deparse(substitute(grp))
  dv.name <- deparse(substitute(dv))

  # copy vectors from data frame
  dv <- x[[dv.name]]
  grp <- x[[grp.name]]
  weight.by <- x[[deparse(substitute(weight.by))]]

  # set value and row labels
  varGrpLabel <- sjlabelled::get_label(grp, def.value = grp.name)
  varCountLabel <- sjlabelled::get_label(dv, def.value = dv.name)

  # first, drop unused labels
  grp <- sjlabelled::drop_labels(grp, drop.na = TRUE)
  # now get valid value labels
  value.labels <- sjlabelled::get_labels(
    grp, attr.only = F, include.values = NULL, include.non.labelled = TRUE
  )

  # convert values to numeric
  dv <- sjmisc::to_value(dv)

  # compute anova statistics for mean table
  if (!is.null(weight.by)) {
    fit <- stats::lm(dv ~ as.factor(grp), weights = weight.by)
  } else {
    fit <- stats::lm(dv ~ as.factor(grp))
  }

  # p-values of means
  means.p <- get_model_pval(fit)[["p.value"]]

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

  class(dat) <- c("sj_grpmean", class(dat))
  dat
}
