#' @title Summary of mean values by group
#' @name grpmean
#'
#' @description Computes mean, sd and se for each sub-group (indicated by \code{grp})
#'                of \code{dv}.
#'
#' @param x A (grouped) data frame.
#' @param dv Name of the dependent variable, for which the mean value, grouped
#'   by \code{grp}, is computed.
#' @param grp Factor with the cross-classifying variable, where \code{dv} is
#'   grouped into the categories represented by \code{grp}. Numeric vectors
#'   are coerced to factors.
#' @param weights Name of variable in \code{x} that indicated the vector of
#'   weights that will be applied to weight all observations. Default is
#'   \code{NULL}, so no weights are used.
#' @param digits Numeric, amount of digits after decimal point when rounding
#'   estimates and values.
#' @param weight.by Deprecated.
#'
#' @inheritParams hdi
#'
#' @return For non-grouped data frames, \code{grpmean()} returns a data frame with
#'   following columns: \code{term}, \code{mean}, \code{N}, \code{std.dev},
#'   \code{std.error} and \code{p.value}. For grouped data frames, returns
#'   a list of such data frames.
#'
#' @details This function performs a One-Way-Anova with \code{dv} as dependent
#'   and \code{grp} as independent variable, by calling
#'   \code{lm(count ~ as.factor(grp))}. Then \code{\link[emmeans]{contrast}}
#'   is called to get p-values for each sub-group. P-values indicate whether
#'   each group-mean is significantly different from the total mean.
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
#' # weighting
#' efc$weight <- abs(rnorm(n = nrow(efc), mean = 1, sd = .5))
#' grpmean(efc, c12hour, e42dep, weights = weight)
#'
#' @importFrom sjlabelled get_label drop_labels get_labels
#' @importFrom stats lm na.omit sd weighted.mean
#' @importFrom purrr map_chr map_df
#' @importFrom sjmisc to_value is_empty
#' @importFrom rlang enquo .data quo_name
#' @export
grpmean <- function(x, dv, grp, weights = NULL, digits = 2, out = c("txt", "viewer", "browser"), weight.by) {

  out <- match.arg(out)

  if (out != "txt" && !requireNamespace("sjPlot", quietly = TRUE)) {
    message("Package `sjPlot` needs to be loaded to print HTML tables.")
    out <- "txt"
  }

  ## TODO remove deprecated argument later

  if (!missing(weight.by)) {
    message("Argument `weight.by` is deprecated. Please use `weights`.")
    weights <- weight.by
  }

  # create quosures
  grp.name <- rlang::quo_name(rlang::enquo(grp))
  dv.name <- rlang::quo_name(rlang::enquo(dv))

  # weights need extra checking, might be NULL
  if (!missing(weights)) {
    .weights <- rlang::quo_name(rlang::enquo(weights))

    w.string <- tryCatch(
      {
        eval(weights)
      },
      error = function(x) { NULL },
      warning = function(x) { NULL },
      finally = function(x) { NULL }
    )

    if (!is.null(w.string) && is.character(w.string)) .weights <- w.string
    if (sjmisc::is_empty(.weights) || .weights == "NULL") .weights <- NULL
  } else
    .weights <- NULL


  # create string with variable names
  vars <- c(grp.name, dv.name, .weights)

  # get data
  x <- suppressMessages(dplyr::select(x, !! vars))

  # set value and row labels
  varGrpLabel <- sjlabelled::get_label(x[[grp.name]], def.value = grp.name)
  varCountLabel <- sjlabelled::get_label(x[[dv.name]], def.value = dv.name)

  # first, drop unused labels
  x[[grp.name]] <- sjlabelled::drop_labels(x[[grp.name]], drop.na = TRUE)

  # now get valid value labels
  value.labels <- sjlabelled::get_labels(
    x[[grp.name]], attr.only = F, values = "n", non.labelled = TRUE
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
        weight.by = .weights,
        digits = digits,
        value.labels = value.labels,
        varCountLabel = varCountLabel,
        varGrpLabel = varGrpLabel
      )

      attr(dummy, "group") <- get_grouped_title(x, grps, i, sep = "\n")

      # save data frame for return value
      dataframes[[length(dataframes) + 1]] <- dummy
    }

    # add class-attr for print-method()
    if (out == "txt")
      class(dataframes) <- c("sj_grpmeans", "list")
    else
      class(dataframes) <- c("sjt_grpmeans", "list")

  } else {
    dataframes <- grpmean_helper(
      x = x,
      dv = dv.name,
      grp = grp.name,
      weight.by = .weights,
      digits = digits,
      value.labels = value.labels,
      varCountLabel = varCountLabel,
      varGrpLabel = varGrpLabel
    )

    # add class-attr for print-method()
    if (out == "txt")
      class(dataframes) <- c("sj_grpmean", class(dataframes))
    else
      class(dataframes) <- c("sjt_grpmean", class(dataframes))
  }

  # save how to print output
  attr(dataframes, "print") <- out

  dataframes
}


#' @importFrom stats pf lm weighted.mean na.omit sd
#' @importFrom sjmisc to_value
#' @importFrom emmeans emmeans contrast
#' @importFrom dplyr pull select n_distinct
#' @importFrom purrr map_chr
#' @importFrom rlang .data
grpmean_helper <- function(x, dv, grp, weight.by, digits, value.labels, varCountLabel, varGrpLabel) {
  # copy vectors from data frame
  dv <- x[[dv]]
  grp <- x[[grp]]

  if (!is.null(weight.by))
    weight.by <- x[[weight.by]]
  else
    weight.by <- 1

  # convert values to numeric
  dv <- sjmisc::to_value(dv)

  # create data frame, for emmeans
  mydf <- stats::na.omit(data.frame(
    dv = dv,
    grp = as.factor(grp),
    weight.by = weight.by
  ))

  # compute anova statistics for mean table
  fit <- stats::lm(dv ~ grp, weights = weight.by, data = mydf)

  # p-values of contrast-means
  means.p <- fit %>%
    emmeans::emmeans(specs = "grp") %>%
    emmeans::contrast(method = "eff") %>%
    summary() %>%
    dplyr::pull("p.value")

  # create string with p-values
  pval <- purrr::map_chr(seq_len(length(means.p)), function(i) {
    if (means.p[i] < 0.001) {
      "<0.001"
    } else {
      sprintf("%.*f", digits, means.p[i])
    }
  })


  ## TODO
  # efc %>%
  #   group_by(c172code, c161sex) %>%
  #   grpmean(c12hour, e42dep)


  # check if value labels length matches group count
  if (dplyr::n_distinct(mydf$grp) != length(value.labels)) {
    # get unique factor levels and check if these are numeric.
    # if so, we match the values from value labels and the remaining
    # factor levels, so we get the correct value labels for printing
    nl <- unique(mydf$grp)
    if (sjmisc::is_num_fac(nl))
      value.labels <- value.labels[names(value.labels) %in% levels(nl)]
    else
      value.labels <- nl
  }


  # create summary
  dat <- mydf %>%
    dplyr::group_by(.data$grp) %>%
    summarise(
      mean = sprintf("%.*f", digits, stats::weighted.mean(.data$dv, w = .data$weight.by, na.rm = TRUE)),
      N = round(sum(.data$weight.by)),
      std.dev = sprintf("%.*f", digits, wtd_sd(.data$dv, .data$weight.by)),
      std.error = sprintf("%.*f", digits, wtd_se(.data$dv, .data$weight.by))
    ) %>%
    mutate(p.value = pval) %>%
    dplyr::select(-.data$grp)

  # finally, add total-row
  dat <- dplyr::bind_rows(
    dat,
    data_frame(
      mean = sprintf("%.*f", digits, stats::weighted.mean(mydf$dv, w = mydf$weight.by, na.rm = TRUE)),
      N = nrow(mydf),
      std.dev = sprintf("%.*f", digits, wtd_sd(mydf$dv, mydf$weight.by)),
      std.error = sprintf("%.*f", digits, wtd_se(mydf$dv, mydf$weight.by)),
      p.value = ""
    )
  )


  # add row labels
  dat <- add_cols(
    dat,
    term = c(unname(value.labels), "Total"),
    .after = -1
  )


  # get anova statistics for mean table
  sum.fit <- summary(fit)

  # r-squared values
  r2 <- sum.fit$r.squared
  r2.adj <- sum.fit$adj.r.squared

  # F-statistics
  fstat <- sum.fit$fstatistic
  pval <- stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = F)


  # copy as attributes
  attr(dat, "r2") <- r2
  attr(dat, "adj.r2") <- r2.adj
  attr(dat, "fstat") <- fstat[1]
  attr(dat, "p.value") <- pval
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
