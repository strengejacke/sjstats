#' @title Student's t test
#' @name t_test
#' @description This function performs a Mann-Whitney-Test (or Wilcoxon rank
#' sum test for _unpaired_ samples.
#'
#' A Mann-Whitney-Test is a non-parametric test for the null hypothesis that two
#' independent samples have identical continuous distributions. It can be used
#' when the two continuous variables are not normally distributed.
#'
#' @param data A data frame.
#' @param select Name of the dependent variable (as string) to be used for the
#' test. `select` can also be a character vector, specifying the names of
#' multiple continuous variables. In this case, `by` is ignored and variables
#' specified in `select` are used to compute the test. This can be useful if
#' the data is in wide-format and no grouping variable is available.
#' @param by Name of the grouping variable to be used for the test. If `by` is
#' not a factor, it will be coerced to a factor. For `chi_squared_test()`, if
#' `probabilities` is provided, `by` must be `NULL`.
#' @param weights Name of an (optional) weighting variable to be used for the test.
#' @param alternative A character string specifying the alternative hypothesis,
#' must be one of `"two.sided"` (default), `"greater"` or `"less"`. See `?t.test()`.
#'
#' @return A data frame with test results. The function returns p and Z-values
#' as well as effect size r and group-rank-means.
#'
#' @examples
#' data(efc)
#' # Mann-Whitney-U-Tests for elder's age by elder's sex.
#' t_test(efc, "e17age", by = "e16sex")
#'
#' # when data is in wide-format, specify all relevant continuous
#' # variables in `select` and omit `by`
#' set.seed(123)
#' wide_data <- data.frame(scale1 = runif(20), scale2 = runif(20))
#' t_test(wide_data, select = c("scale1", "scale2"))
#'
#' # same as if we had data in long format, with grouping variable
#' long_data <- data.frame(
#'   scales = c(wide_data$scale1, wide_data$scale2),
#'   groups = rep(c("A", "B"), each = 20)
#' )
#' t_test(long_data, select = "scales", by = "groups")
#' @export
t_test <- function(data,
                   select = NULL,
                   by = NULL,
                   weights = NULL,
                   paired = FALSE,
                   mu = 0,
                   alternative = "two.sided") {
  insight::check_if_installed("datawizard")
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))

  # sanity checks
  .sanitize_htest_input(data, select, by, weights)

  # does select indicate more than one variable?
  if (length(select) > 1) {
    # sanity check - may only specify two variable names
    if (length(select) > 2) {
      insight::format_error("You may only specify two variables for Student's t test.")
    }
    if (!is.null(by)) {
      insight::format_error("If `select` specifies more than one variable, `by` must be `NULL`.")
    }
    # we convert the data into long format, and create a grouping variable
    data <- datawizard::data_to_long(data[select], names_to = "group", values_to = "scale")
    by <- select[2]
    select <- select[1]
    # after converting to long, we have the "grouping" variable first in the data
    colnames(data) <- c(by, select)
  }

  # get data
  dv <- data[[select]]

  # for two-sample t-test...
  if (!is.null(by)) {
    grp <- data[[by]]
    # coerce to factor
    grp <- datawizard::to_factor(grp)
    # only two groups allowed
    if (insight::n_unique(grp) > 2) {
      insight::format_error("Only two groups are allowed for Student's t test.") # nolint
    }
    # value labels
    group_labels <- names(attr(data[[by]], "labels", exact = TRUE))
    if (is.null(group_labels)) {
      group_labels <- levels(droplevels(grp))
    }
    data_name <- paste(select, "by", by)
  } else {
    group_labels <- select
    data_name <- select
  }

  if (is.null(weights)) {
    .calculate_ttest(dv, grp, mu, paired, alternative, group_labels, data_name)
  } else {
    .calculate_weighted_ttest(dv, grp, mu, paired, alternative, data[[weights]], group_labels, data_name)
  }
}


# Mann-Whitney-Test for two groups --------------------------------------------

.calculate_ttest <- function(dv, grp, mu, alternative, group_labels, data_name) {
  insight::check_if_installed("effectsize")
  # prepare data
  if (is.null(grp)) {
    tdat <- data.frame(dv)
    t_formula <- as.formula("dv ~ 1")
  } else {
    tdat <- data.frame(dv, grp)
    t_formula <- as.formula("dv ~ grp")
  }
  # perfom wilcox test
  htest <- stats::t.test(
    t_formula,
    data = tdat,
    alternative = alternative,
    mu = mu,
    paired = paired
  )
  test_statistic <- htest$statistic
  if (nrow(tdat) > 20) {
    effect_size <- stats::setNames(
      effectsize::cohens_d(
        t_formula,
        data = tdat,
        alternative = alternative,
        mu = mu,
        paired = paired
      )$Cohens_d,
      "Cohens_d"
    )
  } else {
    effect_size <- stats::setNames(
      effectsize::hedges_g(
        t_formula,
        data = tdat,
        alternative = alternative,
        mu = mu,
        paired = paired
      )$Hedges_g,
      "Hedges_g"
    )
  }

  # return result
  out <- data.frame(
    data = data_name,
    statistic_name = "t",
    statistic = test_statistic,
    effect_size_name = names(effect_size),
    effect_size = as.numeric(effect_size),
    p = as.numeric(htest$p.value),
    df = as.numeric(htest$parameter),
    method = htest$method,
    alternative = alternative,
    stringsAsFactors = FALSE
  )
  class(out) <- c("sj_htest_t", "data.frame")
  attr(out, "means") <- stats::setNames(htest$estimate, c("Mean Group 1", "Mean Group 2"))
  attr(out, "n_groups") <- stats::setNames(
    c(as.numeric(table(grp))),
    c("N Group 1", "N Group 2")
  )
  attr(out, "paired") <- isTRUE(paired)
  attr(out, "weighted") <- FALSE
  out
}


# Weighted Mann-Whitney-Test for two groups ----------------------------------

.calculate_weighted_ttest <- function(dv, grp, mu, paired, alternative, weights, group_labels, data_name) {
  insight::check_if_installed("datawizard")
  if (is.null(grp)) {
    x_values <- dv
    x_weights <- weights
    y_values <- NULL
  } else {
    dat <- stats::na.omit(data.frame(dv, grp, weights))
    colnames(dat) <- c("y", "g", "w")
    # unique groups
    groups <- unique(dat$grp)
    # values for sample 1
    x_values <- dat$y[dat$g == groups[1]]
    x_weights <- dat$w[dat$g == groups[1]]
    # values for sample 2
    y_values <- dat$y[dat$g == groups[2]]
    y_weights <- dat$w[dat$g == groups[2]]
    # paired t-test?
    if (paired) {
      x_values <- x_values - y_values
      y_values <- NULL
    }
  }

  mu_x <- stats::weighted.mean(x_values, x_weights)
  var_x <- datawizard::weighted_sd(x_values, x_weights)^2
  se_x <- sqrt(var_x / length(x_values))

  if (paired || is.null(y_values)) {
    # paired
    se <- se_x
    dof <- length(x_values) - 1
    test_statistic <- (mu_x - mu) / se
    estimate <- stats::setNames(mu_x, if (paired) "mean of the differences" else "mean of x")
    method <- if (paired) "Paired t-test" else "One Sample t-test"
  } else {
    # unpaired t-test
    mu_y <- stats::weighted.mean(y_values, y_weights)
    var_y <- datawizard::weighted_sd(y_values, y_weights)^2
    se_y <- sqrt(var_y / length(y_values))

    se <- sqrt(se_x^2 + se_y^2)
    dof <- se^4 / (se_x^4 / (length(x_values) - 1) + se_y^4 / (length(y_values) - 1))
    test_statistic <- (mu_x - mu_y - mu) / se

    estimate <- c(mu_x, mu_y)
    names(estimate) <- c("mean of x", "mean of y")
    method <- "Two-Sample t-test"
  }

  if (alternative == "less") {
    pval <- stats::pt(test_statistic, dof)
  } else if (alternative == "greater") {
    pval <- stats::pt(test_statistic, dof, lower.tail = FALSE)
  } else {
    pval <- 2 * stats::pt(-abs(test_statistic), dof)
  }

  # return result
  out <- data.frame(
    data = data_name,
    statistic_name = "t",
    statistic = test_statistic,
    effect_size_name = names(effect_size),
    effect_size = as.numeric(effect_size),
    p = pval,
    df = dof,
    method = method,
    alternative = alternative,
    stringsAsFactors = FALSE
  )
  class(out) <- c("sj_htest_t", "data.frame")
  attr(out, "means") <- stats::setNames(htest$estimate, c("Mean Group 1", "Mean Group 2"))
  attr(out, "n_groups") <- stats::setNames(
    c(as.numeric(table(grp))),
    c("N Group 1", "N Group 2")
  )
  attr(out, "paired") <- isTRUE(paired)
  attr(out, "weighted") <- FALSE
  out
}


# methods ---------------------------------------------------------------------

#' @export
print.sj_htest_t <- function(x, ...) {
  # fetch attributes
  group_labels <- attributes(x)$group_labels
  rank_means <- attributes(x)$rank_means
  n_groups <- attributes(x)$n_groups
  weighted <- attributes(x)$weighted

  if (weighted) {
    weight_string <- " (weighted)"
  } else {
    weight_string <- ""
  }

  # same width
  group_labels <- format(group_labels)

  # header
  insight::print_color(sprintf("# Mann-Whitney test%s\n\n", weight_string), "blue")

  # group-1-info
  insight::print_color(
    sprintf(
      "  Group 1: %s (n = %i, rank mean = %s)\n",
      group_labels[1], n_groups[1], insight::format_value(rank_means[1], protect_integers = TRUE)
    ), "cyan"
  )

  # group-2-info
  insight::print_color(
    sprintf(
      "  Group 2: %s (n = %i, rank mean = %s)\n",
      group_labels[2], n_groups[2], insight::format_value(rank_means[2], protect_integers = TRUE)
    ), "cyan"
  )

  cat(sprintf("\n  r = %.3f, Z = %.3f, %s\n\n", x$r, x$z, insight::format_p(x$p)))
}
