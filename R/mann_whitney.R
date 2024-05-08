#' @title Mann-Whitney-U-Test
#' @name mann_whitney
#' @description This function performs a Mann-Whitney-Test (or Wilcoxon rank
#' sum test for _unpaired_ samples, see [`wilcox.test()`] and [`coin::wilcox_test()`]).
#'
#' The function reports U, p and Z-values as well as effect size r and group-rank-means.
#'
#' @param data A data frame.
#' @param select The dependent variable (numeric) to be used for the test.
#' @param by The grouping variable (factor) to be used for the test. If `by` is
#' not a factor, it will be coerced to a factor.
#' @param weights An optional weighting variable (numeric) to be used for the test.
#' @param distribution Indicates how the null distribution of the test statistic
#' should be computed. May be one of `"exact"`, `"approximate"` or `"asymptotic"`
#' (default). See [`coin::wilcox_test()`] for details.
#'
#' @return A data frame.
#'
#' @note This function calls [`coin::wilcox_test()`] to extract effect sizes.
#' Interpretation of effect sizes, as a rule-of-thumb:
#'
#' - small effect >= 0.1
#' - medium effect >= 0.3
#' - large effect >= 0.5
#'
#' @examples
#' data(efc)
#' # Mann-Whitney-U-Tests for elder's age by elder's sex.
#' mann_whitney(efc, "e17age", "e16sex")
#' @export
mann_whitney <- function(data,
                         select = NULL,
                         by = NULL,
                         weights = NULL,
                         distribution = "asymptotic") {
  insight::check_if_installed("datawizard")

  # check if "select" is in data
  if (!select %in% colnames(data)) {
    insight::format_error(
      sprintf("Variable '%s' not found in data frame.", select),
      .misspelled_string(select, colnames(data), "Maybe misspelled?")
    )
  }
  # check if "by" is in data
  if (!by %in% colnames(data)) {
    insight::format_error(
      sprintf("Variable '%s' not found in data frame.", by),
      .misspelled_string(by, colnames(data), "Maybe misspelled?")
    )
  }
  # check if "weights" is in data
  if (!is.null(weights) && !weights %in% colnames(data)) {
    insight::format_error(
      sprintf("Weighting variable '%s' not found in data frame.", weights),
      .misspelled_string(weights, colnames(data), "Maybe misspelled?")
    )
  }

  # get data
  dv <- data[[select]]
  grp <- data[[by]]

  # coerce to factor and character to numeric
  grp <- datawizard::to_factor(grp)

  # length of value range
  labels <- sjlabelled::get_labels(
    data[[by]], attr.only = FALSE, values = NULL, non.labelled = TRUE
  )

  .calculate_mwu(dv, grp, distribution, labels)
}


.calculate_mwu <- function(dv, grp, distribution, labels) {
  insight::check_if_installed("coin")
  # prepare data
  wcdat <- data.frame(dv, grp)
  # perfom wilcox test
  wt <- coin::wilcox_test(dv ~ grp, data = wcdat, distribution = distribution)

  # for rank mean
  group_levels <- levels(grp)

  # compute statistics
  u <- as.numeric(coin::statistic(wt, type = "linear"))
  z <- as.numeric(coin::statistic(wt, type = "standardized"))
  p <- coin::pvalue(wt)
  r <- abs(z / sqrt(length(dv)))
  w <- stats::wilcox.test(dv ~ grp, data = wcdat)$statistic

  # group means
  rank_mean_1 <- mean(rank(dv)[which(grp == group_levels[1])], na.rm = TRUE)
  rank_mean_2 <- mean(rank(dv)[which(grp == group_levels[2])], na.rm = TRUE)

  # compute n for each group
  n_grp1 <- length(stats::na.omit(dv[which(grp == group_levels[1])]))
  n_grp2 <- length(stats::na.omit(dv[which(grp == group_levels[2])]))

  out <- data.frame(
    group1 = group_levels[1],
    group_2 = group_levels[2],
    u = u,
    w = w,
    p = p,
    z = z,
    r = r
  )
  attr(out, "rank_means") <- stats::setNames(
    c(rank_mean_1, rank_mean_2),
    c("Mean Group 1", "Mean Group 2")
  )
  attr(out, "n_groups") <- stats::setNames(
    c(n_grp1, n_grp2),
    c("N Group 1", "N Group 2")
  )
  attr(out, "group_labels") <- labels
  out
}


.compute_weighted_mwu <- function(dv, grp, weights, labels) {
  # check if pkg survey is available
  insight::check_if_installed("survey")

  dat <- data.frame(dv, grp, weights)
  colnames(dat) <- c("x", "g", "w")

  if (insight::n_unqiue(dat$g) > 2) {
    m <- "Weighted Kruskal-Wallis test"
    method <- "KruskalWallis"
  } else {
    m <- "Weighted Mann-Whitney-U test"
    method <- "wilcoxon"
  }

  design <- survey::svydesign(ids = ~0, data = dat, weights = ~w)
  out <- survey::svyranktest(formula = x ~ g, design, test = method)

  # for rank mean
  group_levels <- levels(grp)

  design_mean1 <- survey::svydesign(
    ids = ~0,
    data = dat[dat$grp == group_levels[1]],
    weights = ~w
  )
  rank_mean_1 <- survey::svymean(~x, design_mean1)

  design_mean2 <- survey::svydesign(
    ids = ~0,
    data = dat[dat$grp == group_levels[2]],
    weights = ~w
  )
  rank_mean_2 <- survey::svymean(~x, design_mean2)

  out$method <- m
  attr(out, "rank_means") <- stats::setNames(
    c(rank_mean_1, rank_mean_2),
    c("Mean Group 1", "Mean Group 2")
  )

  out
}


.misspelled_string <- function(source, searchterm, default_message = NULL) {
  if (is.null(searchterm) || length(searchterm) < 1) {
    return(default_message)
  }
  # used for many matches
  more_found <- ""
  # init default
  msg <- ""
  # remove matching strings
  same <- intersect(source, searchterm)
  searchterm <- setdiff(searchterm, same)
  source <- setdiff(source, same)
  # guess the misspelled string
  possible_strings <- unlist(lapply(searchterm, function(s) {
    source[.fuzzy_grep(source, s)] # nolint
  }), use.names = FALSE)
  if (length(possible_strings)) {
    msg <- "Did you mean "
    if (length(possible_strings) > 1) {
      # make sure we don't print dozens of alternatives for larger data frames
      if (length(possible_strings) > 5) {
        more_found <- sprintf(
          " We even found %i more possible matches, not shown here.",
          length(possible_strings) - 5
        )
        possible_strings <- possible_strings[1:5]
      }
      msg <- paste0(msg, "one of ", toString(paste0("\"", possible_strings, "\"")))
    } else {
      msg <- paste0(msg, "\"", possible_strings, "\"")
    }
    msg <- paste0(msg, "?", more_found)
  } else {
    msg <- default_message
  }
  # no double white space
  insight::trim_ws(msg)
}


.fuzzy_grep <- function(x, pattern, precision = NULL) {
  if (is.null(precision)) {
    precision <- round(nchar(pattern) / 3)
  }
  if (precision > nchar(pattern)) {
    return(NULL)
  }
  p <- sprintf("(%s){~%i}", pattern, precision)
  grep(pattern = p, x = x, ignore.case = FALSE)
}
