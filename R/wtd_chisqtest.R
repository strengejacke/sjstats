#' @rdname weighted_sd
#' @export
weighted_chisqtest <- function(data, ...) {
  UseMethod("weighted_chisqtest")
}

#' @importFrom dplyr select
#' @importFrom stats na.omit chisq.test as.formula
#' @rdname weighted_sd
#' @export
weighted_chisqtest.default <- function(data, x, y, weights, ...) {
  x.name <- deparse(substitute(x))
  y.name <- deparse(substitute(y))
  w.name <- deparse(substitute(weights))

  if (w.name == "NULL") {
    w.name <- "weights"
    data$weights <- 1
  }

  # create string with variable names
  vars <- .compact_character(c(x.name, y.name, w.name))

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- stats::na.omit(dat)

  colnames(dat)[ncol(dat)] <- ".weights"

  # check if we have chisq-test for given probabilities
  dot_args <- list(...)
  if ("p" %in% names(dot_args)) {
    .weighted_chisq_for_prob(dat, x.name, prob = dot_args[["p"]])
  } else {
    crosstable_statistics(data = dat, statistics = "auto", weights = ".weights", ...)
  }
}


#' @importFrom stats xtabs
#' @rdname weighted_sd
#' @export
weighted_chisqtest.formula <- function(formula, data, ...) {
  vars <- all.vars(formula)
  dot_args <- list(...)

  if (length(vars) < 3 && !"p" %in% names(dot_args)) {
    vars <- c(vars, ".weights")
    data$.weights <- 1
  }

  if ("p" %in% names(dot_args)) {
    dat <- data[vars]
    colnames(dat)[ncol(dat)] <- ".weights"
    .weighted_chisq_for_prob(dat, names(dat)[1], prob = dot_args[["p"]])
  } else {
    tab <- as.table(round(stats::xtabs(data[[vars[3]]] ~ data[[vars[1]]] + data[[vars[2]]])))
    class(tab) <- "table"
    crosstable_statistics(data = tab, statistics = "auto", weights = NULL, ...)
  }
}



.weighted_chisq_for_prob <- function(dat, x.name, prob) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package `survey` needed to for this function to work. Please install it.", call. = FALSE)
  }

  if (abs(sum(prob) - 1) > sqrt(.Machine$double.eps)) {
    prob <- prob / sum(prob)
  }

  dat$sj_subject_id <- 1:nrow(dat)
  dat$sj_weights <- dat$.weights
  design <- survey::svydesign(id = ~sj_subject_id, weights = ~sj_weights, data = dat)
  stable <- survey::svytable(stats::as.formula(paste0("~", x.name)), design)
  out <- stats::chisq.test(stable, p = prob)

  structure(class = "sj_xtab_stat2", list(
    estimate = out$statistic,
    p.value = out$p.value,
    stat.name = "Chi-squared",
    stat.html = "&chi;<sup>2</sup>",
    df = out$parameter,
    method = "Weighted chi-squared test for given probabilities",
    method.html = "Weighted &chi;<sup>2</sup> test for given probabilities",
    method.short = "Chi-squared"
  ))
}
