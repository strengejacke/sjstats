#' @rdname weighted_sd
#' @export
weighted_mannwhitney <- function(data, ...) {
  UseMethod("weighted_mannwhitney")
}


#' @importFrom dplyr select
#' @rdname weighted_sd
#' @export
weighted_mannwhitney.default <- function(data, x, grp, weights, ...) {
  x.name <- deparse(substitute(x))
  g.name <- deparse(substitute(grp))
  w.name <- deparse(substitute(weights))

  # create string with variable names
  vars <- c(x.name, g.name, w.name)

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  weighted_mannwhitney_helper(dat)
}


#' @importFrom dplyr select
#' @rdname weighted_sd
#' @export
weighted_mannwhitney.formula <- function(formula, data, ...) {
  vars <- all.vars(formula)

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  weighted_mannwhitney_helper(dat)
}

weighted_mannwhitney_helper <- function(dat, vars) {
  # check if pkg survey is available
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package `survey` needed to for this function to work. Please install it.", call. = FALSE)
  }

  x.name <- colnames(dat)[1]
  group.name <- colnames(dat)[2]

  colnames(dat) <- c("x", "g", "w")

  if (dplyr::n_distinct(dat$g, na.rm = TRUE) > 2) {
    m <- "Weighted Kruskal-Wallis test"
    method <- "KruskalWallis"
  } else {
    m <- "Weighted Mann-Whitney-U test"
    method <- "wilcoxon"
  }

  design <- survey::svydesign(ids = ~0, data = dat, weights = ~w)
  mw <- survey::svyranktest(formula = x ~ g, design, test = method)

  attr(mw, "x.name") <- x.name
  attr(mw, "group.name") <- group.name
  class(mw) <- c("sj_wmwu", "list")

  mw$method <- m

  mw
}
