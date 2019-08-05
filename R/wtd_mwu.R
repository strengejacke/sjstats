#' @rdname wtd_sd
#' @export
wtd_mwu <- function(data, ...) {
  UseMethod("wtd_mwu")
}


#' @importFrom dplyr select
#' @rdname wtd_sd
#' @export
wtd_mwu.default <- function(data, x, grp, weights, ...) {
  x.name <- deparse(substitute(x))
  g.name <- deparse(substitute(grp))
  w.name <- deparse(substitute(weights))

  # create string with variable names
  vars <- c(x.name, g.name, w.name)

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  wtd_mwu_helper(dat)
}


#' @importFrom dplyr select
#' @rdname wtd_sd
#' @export
wtd_mwu.formula <- function(formula, data, ...) {
  vars <- all.vars(formula)

  # get data
  dat <- suppressMessages(dplyr::select(data, !! vars))
  dat <- na.omit(dat)

  wtd_mwu_helper(dat)
}

wtd_mwu_helper <- function(dat, vars) {
  # check if pkg survey is available
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package `survey` needed to for this function to work. Please install it.", call. = FALSE)
  }

  x.name <- colnames(dat)[1]
  group.name <- colnames(dat)[2]

  colnames(dat) <- c("x", "g", "w")

  design <- survey::svydesign(ids = ~0, data = dat, weights = ~w)
  mw <- survey::svyranktest(formula = x ~ g, design)

  attr(mw, "x.name") <- x.name
  attr(mw, "group.name") <- group.name
  class(mw) <- c("sj_wmwu", "list")

  if (dplyr::n_distinct(dat$g, na.rm = TRUE) > 2)
    m <- "Weighted Kruskal-Wallis test"
  else
    m <- "Weighted Mann-Whitney-U test"

  mw$method <- m

  mw
}
