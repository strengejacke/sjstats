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


#' @rdname wtd_sd
#' @export
survey_median <- svy_md
