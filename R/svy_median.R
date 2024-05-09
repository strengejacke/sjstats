#' @rdname weighted_se
#' @export
survey_median <- function(x, design) {
  # check if pkg survey is available
  insight::check_if_installed("suvey")

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
