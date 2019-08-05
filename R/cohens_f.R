#' @rdname eta_sq
#' @export
cohens_f <- function(model) {
  es <- aov_stat(model, type = "cohens.f")

  data_frame(
    term = names(es),
    cohens.f = es
  )
}
