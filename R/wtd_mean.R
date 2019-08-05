#' @rdname wtd_sd
#' @export
wtd_mean <- function(x, weights = NULL) {
  UseMethod("wtd_mean")
}

#' @importFrom stats weighted.mean
#' @export
wtd_mean.default <- function(x, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  stats::weighted.mean(x, w = weights, na.rm = TRUE)
}

#' @importFrom stats weighted.mean
#' @importFrom purrr map_dbl
#' @importFrom dplyr select_if
#' @export
wtd_mean.data.frame <- function(x, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  dplyr::select_if(x, is.numeric) %>%
    purrr::map_dbl(~ weighted.mean(.x, w = weights))
}
