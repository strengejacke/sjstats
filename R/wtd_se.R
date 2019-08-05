#' @rdname wtd_sd
#' @export
wtd_se <- function(x, weights = NULL) {
  UseMethod("wtd_se")
}


#' @export
wtd_se.data.frame <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ wtd_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
wtd_se.matrix <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ wtd_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
wtd_se.default <- function(x, weights = NULL) {
  wtd_se_helper(x, weights)
}

wtd_se_helper <- function(x, weights) {
  if (is.null(weights)) weights <- rep(1, length(x))
  sqrt(wtd_var(x, weights) / length(stats::na.omit(x)))
}
