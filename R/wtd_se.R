#' @rdname weighted_sd
#' @export
weighted_se <- function(x, weights = NULL) {
  UseMethod("weighted_se")
}


#' @export
weighted_se.data.frame <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ weighted_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
weighted_se.matrix <- function(x, weights = NULL) {
  se_result <- purrr::map_dbl(x, ~ weighted_se_helper(.x, weights = weights))
  names(se_result) <- colnames(x)

  se_result
}

#' @export
weighted_se.default <- function(x, weights = NULL) {
  weighted_se_helper(x, weights)
}

weighted_se_helper <- function(x, weights) {
  if (is.null(weights)) weights <- rep(1, length(x))
  sqrt(weighted_variance(x, weights) / length(stats::na.omit(x)))
}
