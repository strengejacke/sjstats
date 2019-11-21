#' @rdname weighted_sd
#' @export
weighted_median <- function(x, weights = NULL) {
  UseMethod("weighted_median")
}

#' @export
weighted_median.default <- function(x, weights = NULL) {
  weighted_md_helper(x, w = weights, p = 0.5)
}

#' @importFrom purrr map_dbl
#' @importFrom dplyr select_if
#' @export
weighted_median.data.frame <- function(x, weights = NULL) {
  dplyr::select_if(x, is.numeric) %>%
    purrr::map_dbl(~ weighted_md_helper(.x, w = weights, p = 0.5))
}

weighted_md_helper <- function(x, w, p = .5) {
  if (is.null(w)) w <- rep(1, length(x))

  x[is.na(w)] <- NA
  w[is.na(x)] <- NA

  w <- na.omit(w)
  x <- na.omit(x)

  order <- order(x)
  x <- x[order]
  w <- w[order]

  rw <- cumsum(w) / sum(w)
  md.values <- min(which(rw >= p))

  if (rw[md.values] == p)
    q <- mean(x[md.values:(md.values + 1)])
  else
    q <- x[md.values]

  q
}
