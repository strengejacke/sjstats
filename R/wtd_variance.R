weighted_variance <- function(x, w) {
  if (is.null(w)) w <- rep(1, length(x))

  x[is.na(w)] <- NA
  w[is.na(x)] <- NA

  w <- stats::na.omit(w)
  x <- stats::na.omit(x)

  xbar <- sum(w * x) / sum(w)
  sum(w * ((x - xbar)^2)) / (sum(w) - 1)
}
