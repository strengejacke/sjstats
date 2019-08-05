wtd_var <- function(x, w) {
  if (is.null(w)) w <- rep(1, length(x))

  x[is.na(w)] <- NA
  w[is.na(x)] <- NA

  w <- na.omit(w)
  x <- na.omit(x)

  xbar <- sum(w * x) / sum(w)
  sum(w * ((x - xbar)^2)) / (sum(w) - 1)
}
