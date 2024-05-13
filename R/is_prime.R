#' @title Find prime numbers
#' @name is_prime
#'
#' @description This functions checks whether a number is, or numbers in a
#' vector are prime numbers.
#'
#' @param x An integer, or a vector of integers.
#'
#' @return `TRUE` for each prime number in `x`, `FALSE` otherwise.
#'
#' @examples
#' is_prime(89)
#' is_prime(15)
#' is_prime(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
#'
#' @export
is_prime <- function(x) {
  if (is.numeric(x) && !all(x %% 1 == 0, na.rm = TRUE)) {
    insight::format_error("`x` needs to be an integer value.")
  }
  vapply(x, function(.x) .x == 2L || all(.x %% 2L:max(2, floor(sqrt(.x))) != 0), logical(1))
}
