#' @title Find prime numbers
#' @name is_prime
#'
#' @description This functions checks whether a number is, or numbers in a
#'    vector are prime numbers.
#'
#' @param x An integer, or a vector of integers.
#'
#' @return \code{TRUE} for each prime number in \code{x}, \code{FALSE} otherwise.
#'
#' @examples
#' is_prime(89)
#' is_prime(15)
#' is_prime(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
#'
#' @importFrom purrr map_lgl
#' @importFrom sjmisc is_float
#' @export
is_prime <- function(x) {
  if (sjmisc::is_float(x))
    stop("`x` needs to be an integer value.", call. = F)

  purrr::map_lgl(x, ~ .x == 2L || all(.x %% 2L:max(2, floor(sqrt(.x))) != 0))
}
