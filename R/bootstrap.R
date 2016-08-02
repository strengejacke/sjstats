#' @title Generate bootstrap replications
#' @name bootstrap
#'
#' @description Generates \code{n} bootstrap samples of \code{data} and
#'              returns the bootstrapped data frames as list-variable.
#'
#' @param data A data frame.
#' @param n Number of bootstraps to be generated
#' @param size Optional, size of the bootstrap samples. May either be a number
#'          between 1 and \code{nrow(data)} or a value between 0 and 1 to sample
#'          a proportion of observations from \code{data} (see 'Examples').
#'
#' @return A \code{\link[tibble]{data_frame}} with one column: a list-variable
#'           \code{strap}, which contains the bootstrapped samples from \code{data}.
#'
#' @details By default, each bootstrap sample has the same number of observations
#'            as \code{data}. To generate bootstrap samples without resampling
#'            same observations (i.e. sampling without replacement), use
#'            \code{size} to get bootstrapped data with a specific number
#'            of observations.
#'
#' @seealso \code{\link{boot_ci}} to calculate confidence intervals from
#'            bootstrap samples.
#'
#' @examples
#' data(efc)
#' bs <- bootstrap(efc, 5)
#'
#' # now run models for each bootstrapped sample
#' lapply(bs$strap, function(x) lm(neg_c_7 ~ e42dep + c161sex, data = x))
#'
#' # generate bootstrap samples with 600 observations for each sample
#' bs <- bootstrap(efc, 5, 600)
#'
#' # generate bootstrap samples with 70% observations of the original sample size
#' bs <- bootstrap(efc, 5, .7)
#'
#'
#' @importFrom tibble data_frame
#' @export
bootstrap <- function(data, n, size) {
  if (!missing(size) && !is.null(size)) {
    # check for valid range
    if (size < 0 || size > nrow(data))
      stop("`size` must be greater than 0, but not greater than number of rows of `data`.", call. = F)
    # check if we want proportions
    if (size < 1) size <- as.integer(nrow(data) * size)
    # generate bootstraps w/o replacement
    repl <- F
  } else {
    # size = observations
    size <- nrow(data)
    # generate bootstraps with replacement
    repl <- T
  }
  tibble::data_frame(
    strap = replicate(n, resample(data, size, repl), simplify = F)
  )
}


resample <- function(data, size, replace) {
  structure(class = "sj_resample", list(data = data, id = sample(nrow(data), size = size, replace = replace)))
}
