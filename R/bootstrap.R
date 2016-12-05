#' @title Generate nonparametric bootstrap replications
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
#' @return A \code{\link[tibble]{tibble}} with one column: a list-variable
#'           \code{strap}, which contains resample-objects of class \code{sj_resample}.
#'           These resample-objects are lists with three elemenents:
#'           \enumerate{
#'             \item the original data frame, \code{data}
#'             \item the rownmumbers \code{id}, i.e. rownumbers of \code{data}, indicating the resampled rows with replacement
#'             \item the \code{resample.id}, indicating the index of the resample (i.e. the position of the \code{sj_resample}-object in the list \code{strap})
#'           }
#'
#' @details By default, each bootstrap sample has the same number of observations
#'            as \code{data}. To generate bootstrap samples without resampling
#'            same observations (i.e. sampling without replacement), use
#'            \code{size} to get bootstrapped data with a specific number
#'            of observations. However, specifying the \code{size}-argument is much
#'            less memory-efficient than the bootstrap with replacement. Hence,
#'            it is recommended to ignore the \code{size}-argument, if it is
#'            not really needed.
#'
#' @note This function applies nonparametric bootstrapping, i.e. the function
#'       draws samples with replacement.
#'       \cr \cr
#'       There is an \code{as.data.frame}- and a \code{print}-method to get or
#'       print the resampled data frames. See 'Examples'. The \code{as.data.frame}-
#'       method automatically applies whenever coercion is done because a data
#'       frame is required as input. See 'Examples' in \code{\link{boot_ci}}.
#'
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
#' # compute standard error for a simple vector from bootstraps
#' # use the `as.data.frame()`-method to get the resampled
#' # data frame
#' bs <- bootstrap(efc, 100)
#' bs$c12hour <- unlist(lapply(bs$strap, function(x) {
#'   mean(as.data.frame(x)$c12hour, na.rm = TRUE)
#' }))
#' # bootstrapped standard error
#' boot_se(bs, c12hour)
#' # standard error of original variable
#' se(efc$c12hour)
#'
#' @importFrom tibble tibble
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
  # generate bootstrap resamples
  strap <- replicate(n, resample(data, size, repl), simplify = F)
  # add resample ID, may be used for other functions (like 'se()' for 'icc()')
  for (i in seq_len(length(strap))) strap[[i]]$resample.id <- i
  # return tibble
  tibble::tibble(strap)
}


resample <- function(data, size, replace) {
  structure(class = "sj_resample", list(data = data, id = sample(nrow(data), size = size, replace = replace)))
}
