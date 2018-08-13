#' @title Determining distribution parameters
#' @name find_beta
#'
#' @description \code{find_beta()}, \code{find_normal()} and \code{find_cauchy()} find the
#'              shape, mean and standard deviation resp. the location and scale parameters
#'              to describe the beta, normal or cauchy distribution, based on two
#'              percentiles. \code{find_beta2()} finds the shape parameters for a Beta
#'              distribution, based on a probability value and its standard error
#'              or confidence intervals.
#'
#' @param x1 Value for the first percentile.
#' @param p1 Probability of the first percentile.
#' @param x2 Value for the second percentile.
#' @param p2 Probability of the second percentile.
#' @param x Numeric, a probability value between 0 and 1. Typically indicates
#'          a prevalence rate of an outcome of interest; Or an integer value
#'          with the number of observed events. In this case, specify \code{n}
#'          to indicate the toral number of observations.
#' @param se The standard error of \code{x}. Either \code{se} or \code{ci} must
#'          be specified.
#' @param ci The upper limit of the confidence interval of \code{x}. Either
#'          \code{se} or \code{ci} must be specified.
#' @param n Numeric, number of total observations. Needs to be specified, if
#'          \code{x} is an integer (number of observed events), and no
#'          probability. See 'Examples'.
#'
#' @return A list of length two, with the two distribution parameters than can
#'         be used to define the distribution, which (best) describes
#'         the shape for the given input parameters.
#'
#' @details These functions can be used to find parameter for various distributions,
#'          to define prior probabilities for Bayesian analyses. \code{x1},
#'          \code{p1}, \code{x2} and \code{p2} are parameters that describe two
#'          quantiles. Given this knowledge, the distribution parameters are
#'          returned. \cr \cr
#'          Use \code{find_beta2()}, if the known parameters are, e.g. a prevalence
#'          rate or similar probability, and its standard deviation or confidence
#'          interval. In this case. \code{x} should be a probability,
#'          for example a prevalence rate of a certain event. \code{se} then
#'          needs to be the standard error for this probability. Alternatively,
#'          \code{ci} can be specified, which should indicate the upper limit
#'          of the confidence interval od the probability (prevalence rate) \code{x}.
#'          If the number of events out of a total number of trials is known
#'          (e.g. 12 heads out of 30 coin tosses), \code{x} can also be the number
#'          of observed events, while \code{n} indicates the total amount of trials
#'          (in the above example, the function call would be: \code{find_beta2(x = 12, n = 30)}).
#'
#' @references Cook JD. Determining distribution parameters from quantiles. 2010: Department of Biostatistics, Texas (\href{https://www.johndcook.com/quantiles_parameters.pdf}{PDF})
#'
#' @examples
#' # example from blogpost:
#' # https://www.johndcook.com/blog/2010/01/31/parameters-from-percentiles/
#' # 10% of patients respond within 30 days of treatment
#' # and 80% respond within 90 days of treatment
#' find_normal(x1 = 30, p1 = .1, x2 = 90, p2 = .8)
#' find_cauchy(x1 = 30, p1 = .1, x2 = 90, p2 = .8)
#'
#' parms <- find_normal(x1 = 30, p1 = .1, x2 = 90, p2 = .8)
#' curve(
#'   dnorm(x, mean = parms$mean, sd = parms$sd),
#'   from = 0, to = 200
#' )
#'
#' parms <- find_cauchy(x1 = 30, p1 = .1, x2 = 90, p2 = .8)
#' curve(
#'   dcauchy(x, location = parms$location, scale = parms$scale),
#'   from = 0, to = 200
#' )
#'
#'
#' find_beta2(x = .25, ci = .5)
#'
#' shapes <- find_beta2(x = .25, ci = .5)
#' curve(dbeta(x, shapes[[1]], shapes[[2]]))
#'
#' # find Beta distribution for 3 events out of 20 observations
#' find_beta2(x = 3, n = 20)
#'
#' shapes <- find_beta2(x = 3, n = 20)
#' curve(dbeta(x, shapes[[1]], shapes[[2]]))
#'
#' @importFrom stats pbeta approx
#' @importFrom purrr map_dbl
#' @export
find_beta <- function(x1, p1, x2, p2) {
  logK <- seq(-5, 10, length = 200)
  K <- exp(logK)

  m <- purrr::map_dbl(K, ~ betaprior(.x, x1, p1))

  prob2 <- stats::pbeta(x2, K * m, K * (1 - m))
  ind <- ((prob2 > 0) & (prob2 < 1))
  app <- stats::approx(prob2[ind], logK[ind], p2)
  K0 <- exp(app$y)
  m0 <- betaprior(K0, x1, p1)

  s1 <- K0 * m0
  s2 <- K0 * (1 - m0)

  list(shape1 = s1, shape2 = s2)
}


betaprior <- function(K, x, p) {
  m.lo <- 0
  m.hi <- 1
  flag <- TRUE

  while (flag) {
    m0 <- (m.lo + m.hi) / 2
    p0 <- stats::pbeta(x, K * m0, K * (1 - m0))

    if (p0 < p)
      m.hi <- m0
    else
      m.lo <- m0

    if (abs(p0 - p) < 1e-04) flag <- FALSE
  }

  m0
}


#' @rdname find_beta
#' @export
find_beta2 <- function(x, se, ci, n) {
  # check if all required arguments are given
  if (missing(se) && missing(ci) && missing(n)) {
    stop("Either `se` or `ci`, or `n` must be specified.", call. = F)
  }

  # for number of observations, compute variance of beta distribution
  if (!missing(n)) {
    if (!is.integer(x) && x < 1)
      stop("If `n` is given, x` must be an integer value greater than 0.", call. = F)

    # compute 2 SD from beta variance
    bvar <- 2 * sqrt((x * n) / ((x + n)^2 * (x + n + 1)))

    # need to compute proportion
    x <- x / n
    p2 <- .95
    x2 <- x + bvar
  }

  # for standard errors, we assume a 68% quantile
  if (!missing(se)) {
    p2 <- .68
    x2 <- x + se
  }

  # for CI, we assume a 68% quantile
  if (!missing(ci)) {
    p2 <- .95
    x2 <- ci
  }

  # the probability is assumed to be the median
  p1 <- .5
  x1 <- x

  find_beta(x1, p1, x2, p2)
}


#' @importFrom stats qcauchy
#' @rdname find_beta
#' @export
find_cauchy <- function(x1, p1, x2, p2) {
  # find location paramater
  l <- (x1 * stats::qcauchy(p2) ^ -1 - x2 * stats::qcauchy(p1) ^ -1) / (stats::qcauchy(p2) ^ -1 - stats::qcauchy(p1) ^ -1)
  s <- (x2 - x1) / (stats::qcauchy(p2) ^ -1 - stats::qcauchy(p1) ^ -1)

  list(location = l, scale = s)
}



#' @importFrom stats qnorm
#' @rdname find_beta
#' @export
find_normal <- function(x1, p1, x2, p2) {
  # find location paramater
  mw <- (x1 * stats::qnorm(p2) ^ -1 - x2 * stats::qnorm(p1) ^ -1) / (stats::qnorm(p2) ^ -1 - stats::qnorm(p1) ^ -1)
  stddev <- (x2 - x1) / (stats::qnorm(p2) ^ -1 - stats::qnorm(p1) ^ -1)

  list(mean = mw, sd = stddev)
}

