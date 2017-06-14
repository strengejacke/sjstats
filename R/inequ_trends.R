#' @title Compute trends in status inequalities
#' @name inequ_trend
#'
#' @description This method computes the proportional change of absolute
#'              (rate differences) and relative (rate ratios) inequalities
#'              of prevalence rates for two different status groups, as proposed
#'              by Mackenbach et al. (2015).
#'
#' @param data A data frame that contains the variables with prevalence rates for both low
#'          and high status groups (see 'Examples').
#' @param prev.low The name of the variable with the prevalence rates for
#'          the low status groups.
#' @param prev.hi The name of the variable with the prevalence rates for
#'          the hi status groups.
#'
#' @return A data frame with the prevalence rates as well as the values for the
#'           proportional change in absolute (\code{rd}) and relative (\code{rr})
#'           ineqqualities.
#'
#' @references Mackenbach JP, Martikainen P, Menvielle G, de Gelder R. 2015. The Arithmetic of Reducing Relative and Absolute Inequalities in Health: A Theoretical Analysis Illustrated with European Mortality Data. Journal of Epidemiology and Community Health 70(7): 730–36. \doi{10.1136/jech-2015-207018}
#'
#' @details Given the time trend of prevalence rates of an outcome for two status
#'          groups (e.g. the mortality rates for people with lower and higher
#'          socioeconomic status over 40 years), this function computes the
#'          proportional change of absolute and relative inequalities, expressed
#'          in changes in rate differences and rate ratios. The function implements
#'          the algorithm proposed by \emph{Mackenbach et al. 2015}.
#'
#' @examples
#' # This example reproduces Fig. 1 of Mackenbach et al. 2015, p.5
#'
#' # 40 simulated time points, with an initial rate ratio of 2 and
#' # a rate difference of 100 (i.e. low status group starts with a
#' # prevalence rate of 200, the high status group with 100)
#'
#' # annual decline of prevalence is 1% for the low, and 3% for the
#' # high status group
#'
#' n <- 40
#' time <- seq(1, n, by = 1)
#' lo <- rep(200, times = n)
#' for (i in 2:n) lo[i] <- lo[i - 1] * .99
#'
#' hi <- rep(100, times = n)
#' for (i in 2:n) hi[i] <- hi[i - 1] * .97
#'
#' prev.data <- data.frame(lo, hi)
#'
#' # print values
#' inequ_trend(prev.data, lo, hi)
#'
#' # plot trends - here we see that the relative inequalities
#' # are increasing over time, while the absolute inequalities
#' # are first increasing as well, but later are decreasing
#' # (while rel. inequ. are still increasing)
#' plot(inequ_trend(prev.data, lo, hi))
#'
#' @importFrom dplyr select
#' @export
inequ_trend <- function(data, prev.low, prev.hi) {
  # prepare data for prevalence rates for low and hi status groups
  if (is.null(data) || missing(data)) {
    dat <- data.frame(prev.low, prev.hi)
  } else {
    # get variable names
    low <- deparse(substitute(prev.low))
    high <- deparse(substitute(prev.hi))
    dat <- dplyr::select_(data, low, high)

  }

  # ensure common column names
  colnames(dat) <- c("lo", "hi")

  # trends in rate ratios

  # compute relative inequality for first time point, needed
  # as reference to compute proportional change over time
  dat$rr <- dat$lo[1] / dat$hi[1]

  # compute proportional change of relative inequalities over time
  for (t in 2:nrow(dat)) {
    delta.low <- (dat$lo[t] - dat$lo[t - 1]) / dat$lo[t - 1]
    delta.hi <- (dat$hi[t] - dat$hi[t - 1]) / dat$hi[t - 1]
    dat$rr[t] <- dat$rr[t - 1] * ((1 + delta.low) / (1 + delta.hi))
  }

  # trends in rate difference

  # compute absolute inequality for first time point, needed
  # as reference to compute proportional change over time
  dat$rd <- dat$lo[1] - dat$hi[1]

  # compute proportional change of absolute inequalities over time
  for (t in 2:nrow(dat)) {
    delta.low <- (dat$lo[t] - dat$lo[t - 1]) / dat$lo[t - 1]
    delta.hi <- (dat$hi[t] - dat$hi[t - 1]) / dat$hi[t - 1]
    dat$rd[t] <- dat$rd[t - 1] + (dat$lo[t - 1 ] * delta.low - dat$hi[t - 1] * delta.hi)
  }

  # return
  structure(class = "sj_inequ_trend", list(data = dat))
}
