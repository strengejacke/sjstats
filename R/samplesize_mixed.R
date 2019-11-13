#' @title Sample size for linear mixed models
#' @name samplesize_mixed
#'
#' @description Compute an approximated sample size for linear mixed models
#'                (two-level-designs), based on power-calculation for standard
#'                design and adjusted for design effect for 2-level-designs.
#'
#' @param eff.size Effect size.
#' @param df.n Optional argument for the degrees of freedom for numerator. See 'Details'.
#' @param power Power of test (1 minus Type II error probability).
#' @param sig.level Significance level (Type I error probability).
#' @param k Number of cluster groups (level-2-unit) in multilevel-design.
#' @param n Optional, number of observations per cluster groups
#'       (level-2-unit) in multilevel-design.
#' @param icc Expected intraclass correlation coefficient for multilevel-model.
#'
#' @return A list with two values: The number of subjects per cluster, and the
#'           total sample size for the linear mixed model.
#'
#' @references Cohen J. 1988. Statistical power analysis for the behavioral sciences (2nd ed.). Hillsdale,NJ: Lawrence Erlbaum.
#'             \cr \cr
#'             Hsieh FY, Lavori PW, Cohen HJ, Feussner JR. 2003. An Overview of Variance Inflation Factors for Sample-Size Calculation. Evaluation and the Health Professions 26: 239-257. \doi{10.1177/0163278703255230}
#'             \cr \cr
#'             Snijders TAB. 2005. Power and Sample Size in Multilevel Linear Models. In: Everitt BS, Howell DC (Hrsg.). Encyclopedia of Statistics in Behavioral Science. Chichester, UK: John Wiley and Sons, Ltd. \doi{10.1002/0470013192.bsa492}
#'
#' @details The sample size calculation is based on a power-calculation for the
#'          standard design. If \code{df.n} is not specified, a power-calculation
#'          for an unpaired two-sample t-test will be computed (using
#'          \code{\link[pwr]{pwr.t.test}} of the \CRANpkg{pwr}-package).
#'          If \code{df.n} is given, a power-calculation for general linear models
#'          will be computed (using \code{\link[pwr]{pwr.f2.test}} of the
#'          \pkg{pwr}-package). The sample size of the standard design
#'          is then adjusted for the design effect of two-level-designs (see
#'          \code{\link{design_effect}}). Thus, the sample size calculation is appropriate
#'          in particular for two-level-designs (see \cite{Snijders 2005}). Models that
#'          additionally include repeated measures (three-level-designs) may work
#'          as well, however, the computed sample size may be less accurate.
#'
#' @examples
#' # Sample size for multilevel model with 30 cluster groups and a small to
#' # medium effect size (Cohen's d) of 0.3. 27 subjects per cluster and
#' # hence a total sample size of about 802 observations is needed.
#' samplesize_mixed(eff.size = .3, k = 30)
#'
#' # Sample size for multilevel model with 20 cluster groups and a medium
#' # to large effect size for linear models of 0.2. Five subjects per cluster and
#' # hence a total sample size of about 107 observations is needed.
#' samplesize_mixed(eff.size = .2, df.n = 5, k = 20, power = .9)
#'
#'
#' @export
samplesize_mixed <- function(eff.size, df.n = NULL, power = .8, sig.level = .05, k, n, icc = 0.05) {
  if (!requireNamespace("pwr", quietly = TRUE)) {
    stop("Package `pwr` needed for this function to work. Please install it.", call. = FALSE)
  }

  # compute sample size for standard design
  if (is.null(df.n))
    # if we have no degrees of freedom specified, use t-test
    obs <- 2 * pwr::pwr.t.test(d = eff.size, sig.level = sig.level, power = power)$n
  else
    # we have df, so power-calc for linear models
    obs <- pwr::pwr.f2.test(u = df.n, f2 = eff.size, sig.level = sig.level, power = power)$v + df.n + 1

  # if we have no information on the number of observations per cluster,
  # compute this number now
  if (missing(n) || is.null(n)) {
    n <- (obs * (1 - icc)) / (k - (obs * icc))
    if (n < 1) {
      warning("Minimum required number of subjects per cluster is negative and was adjusted to be positive. You may reduce the requirements for the multi-level structure (i.e. reduce `k` or `icc`), or you can increase the effect-size.", call. = FALSE)
      n <- 1
    }
  }

  # adjust standard design by design effect
  total.n <- obs * design_effect(n = n, icc = icc)


  # sample size for each group and total n
  smpsz <- list(round(total.n / k), round(total.n))

  # name list
  names(smpsz) <- c("Subjects per Cluster", "Total Sample Size")
  smpsz
}


#' @rdname samplesize_mixed
#' @export
smpsize_lmm <- samplesize_mixed
