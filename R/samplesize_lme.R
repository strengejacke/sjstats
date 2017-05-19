#' @title Sample size for linear mixed models
#' @name smpsize_lmm
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
#' @param icc Expected intraclass correlation coefficient for multilevel-model.
#'
#' @return A list with two values: The number of subjects per cluster, and the
#'           total sample size for the linear mixed model.
#'
#' @references Cohen J. 1988. Statistical power analysis for the behavioral sciences (2nd ed.). Hillsdale,NJ: Lawrence Erlbaum.
#'             \cr \cr
#'             Hsieh FY, Lavori PW, Cohen HJ, Feussner JR. 2003. An Overview of Variance Inflation Factors for Sample-Size Calculation. Evaluation & the Health Professions 26: 239–257. \doi{10.1177/0163278703255230}
#'             \cr \cr
#'             Snijders TAB. 2005. Power and Sample Size in Multilevel Linear Models. In: Everitt BS, Howell DC (Hrsg.). Encyclopedia of Statistics in Behavioral Science. Chichester, UK: John Wiley & Sons, Ltd. \doi{10.1002/0470013192.bsa492}
#'
#' @details The sample size calculation is based on a power-calculation for the
#'          standard design. If \code{df.n} is not specified, a power-calculation
#'          for an unpaired two-sample t-test will be computed (using
#'          \code{\link[pwr]{pwr.t.test}} of the \CRANpkg{pwr}-package).
#'          If \code{df.n} is given, a power-calculation for general linear models
#'          will be computed (using \code{\link[pwr]{pwr.f2.test}} of the
#'          \pkg{pwr}-package). The sample size of the standard design
#'          is then adjusted for the design effect of two-level-designs (see
#'          \code{\link{deff}}). Thus, the sample size calculation is appropriate
#'          in particular for two-level-designs (see \cite{Snijders 2005}). Models that
#'          additionally include repeated measures (three-level-designs) may work
#'          as well, however, the computed sample size may be less accurate.
#'
#' @examples
#' # Sample size for multilevel model with 30 cluster groups and a small to
#' # medium effect size (Cohen's d) of 0.3. 29 subjects per cluster and
#' # hence a total sample size of about 859 observations is needed.
#' smpsize_lmm(eff.size = .3, k = 30)
#'
#' # Sample size for multilevel model with 20 cluster groups and a medium
#' # to large effect size for linear models of 0.2. Nine subjects per cluster and
#' # hence a total sample size of about 172 observations is needed.
#' smpsize_lmm(eff.size = .2, df.n = 5, k = 20, power = .9)
#'
#' @export
smpsize_lmm <- function(eff.size, df.n = NULL, power = .8, sig.level = .05, k, icc = 0.05) {
  if (!requireNamespace("pwr", quietly = TRUE)) {
    stop("Package `pwr` needed for this function to work. Please install it.", call. = FALSE)
  }

  # compute sample size for standard design
  if (is.null(df.n))
    # if we have no degrees of freedom specified, use t-test
    n <- 2 * pwr::pwr.t.test(d = eff.size, sig.level = sig.level, power = power)$n
  else
    # we have df, so power-calc for linear models
    n <- pwr::pwr.f2.test(u = df.n, f2 = eff.size, sig.level = sig.level, power = power)$v + df.n + 1

  # adjust standard design by design effect
  total.n <- n * deff(n = k, icc = icc)

  # sample size for each group and total n
  smpsz <- list(round(total.n / k), round(total.n))

  # name list
  names(smpsz) <- c("Subjects per Cluster", "Total Sample Size")
  smpsz
}


#' @title Design effects for two-level mixed models
#' @name deff
#'
#' @description Compute the design effect (also called \emph{Variance Inflation Factor})
#'              for mixed models with two-level design.
#'
#' @param n Average number of observations per grouping cluster (i.e. level-2 unit).
#' @param icc Assumed intraclass correlation coefficient for multilevel-model.
#'
#' @return The design effect (Variance Inflation Factor) for the two-level model.
#'
#' @references Bland JM. 2000. Sample size in guidelines trials. Fam Pract. (17), 17-20.
#'             \cr \cr
#'             Hsieh FY, Lavori PW, Cohen HJ, Feussner JR. 2003. An Overview of Variance Inflation Factors for Sample-Size Calculation. Evaluation & the Health Professions 26: 239–257. \doi{10.1177/0163278703255230}
#'             \cr \cr
#'             Snijders TAB. 2005. Power and Sample Size in Multilevel Linear Models. In: Everitt BS, Howell DC (Hrsg.). Encyclopedia of Statistics in Behavioral Science. Chichester, UK: John Wiley & Sons, Ltd. \doi{10.1002/0470013192.bsa492}
#'             \cr \cr
#'             Thompson DM, Fernald DH, Mold JW. 2012. Intraclass Correlation Coefficients Typical of Cluster-Randomized Studies: Estimates From the Robert Wood Johnson Prescription for Health Projects. The Annals of Family Medicine;10(3):235–40. \doi{10.1370/afm.1347}
#'
#' @details The formula for the design effect is simply \code{(1 + (n - 1) * icc)}.
#'
#' @examples
#' # Design effect for two-level model with 30 cluster groups
#' # and an assumed intraclass correlation coefficient of 0.05.
#' deff(n = 30)
#'
#' # Design effect for two-level model with 24 cluster groups
#' # and an assumed intraclass correlation coefficient of 0.2.
#' deff(n = 24, icc = 0.2)
#'
#' @export
deff <- function(n, icc = 0.05) {
  return(1 + (n - 1) * icc)
}


#' @title Standard error of sample mean for mixed models
#' @name se_ybar
#'
#' @description Compute the standard error for the sample mean for mixed models,
#'                regarding the extent to which clustering affects the standard errors.
#'                May be used as part of the multilevel power calculation for cluster sampling
#'                (see \cite{Gelman and Hill 2007, 447ff}).
#'
#' @param fit Fitted mixed effects model (\code{\link[lme4]{merMod}}-class).
#'
#' @return The standard error of the sample mean of \code{fit}.
#'
#' @references Gelman A, Hill J. 2007. Data analysis using regression and multilevel/hierarchical models. Cambridge, New York: Cambridge University Press
#'
#' @examples
#' library(lme4)
#' fit <- lmer(Reaction ~ 1 + (1 | Subject), sleepstudy)
#' se_ybar(fit)
#'
#' @export
se_ybar <- function(fit) {
  # get model icc
  icc <- icc(fit)

  # get group variances
  tau.00 <- unname(attr(icc, "tau.00", exact = T))

  # total variance
  tot_var <- sum(tau.00, attr(icc, "sigma_2", exact = T))

  # get number of groups
  m.cnt <- unlist(lapply(getME(fit, "flist"), nlevels))

  # compute standard error of sample mean
  y.se <- c()

  for (i in seq_len(length(m.cnt))) {
    y.se <- c(y.se, sqrt((tot_var / nobs(fit)) * deff(n = m.cnt[i], icc = icc[i])))
  }

  y.se
}
