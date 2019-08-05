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
#' @importFrom lme4 ngrps
#' @importFrom stats nobs
#' @importFrom purrr map_dbl
#' @export
se_ybar <- function(fit) {
  # get model icc
  vars <- insight::get_variance(fit, verbose = FALSE)

  # get group variances
  tau.00 <- unname(vars$var.intercept)

  # total variance
  tot_var <- sum(tau.00, vars$var.residual)

  # get number of groups
  m.cnt <- lme4::ngrps(fit)

  # compute number of observations per group (level-2-unit)
  obs <- round(stats::nobs(fit) / m.cnt)

  # compute simple icc
  icc <- tau.00 / tot_var

  # compute standard error of sample mean
  se <- purrr::map_dbl(seq_len(length(m.cnt)), ~ sqrt((tot_var / stats::nobs(fit)) * design_effect(n = obs[.x], icc = icc[.x])))

  # give names for se, so user sees, which random effect has what impact
  names(se) <- names(m.cnt)
  se
}
