#' @title Get p-values for merMod objects
#' @name merMod_p
#'
#' @description This function computes p-values for mixed effects models
#'                (\code{merMod}-objects) that have been fitted with the
#'                \pkg{lme4}-package.
#'
#' @param x A fitted (generalized) linear (mixed) model (\code{merMod}-object).
#' @param p.kr Logical, if \code{TRUE}, the computation of p-values is based on
#'         conditional F-tests with Kenward-Roger approximation for the df (see
#'         'Details').
#'
#' @return A named vector with p-values for the model coefficients.
#'
#' @details For linear mixed models (\code{lmerMod}-objects), the computation of
#'         p-values (if \code{p.kr = TRUE}) is based on conditional F-tests
#'         with Kenward-Roger approximation for the df, using the
#'         \pkg{pbkrtest}-package. If \pkg{pbkrtest} is not available or
#'         \code{p.kr = FALSE}, or if \code{x} is a \code{glmerMod}-object,
#'         computation of p-values is based on normal-distribution assumption,
#'         treating the t-statistics as Wald z-statistics.
#'         \cr \cr
#'         If p-values already have been computed (e.g. for \code{merModLmerTest})-objects
#'         from the \pkg{lmerTest}-package), these will be returned.
#'
#' @examples
#' library(lme4)
#' # fit model
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' merMod_p(fit)
#' merMod_p(fit, p.kr = FALSE)
#'
#' @importFrom stats coef pt pnorm
#' @export
merMod_p <- function(x, p.kr = TRUE) {
  # retrieve sigificance level of independent variables (p-values)
  if (any(class(x) == "merModLmerTest") && requireNamespace("lmerTest", quietly = TRUE)) {
    cs <- suppressWarnings(stats::coef(lmerTest::summary(x)))
  } else {
    cs <- stats::coef(summary(x))
  }
  # remeber coef-names
  coef_names <- rownames(cs)
  # check if we have p-values in summary
  if (ncol(cs) >= 4) {
    # do we have a p-value column?
    pvcn <- which(colnames(cs) == "Pr(>|t|)")
    # if not, default to 4
    if (length(pvcn) == 0) pvcn <- 4
    pv <- cs[, pvcn]
  } else if (any(class(x) == "lmerMod") && requireNamespace("pbkrtest", quietly = TRUE) && p.kr) {
    # compute Kenward-Roger-DF for p-statistic. Code snippet adapted from
    # http://mindingthebrain.blogspot.de/2014/02/three-ways-to-get-parameter-specific-p.html
    message("Computing p-values via Kenward-Roger approximation. Use `p.kr = FALSE` if computation takes too long.")
    #first coefficients need to be data frame
    cs <- as.data.frame(cs)
    # get KR DF
    df.kr <- suppressMessages(pbkrtest::get_Lb_ddf(x, lme4::fixef(x)))
    # compute p-values, assuming an approximate t-dist
    pv <- 2 * stats::pt(abs(cs$`t value`), df = df.kr, lower.tail = FALSE)
    # name vector
    names(pv) <- coef_names
  } else {
    message("Computing p-values via Wald-statistics approximation (treating t as Wald z).")
    pv <- 2 * stats::pnorm(abs(cs[, 3]), lower.tail = FALSE)
  }
  return(pv)
}
