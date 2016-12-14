#' @title Get p-values from regression model objects
#' @name get_model_pval
#'
#' @description This function returns the p-values for fitted model objects.
#'
#' @param fit A fitted model object of \code{lm}, \code{glm}, \code{merMod},
#'        \code{merModLmerTest}, \code{pggls} or \code{gls}. Other classes may
#'        work as well.
#' @param p.kr Logical, if \code{TRUE}, the computation of p-values is based on
#'         conditional F-tests with Kenward-Roger approximation for the df (see
#'         'Details').
#'
#' @return A \code{\link{tibble}} with the model coefficients' names (\code{term}),
#'         p-values (\code{p.value}) and standard errors (\code{std.error}).
#'
#' @note \code{merMod_p()} will become deprecated soon, please use
#'       \code{get_model_pval()} instead.
#'
#' @details For linear mixed models (\code{lmerMod}-objects), the computation of
#'         p-values (if \code{p.kr = TRUE}) is based on conditional F-tests
#'         with Kenward-Roger approximation for the df, using the
#'         \CRANpkg{pbkrtest}-package. If \pkg{pbkrtest} is not available or
#'         \code{p.kr = FALSE}, or if \code{x} is a \code{glmerMod}-object,
#'         computation of p-values is based on normal-distribution assumption,
#'         treating the t-statistics as Wald z-statistics.
#'         \cr \cr
#'         If p-values already have been computed (e.g. for \code{merModLmerTest})-objects
#'         from the \CRANpkg{lmerTest}-package), these will be returned.
#'
#' @examples
#' data(efc)
#' # linear model fit
#' fit <- lm(neg_c_7 ~ e42dep + c172code, data = efc)
#' get_model_pval(fit)
#'
#' # Generalized Least Squares fit
#' library(nlme)
#' fit <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
#'            correlation = corAR1(form = ~ 1 | Mare))
#' get_model_pval(fit)
#'
#' # lme4-fit
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' get_model_pval(fit, p.kr = TRUE)
#'
#' @importFrom stats coef
#' @importFrom tibble tibble
#' @export
get_model_pval <- function(fit, p.kr = FALSE) {
  # retrieve sigificance level of independent variables (p-values)
  if (inherits(fit, "pggls")) {
    p <- summary(fit)$CoefTable[, 4]
    se <- summary(fit)$CoefTable[, 2]
  } else if (inherits(fit, "gls")) {
    p <- summary(fit)$tTable[, 4]
    se <- summary(fit)$tTable[, 2]
  } else if (is_merMod(fit)) {
    p <- merMod_p(fit, p.kr)
    se <- stats::coef(summary(fit))[, 2]
  } else if (inherits(fit, c("pglm", "maxLik"))) {
    # p-values
    p <- summary(fit)$estimate[, 4]
    se <- summary(fit)$estimate[, 2]
  } else {
    p <- stats::coef(summary(fit))[, 4]
    se <- stats::coef(summary(fit))[, 2]
  }
  tibble::tibble(term = names(p), p.value = as.vector(p), std.error = as.vector(se))
}


#' @rdname get_model_pval
#' @importFrom stats coef pt pnorm
#' @export
merMod_p <- function(fit, p.kr = TRUE) {
  # retrieve sigificance level of independent variables (p-values)
  if (inherits(fit, "merModLmerTest") && requireNamespace("lmerTest", quietly = TRUE)) {
    cs <- suppressWarnings(stats::coef(lmerTest::summary(fit)))
  } else {
    cs <- stats::coef(summary(fit))
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
  } else if (inherits(fit, "lmerMod") && requireNamespace("pbkrtest", quietly = TRUE) && p.kr) {
    # compute Kenward-Roger-DF for p-statistic. Code snippet adapted from
    # http://mindingthebrain.blogspot.de/2014/02/three-ways-to-get-parameter-specific-p.html
    message("Computing p-values via Kenward-Roger approximation. Use `p.kr = FALSE` if computation takes too long.")
    #first coefficients need to be data frame
    cs <- as.data.frame(cs)
    # get KR DF
    df.kr <- suppressMessages(pbkrtest::get_Lb_ddf(fit, lme4::fixef(fit)))
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
