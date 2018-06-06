#' @title Get p-values from regression model objects
#' @name p_value
#'
#' @description This function returns the p-values for fitted model objects.
#'
#' @param fit A fitted model object of class \code{lm}, \code{glm}, \code{merMod},
#'        \code{merModLmerTest}, \code{pggls} or \code{gls}. Other classes may
#'        work as well.
#' @param p.kr Logical, if \code{TRUE}, the computation of p-values is based on
#'         conditional F-tests with Kenward-Roger approximation for the df (see
#'         'Details').
#'
#' @return A \code{data.frame} with the model coefficients' names (\code{term}),
#'         p-values (\code{p.value}) and standard errors (\code{std.error}).
#'
#' @details For linear mixed models (\code{lmerMod}-objects), the computation of
#'         p-values (if \code{p.kr = TRUE}) is based on conditional F-tests
#'         with Kenward-Roger approximation for the df, using the
#'         \CRANpkg{pbkrtest}-package. If \pkg{pbkrtest} is not available or
#'         \code{p.kr = FALSE}, or if \code{x} is a \code{glmerMod}-object,
#'         computation of p-values is based on normal-distribution assumption,
#'         treating the t-statistics as Wald z-statistics.
#'         \cr \cr
#'         If p-values already have been computed (e.g. for \code{merModLmerTest}-objects
#'         from the \CRANpkg{lmerTest}-package), these will be returned.
#'         \cr \cr
#'         The \code{print()}-method has a \code{summary}-argument, that - in
#'         case \code{p.kr = TRUE} - also prints information on the approximated
#'         degrees of freedom (see 'Examples').
#'
#' @examples
#' data(efc)
#' # linear model fit
#' fit <- lm(neg_c_7 ~ e42dep + c172code, data = efc)
#' p_value(fit)
#'
#' # Generalized Least Squares fit
#' library(nlme)
#' fit <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
#'            correlation = corAR1(form = ~ 1 | Mare))
#' p_value(fit)
#'
#' # lme4-fit
#' library(lme4)
#' sleepstudy$mygrp <- sample(1:45, size = 180, replace = TRUE)
#' fit <- lmer(Reaction ~ Days + (1 | mygrp) + (1 | Subject), sleepstudy)
#' pv <- p_value(fit, p.kr = TRUE)
#'
#' # normal output
#' pv
#'
#' # add information on df and t-statistic
#' print(pv, summary = TRUE)
#'
#' @importFrom stats coef pnorm
#' @importFrom broom tidy
#' @importFrom dplyr select
#' @export
p_value <- function(fit, p.kr = FALSE) {
  # retrieve sigificance level of independent variables (p-values)
  if (inherits(fit, "pggls")) {
    p <- summary(fit)$CoefTable[, 4]
    se <- summary(fit)$CoefTable[, 2]
  } else if (inherits(fit, "gls")) {
    p <- summary(fit)$tTable[, 4]
    se <- summary(fit)$tTable[, 2]
  } else if (inherits(fit, "Zelig-relogit")) {
    if (!requireNamespace("Zelig", quietly = T))
      stop("Package `Zelig` required. Please install", call. = F)
    p <- unlist(Zelig::get_pvalue(fit))
    se <- unlist(Zelig::get_se(fit))
  } else if (is_merMod(fit)) {
    p <- merMod_p(fit, p.kr)
    sekr <- attr(p, "se.kr", exact = TRUE)
    if (!is.null(sekr))
      se <- sekr
    else
      se <- stats::coef(summary(fit))[, 2]
  } else if (inherits(fit, c("pglm", "maxLik"))) {
    p <- summary(fit)$estimate[, 4]
    se <- summary(fit)$estimate[, 2]
  } else if (inherits(fit, "gam")) {
    sm <- summary(fit)
    lc <- length(sm$p.coeff)
    p <- sm$p.pv[1:lc]
    se <- sm$se[1:lc]
  } else if (inherits(fit, "polr")) {
    smry <- suppressMessages(as.data.frame(stats::coef(summary(fit))))
    tstat <- smry[[3]]
    se <- smry[[2]]
    p <- 2 * stats::pnorm(abs(tstat), lower.tail = FALSE)
    names(p) <- rownames(smry)
  } else if (inherits(fit, "multinom")) {
    return(fit %>%
      broom::tidy() %>%
      dplyr::select(.data$term, .data$p.value, .data$std.error))
  } else {
    p <- stats::coef(summary(fit))[, 4]
    se <- stats::coef(summary(fit))[, 2]
  }

  res <- data.frame(
    term = names(p),
    p.value = as.vector(p),
    std.error = as.vector(se)
  )

  class(res) <- c("sj_pval", class(res))

  attr(res, "df.kr") <- attr(p, "df.kr", exact = TRUE)
  attr(res, "se.kr") <- attr(p, "se.kr", exact = TRUE)
  attr(res, "t.kr") <- attr(p, "t.kr", exact = TRUE)

  res
}




#' @importFrom stats coef pt pnorm
merMod_p <- function(fit, p.kr) {
  # retrieve sigificance level of independent variables (p-values)
  cs <- stats::coef(summary(fit))

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
    df.kr <- get_kr_df(fit)
    se.kr <- get_kr_se(fit)
    t.kr <- lme4::fixef(fit) / se.kr
    # compute p-values, assuming an approximate t-dist
    pv <- 2 * stats::pt(abs(t.kr), df = df.kr, lower.tail = FALSE)
    # name vector
    names(pv) <- coef_names
    attr(pv, "df.kr") <- df.kr
    attr(pv, "se.kr") <- se.kr
    attr(pv, "t.kr") <- t.kr
  } else {
    message("Computing p-values via Wald-statistics approximation (treating t as Wald z).")
    pv <- 2 * stats::pnorm(abs(cs[, 3]), lower.tail = FALSE)
  }


  if (is.null(names(pv))) {
    if (length(coef_names) == length(pv)) names(pv) <- coef_names
  }

  pv
}


#' @importFrom lme4 fixef
get_kr_df <- function(x) {
  if (!requireNamespace("pbkrtest", quietly = TRUE))
    stop("Package `pbkrtest` required for Kenward-Rogers approximation.", call. = FALSE)

  L <- diag(rep(1, length(lme4::fixef(x))))
  L <- as.data.frame(L)
  out <- suppressMessages(purrr::map_dbl(L, pbkrtest::get_Lb_ddf, object = x))
  names(out) <- names(lme4::fixef(x))
  out
}


#' @importFrom lme4 fixef
#' @importFrom purrr map_dbl
get_kr_se <- function(x) {
  if (!requireNamespace("pbkrtest", quietly = TRUE))
    stop("Package `pbkrtest` required for Kenward-Roger approximation.", call. = FALSE)

  vcov_adj <- pbkrtest::vcovAdj(x)

  fe <- lme4::fixef(x)
  le <- length(fe)
  Lmat <- diag(le)

  se <- purrr::map_dbl(1:le, ~ sqrt(qform(Lmat[.x, ], as.matrix(vcov_adj))))
  names(se) <- names(fe)

  se
}


qform <- function(x, A) sum(x * (A %*% x))
