#' @title Check overdispersion of GL(M)M's
#' @name overdisp
#' @description \code{overdisp()} checks generalized linear (mixed) models for
#'              overdispersion, while \code{zero_count()} checks whether models
#'              from poisson-families are over- or underfitting zero-counts in
#'              the outcome.
#'
#' @param x Fitted GLMM (\code{\link[lme4]{merMod}}-class) or \code{glm} model.
#' @param trafo A specification of the alternative, can be numeric or a
#'          (positive) function or \code{NULL} (the default). See 'Details'
#'          in \code{\link[AER]{dispersiontest}} in package \CRANpkg{AER}. Does not
#'          apply to \code{merMod} objects.
#'
#' @return For \code{overdisp()}, information on the overdispersion test; for
#'         \code{zero_count()}, the amount of predicted and observed zeros in
#'         the outcome, as well as the ratio between these two values.
#'
#' @note For the overdispersion-test, the interpretation of the returned p-value
#'       differs between GLM and GLMM. For GLMs, a p-value < .05 indicates
#'       overdispersion, while for GLMMs, a p-value > .05 indicates overdispersion.
#'       \cr \cr
#'       For \code{zero_count()}, a model that is underfitting zero-counts
#'       indicates a zero-inflation in the data, i.e. it is recommended to
#'       use negative binomial or zero-inflated models then.
#'
#' @details For \code{merMod}- and \code{glmmTMB}-objects, \code{overdisp()} is
#'          based on the code in the \href{http://glmm.wikidot.com/faq}{DRAFT r-sig-mixed-models FAQ},
#'          section \emph{How can I deal with overdispersion in GLMMs?}.
#'          Note that this function only returns an \emph{approximate} estimate
#'          of an overdispersion parameter, and is probably inaccurate for
#'          zero-inflated mixed models (fitted with \code{glmmTMB}).
#'          \cr \cr
#'          For \code{glm}'s, \code{overdisp()} simply wraps the \code{dispersiontest}
#'          from the \pkg{AER}-package.
#'
#' @references \href{http://glmm.wikidot.com/faq}{DRAFT r-sig-mixed-models FAQ}
#'             \cr \cr
#'             Bolker B et al. (2017): \href{http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html}{GLMM FAQ.}
#'
#' @examples
#' library(sjmisc)
#' data(efc)
#'
#' # response has many zero-counts, poisson models
#' # might be overdispersed
#' barplot(table(efc$tot_sc_e))
#'
#' fit <- glm(tot_sc_e ~ neg_c_7 + e42dep + c160age,
#'            data = efc, family = poisson)
#' overdisp(fit)
#' zero_count(fit)
#'
#' library(lme4)
#' efc$e15relat <- to_factor(efc$e15relat)
#' fit <- glmer(tot_sc_e ~ neg_c_7 + e42dep + c160age + (1 | e15relat),
#'              data = efc, family = poisson)
#' overdisp(fit)
#' zero_count(fit)
#'
#'
#' @importFrom stats df.residual residuals pchisq
#' @export
overdisp <- function(x, trafo = NULL) {
  if (inherits(x, c("merMod", "glmerMod", "glmmTMB")))
    return(overdisp.lme4(x))
  else
    return(overdisp.default(x, trafo))
}


overdisp.default <- function(x, trafo) {
  # check if suggested package is available
  if (!requireNamespace("AER", quietly = TRUE)) {
    stop("Package `AER` needed for this function to work. Please install it.", call. = FALSE)
  }

  result <- AER::dispersiontest(x, trafo = trafo, alternative = "greater")
  print(result)

  if (result$p.value < 0.05)
    message("Overdispersion detected.")
  else
    message("No overdispersion detected.")

  invisible(result)
}


overdisp.lme4 <- function(x) {
  rdf <- stats::df.residual(x)
  rp <- stats::residuals(x, type = "pearson")
  Pearson.chisq <- sum(rp ^ 2)
  prat <- Pearson.chisq / rdf
  pval <- stats::pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)

  structure(class = "sjstats_ovderdisp",
            list(
              chisq = Pearson.chisq,
              ratio = prat,
              rdf = rdf,
              p = pval
            ))
}


#' @rdname overdisp
#' @importFrom stats predict dpois family
#' @export
zero_count <- function(x) {
  # check if we have poisson
  if (!stats::family(x)$family %in% c("poisson", "quasipoisson"))
    stop("`x` must be from poisson-family.", call. = F)

  # get predictions of outcome
  mu <- predict(x, type = "response")

  # get predicted zero-counts
  pred.zero <- round(sum(stats::dpois(x = 0, lambda = mu)))

  # get actual zero of response
  obs.zero <- sum(resp_val(x) == 0)

  # proportion
  structure(class = "sjstats_zcf", list(
    predicted.zeros = pred.zero,
    observed.zeros = obs.zero,
    ratio = pred.zero / obs.zero)
  )
}
