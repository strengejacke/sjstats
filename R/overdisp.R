#' @title Check overdispersion of GL(M)M's
#' @name overdisp
#' @description \code{overdisp()} checks generalized linear (mixed) models for
#'              overdispersion, while \code{zero_count()} checks whether models
#'              from Poisson-families are over- or underfitting zero-counts in
#'              the outcome.
#'
#' @param x Fitted model of class \code{merMod}, \code{glmmTMB}, \code{glm},
#'    or \code{glm.nb} (package \pkg{MASS}).
#' @param tolerance The tolerance for the ratio of observed and predicted
#'          zeros to considered as over- or underfitting zero-counts. A ratio
#'          between 1 +/- \code{tolerance} is considered as OK, while a ratio
#'          beyond or below this treshold would indicate over- or underfitting.
#' @param ... Currently not used.
#'
#' @return For \code{overdisp()}, information on the overdispersion test; for
#'         \code{zero_count()}, the amount of predicted and observed zeros in
#'         the outcome, as well as the ratio between these two values.
#'
#' @note For overdispersoion test, a p-value < .05 indicates overdispersion.
#'       \cr \cr
#'       For \code{zero_count()}, a model that is underfitting zero-counts
#'       indicates a zero-inflation in the data, i.e. it is recommended to
#'       use negative binomial or zero-inflated models then.
#'
#' @details For \code{merMod}- and \code{glmmTMB}-objects, \code{overdisp()} is
#'          based on the code in the \href{http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html}{GLMM FAQ},
#'          section \emph{How can I deal with overdispersion in GLMMs?}.
#'          Note that this function only returns an \emph{approximate} estimate
#'          of an overdispersion parameter, and is probably inaccurate for
#'          zero-inflated mixed models (fitted with \code{glmmTMB}).
#'          \cr \cr
#'          The same code as above for mixed models is also used to check
#'          overdispersion for negative binomial models.
#'          \cr \cr
#'          For Poisson-models, the overdispersion test is based on the code
#'          from Gelman and Hill (2007), page 115.
#'
#' @references Bolker B et al. (2017): \href{http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html}{GLMM FAQ.}
#'  \cr \cr
#'  Gelman A, Hill J (2007) Data Analysis Using Regression and Multilevel/Hierarchical Models. Cambridge, New York: Cambridge University Press
#'
#' @examples
#' library(sjmisc)
#' data(efc)
#'
#' # response has many zero-counts, Poisson models
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
#' @export
overdisp <- function(x, ...) {
  UseMethod("overdisp")
}


#' @importFrom stats fitted nobs coef pchisq
#' @export
overdisp.glm <- function(x, ...) {
  # check if we have poisson
  if (!stats::family(x)$family %in% c("poisson", "quasipoisson"))
    stop("Model must be from Poisson-family.", call. = F)

  yhat <- stats::fitted(x)

  n <- stats::nobs(x)
  k <- length(stats::coef(x))

  zi <- (resp_val(x) - yhat) / sqrt(yhat)
  chisq <- sum(zi^2)
  ratio <-  chisq / (n - k)
  p.value <- stats::pchisq(chisq, df = n - k, lower.tail = FALSE)

  structure(
    class = "sj_overdisp",
    list(
      chisq = chisq,
      ratio = ratio,
      rdf = n - k,
      p = p.value
    )
  )
}


#' @export
overdisp.negbin <- function(x, ...) {
  overdisp.lme4(x)
}


#' @export
overdisp.merMod <- function(x, ...) {
  overdisp.lme4(x)
}


#' @export
overdisp.glmmTMB <- function(x, ...) {
  overdisp.lme4(x)
}


#' @importFrom stats df.residual residuals pchisq
overdisp.lme4 <- function(x) {
  rdf <- stats::df.residual(x)
  rp <- stats::residuals(x, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- stats::pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)

  structure(class = "sj_overdisp",
            list(
              chisq = Pearson.chisq,
              ratio = prat,
              rdf = rdf,
              p = pval
            ))
}


#' @rdname overdisp
#' @importFrom stats fitted dpois family
#' @export
zero_count <- function(x, tolerance = .05) {
  # check if we have poisson
  if (!stats::family(x)$family %in% c("poisson", "quasipoisson"))
    stop("Model must be from Poisson-family.", call. = F)

  # get predictions of outcome
  mu <- stats::fitted(x)

  # get predicted zero-counts
  pred.zero <- round(sum(stats::dpois(x = 0, lambda = mu)))

  # get actual zero of response
  obs.zero <- sum(resp_val(x) == 0)

  # proportion
  structure(class = "sj_zcf", list(
    predicted.zeros = pred.zero,
    observed.zeros = obs.zero,
    ratio = pred.zero / obs.zero,
    tolerance = tolerance)
  )
}
