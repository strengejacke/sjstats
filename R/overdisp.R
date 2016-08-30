#' @title Check overdispersion of GL(M)M's
#' @name overdisp
#' @description This function checks generalized linear (mixed) models for
#'                overdispersion.
#'
#' @param x Fitted GLMM (\code{\link[lme4]{merMod}}-class) or \code{glm} model.
#' @param trafo A specification of the alternative, can be numeric or a
#'          (positive) function or \code{NULL} (the default). See 'Details'
#'          in \code{\link[AER]{dispersiontest}} in package \CRANpkg{AER}. Does not
#'          apply to \code{merMod} objects.
#'
#' @return Information on the overdispersion test.
#'
#' @note The interpretation of the returned p-value differs between GLM and
#'         GLMM. For GLMs, a p-value < .05 indicates overdispersion, while
#'         for GLMMs, a p-value > .05 indicates overdispersion.
#'
#' @details For \code{merMod}-objects, this function is based on the code in the
#'            \href{http://glmm.wikidot.com/faq}{DRAFT r-sig-mixed-models FAQ},
#'            section \emph{How can I deal with overdispersion in GLMMs?}.
#'            Note that this function only returns an \emph{approximate} estimate
#'            of an overdispersion parameter.
#'            \cr \cr
#'            For \code{glm}'s, this function simply wraps the \code{dispersiontest}
#'            from the \pkg{AER}-package.
#'
#' @references \href{http://glmm.wikidot.com/faq}{DRAFT r-sig-mixed-models FAQ}
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
#'
#' library(lme4)
#' efc$e15relat <- to_factor(efc$e15relat)
#' fit <- glmer(tot_sc_e ~ neg_c_7 + e42dep + c160age + (1 | e15relat),
#'              data = efc, family = poisson)
#' overdisp(fit)
#'
#'
#' @importFrom stats df.residual residuals pchisq
#' @export
overdisp <- function(x, trafo = NULL) {
  if (sjmisc::str_contains(class(x), "merMod", ignore.case = TRUE))
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
  return(invisible(result))
}


overdisp.lme4 <- function(x) {
  # check object class
  if (any(class(x) == "glmerMod")) {
    rdf <- stats::df.residual(x)
    rp <- stats::residuals(x, type = "pearson")
    Pearson.chisq <- sum(rp ^ 2)
    prat <- Pearson.chisq / rdf
    pval <- stats::pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
    cat(sprintf("\n        Overdispersion test\n\ndispersion ratio = %.4f, p-value = %.4f\n\n",
                prat, pval))
    if (pval > 0.05)
      message("No overdispersion detected.")
    else
      message("Overdispersion detected.")
    return(invisible(list(
      chisq = Pearson.chisq,
      ratio = prat,
      rdf = rdf,
      p = pval
    )))
  } else {
    warning("This method currently only supports `glmer` fitted models.", call. = F)
  }
}
