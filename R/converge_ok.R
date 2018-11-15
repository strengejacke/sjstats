#' @title Convergence test for mixed effects models
#' @name converge_ok
#'
#' @description \code{converge_ok()} provides an alternative convergence test for
#'                \code{\link[lme4]{merMod}}-objects; \code{is_singular()} checks
#'                post-fitting convergence warnings. If the model fit is singular,
#'                warning about negative eigenvalues of the Hessian can most likely
#'                be ignored.
#'
#' @param x A \code{merMod}-object. For \code{is_singluar()}, may also be a
#'          \code{glmmTMB}-object.
#' @param tolerance Indicates up to which value the convergence result is
#'          accepted. The smaller \code{tolerance} is, the stricter the test
#'          will be.
#' @param ... Currently not used.
#'
#' @return For \code{converge_ok()}, a logical vector, which is \code{TRUE} if
#'           convergence is fine and \code{FALSE} if convergence is suspicious.
#'           Additionally, the convergence value is returned as return value's name.
#'           \code{is_singluar()} returns \code{TRUE} if the model fit is singular.
#'
#' @details \code{converge_ok()} provides an alternative convergence test for
#'   \code{\link[lme4]{merMod}}-objects, as discussed
#'   \href{https://github.com/lme4/lme4/issues/120}{here}
#'   and suggested by Ben Bolker in
#'   \href{https://github.com/lme4/lme4/issues/120#issuecomment-39920269}{this comment}.
#'   \cr \cr
#'   If a model is "singular", this means that some dimensions of the variance-covariance
#'   matrix have been estimated as exactly zero. \code{is_singular()} checks if
#'   a model fit is singular, and can be used in case of post-fitting convergence
#'   warnings, such as warnings about negative eigenvalues of the Hessian. If the fit
#'   is singular (i.e. \code{is_singular()} returns \code{TRUE}), these warnings
#'   can most likely be ignored.
#'   \cr \cr
#'   There is no gold-standard about how to deal with singularity and which
#'   random-effects specification to choose. Beside using fully Bayesian methods
#'   (with informative priors), proposals in a frequentist framework are:
#'   \itemize{
#'   \item avoid fitting overly complex models, such that the variance-covariance matrices can be estimated precisely enough (\cite{Matuschek et al. 2017})
#'   \item use some form of model selection to choose a model that balances predictive accuracy and overfitting/type I error (\cite{Bates et al. 2015}, \cite{Matuschek et al. 2017})
#'   \item \dQuote{keep it maximal}, i.e. fit the most complex model consistent with the experimental design, removing only terms required to allow a non-singular fit (\cite{Barr et al. 2013})
#'   }
#'
#' @references \itemize{
#'   \item Bates D, Kliegl R, Vasishth S, Baayen H. Parsimonious Mixed Models. arXiv:1506.04967, June 2015.
#'   \item Barr DJ, Levy R, Scheepers C, Tily HJ. Random effects structure for confirmatory hypothesis testing: Keep it maximal. Journal of Memory and Language, 68(3):255–278, April 2013.
#'   \item Matuschek H, Kliegl R, Vasishth S, Baayen H, Bates D. Balancing type I error and power in linear mixed models. Journal of Memory and Language, 94:305–315, 2017.
#'   }
#'
#' @examples
#' library(sjmisc)
#' library(lme4)
#' data(efc)
#' # create binary response
#' efc$hi_qol <- dicho(efc$quol_5)
#' # prepare group variable
#' efc$grp = as.factor(efc$e15relat)
#' # data frame for fitted model
#' mydf <- data.frame(hi_qol = as.factor(efc$hi_qol),
#'                    sex = as.factor(efc$c161sex),
#'                    c12hour = as.numeric(efc$c12hour),
#'                    neg_c_7 = as.numeric(efc$neg_c_7),
#'                    grp = efc$grp)
#' # fit glmer
#' fit <- glmer(hi_qol ~ sex + c12hour + neg_c_7 + (1|grp),
#'              data = mydf, family = binomial("logit"))
#'
#' converge_ok(fit)
#'
#' @importFrom Matrix solve
#' @export
converge_ok <- function(x, tolerance = 0.001) {
  # check for package availability
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package `Matrix` needed for this function to work. Please install it.", call. = FALSE)

  # is 'x' an lmer object?
  if (!is_merMod(x))
    warning("`x` must be a `merMod` object.", call. = F)


  relgrad <- with(x@optinfo$derivs, Matrix::solve(Hessian, gradient))

  # copy logical value, TRUE if convergence is OK
  retval <- max(abs(relgrad)) < tolerance
  # copy convergence value
  names(retval) <- max(abs(relgrad))

  # return result
  retval
}


#' @rdname converge_ok
#' @export
is_singular <- function(x, tolerance = 1e-5, ...) {
  UseMethod("is_singular")
}

#' @importFrom lme4 getME
#' @export
is_singular.merMod <- function(x, tolerance = 1e-5, ...) {
  theta <- lme4::getME(x, "theta")
  # diagonal elements are identifiable because they are fitted
  #  with a lower bound of zero ...
  diag.element <- lme4::getME(x, "lower") == 0
  any(abs(theta[diag.element]) < tolerance)
}

#' @importFrom lme4 VarCorr
#' @export
is_singular.glmmTMB <- function(x, tolerance = 1e-5, ...) {
  vc <- collapse_cond(lme4::VarCorr(x))
  any(sapply(vc, function(.x) any(abs(diag(.x)) < tolerance)))
}
