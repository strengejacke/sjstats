#' @title Convergence test for mixed effects models
#' @name converge_ok
#'
#' @description This function provides an alternative convergence test for
#'                \code{\link[lme4]{merMod}}-objects.
#'
#' @param x A \code{merMod}-object.
#' @param tolerance Indicates up to which value the convergence result is
#'          accepted. The smaller \code{tolerance} is, the stricter the test
#'          will be.
#'
#' @return Logical vector, \code{TRUE} if convergence is fine, \code{FALSE}
#'           if convergence is suspicious. Additionally, the convergence
#'           value is returned as return value's name.
#'
#' @details This function provides an alternative convergence test for
#'                \code{\link[lme4]{merMod}}-objects, as discussed
#'                \href{https://github.com/lme4/lme4/issues/120}{here}
#'                and suggested by Ben Bolker in
#'                \href{https://github.com/lme4/lme4/issues/120#issuecomment-39920269}{this comment}.
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
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package `Matrix` needed for this function to work. Please install it.", call. = FALSE)
  }

  # is 'x' an lmer object?
  if (is_merMod(x)) {
    relgrad <- with(x@optinfo$derivs, Matrix::solve(Hessian, gradient))
    # copy logical value, TRUE if convergence is OK
    retval <- max(abs(relgrad)) < tolerance
    # copy convergence value
    names(retval) <- max(abs(relgrad))
    # return result
    return(retval)
  } else {
    warning("`x` must be a `merMod` object.", call. = F)
  }
}
