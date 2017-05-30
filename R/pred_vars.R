#' @title Get predictor and response variables from models
#' @name pred_vars
#'
#' @description \code{pred_vars()} and \code{resp_var()} return the names of a
#'                model's predictor or response variables as character vector.
#'                \cr \cr
#'                \code{resp_val()} returns the values of the model's response
#'                vector.
#'
#' @param x A fitted model.
#'
#' @return The name(s) of the response or predictor variables from \code{x}
#'           as character vector; or the values from \code{x}'s response vector.
#'
#'
#' @examples
#' data(efc)
#' fit <- lm(neg_c_7 ~ e42dep + c161sex, data = efc)
#'
#' pred_vars(fit)
#' resp_var(fit)
#'
#' resp_val(fit)
#'
#' @importFrom stats formula
#' @export
pred_vars <- function(x) {
  all.vars(stats::formula(x)[[3L]])
}

#' @rdname pred_vars
#' @export
resp_var <- function(x) {
  deparse(stats::formula(x)[[2L]])
}

#' @rdname pred_vars
#' @importFrom nlme getResponse
#' @importFrom stats model.frame
#' @export
resp_val <- function(x) {
  if (inherits(x, "lme"))
    as.vector(nlme::getResponse(x))
  else
    as.vector(stats::model.frame(x)[[resp_var(x)]])
}
