#' @title Get inverse link function from model objects
#' @name link_inverse
#'
#' @description Get inverse link function from model objects.
#'
#' @param model A model object.
#'
#' @return If known, the inverse link function from \code{model}; else \code{NULL}
#'         for those models where the inverse link function can't be identified.
#'
#' @examples
#' data(efc)
#' m <- lm(neg_c_7 ~ e42dep + c172code, efc)
#' link_inverse(m)(2.3)
#'
#' # example from ?stats::glm
#' counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
#' outcome <- gl(3, 1, 9)
#' treatment <- gl(3, 3)
#' m <- glm(counts ~ outcome + treatment, family = poisson())
#'
#' link_inverse(m)(.3)
#' # same as
#' exp(.3)
#'
#'
#' @importFrom stats family binomial gaussian
#' @export
link_inverse <- function(model) {

  # handle glmmTMB models
  if (inherits(model, "glmmTMB")) {
    ff <- stats::family(model)

    if ("linkinv" %in% names(ff))
      return(ff$linkinv)
    else
      return(match.fun("exp"))
  }

  # do we have glm? if so, get link family. make exceptions
  # for specific models that don't have family function
  if (inherits(model, c("truncreg", "coxph"))) {
    il <- NULL
  } else if (inherits(model, c("hurdle", "zeroinfl"))) {
    il <- model$linkinv
  } else if (inherits(model, c("lme", "plm", "gls", "lm")) && !inherits(model, "glm")) {
    il <- stats::gaussian(link = "identity")$linkinv
  } else if (inherits(model, "betareg")) {
    il <- model$link$mean$linkinv
  } else if (inherits(model, c("lrm", "polr"))) {
    # "lrm"-object from pkg "rms" have no family method
    # so we construct a logistic-regression-family-object
    il <- stats::binomial(link = "logit")$linkinv
  } else {
    # get family info
    il <- stats::family(model)$linkinv
  }

  il
}
