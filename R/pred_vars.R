#' @title Access information from model objects
#' @name pred_vars
#'
#' @description Several functions to retrieve information from model objects,
#'    like variable names, link-inverse function, model frame etc.
#'
#' @param x A fitted model.
#' @param fe.only Logical, if \code{TRUE} (default) and \code{x} is a mixed effects
#'    model, returns the model frame for fixed effects only.
#'
#' @return For \code{pred_vars()} and \code{resp_var()}, the name(s) of the
#'    response or predictor variables from \code{x} as character vector.
#'    \code{resp_val()} returns the values from \code{x}'s response vector.
#'    \code{link_inverse()} returns, if known, the inverse link function from
#'    \code{x}; else \code{NULL} for those models where the inverse link function
#'    can't be identified. \code{model_frame()} is similar to \code{model.frame()},
#'    but should also work for model objects that don't have a S3-generic for
#'    \code{model.frame()}. \code{var_names()} returns the "cleaned" variable
#'    names, i.e. things like \code{s()} for splines or \code{log()} are
#'    removed.
#'
#'
#' @examples
#' data(efc)
#' fit <- lm(neg_c_7 ~ e42dep + c161sex, data = efc)
#'
#' pred_vars(fit)
#' resp_var(fit)
#' resp_val(fit)
#'
#' link_inverse(fit)(2.3)
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
#' @importFrom stats formula terms
#' @export
pred_vars <- function(x) {
  if (inherits(x, "brmsfit"))
    av <- all.vars(stats::formula(x)$formula[[3L]])
  else
    av <- all.vars(stats::formula(x)[[3L]])

  if (length(av) == 1 && av == ".")
    av <- all.vars(stats::terms(x))

  av
}

#' @rdname pred_vars
#' @export
resp_var <- function(x) {
  if (inherits(x, "brmsfit"))
    deparse(stats::formula(x)$formula[[2L]])
  else
    deparse(stats::formula(x)[[2L]])
}


#' @rdname pred_vars
#' @importFrom nlme getResponse
#' @export
resp_val <- function(x) {
  if (inherits(x, c("lme", "gls")))
    as.vector(nlme::getResponse(x))
  else
    as.vector(model_frame(x)[[resp_var(x)]])
}


#' @rdname pred_vars
#' @importFrom stats family binomial gaussian
#' @export
link_inverse <- function(x) {

  # handle glmmTMB models
  if (inherits(x, "glmmTMB")) {
    ff <- stats::family(x)

    if ("linkinv" %in% names(ff))
      return(ff$linkinv)
    else
      return(match.fun("exp"))
  }

  # do we have glm? if so, get link family. make exceptions
  # for specific models that don't have family function
  if (inherits(x, c("truncreg", "coxph"))) {
    il <- NULL
  } else if (inherits(x, c("hurdle", "zeroinfl"))) {
    il <- x$linkinv
  } else if (inherits(x, c("lme", "plm", "gls", "lm")) && !inherits(x, "glm")) {
    il <- stats::gaussian(link = "identity")$linkinv
  } else if (inherits(x, "betareg")) {
    il <- x$link$mean$linkinv
  } else if (inherits(x, "vgam")) {
    il <- x@family@linkinv
  } else if (inherits(x, c("lrm", "polr"))) {
    # "lrm"-object from pkg "rms" have no family method
    # so we construct a logistic-regression-family-object
    il <- stats::binomial(link = "logit")$linkinv
  } else {
    # get family info
    il <- stats::family(x)$linkinv
  }

  il
}


#' @rdname pred_vars
#' @importFrom stats model.frame formula getCall
#' @importFrom prediction find_data
#' @importFrom purrr map_lgl
#' @importFrom dplyr select bind_cols one_of
#' @export
model_frame <- function(x, fe.only = TRUE) {
  if (inherits(x, c("merMod", "lmerMod", "glmerMod", "nlmerMod", "merModLmerTest")))
    fitfram <- stats::model.frame(x, fixed.only = fe.only)
  else if (inherits(x, "lme"))
    fitfram <- x$data
  else if (inherits(x, c("vgam", "gee", "gls")))
    fitfram <- prediction::find_data(x)
  else
    fitfram <- stats::model.frame(x)

  # check if we have any matrix columns, e.g. from splines
  mc <- purrr::map_lgl(fitfram, is.matrix)

  # if we have any matrix columns, we remove them from original
  # model frame and convert them to regular data frames, give
  # proper column names and bind them back to the original model frame
  if (any(mc)) {
    fitfram <- dplyr::select(fitfram, -which(mc))
    spline.term <- var_names(names(which(mc)))
    # try to get model data from environment
    md <- eval(stats::getCall(x)$data, environment(stats::formula(x)))
    # bind spline terms to model frame
    fitfram <- dplyr::bind_cols(fitfram, dplyr::select(md, dplyr::one_of(spline.term)))
  }

  # clean variable names
  colnames(fitfram) <- var_names(colnames(fitfram))

  fitfram
}


#' @rdname pred_vars
#' @importFrom purrr map_chr
#' @export
var_names <- function(x) {
  if (is.character(x))
    get_vn_helper(x)
  else
    get_vn_helper(colnames(model_frame(x)))
}


get_vn_helper <- function(x) {
  # for gam-smoothers/loess, remove s()- and lo()-function in column name
  # for survival, remove strata()
  pattern <- c("log", "s", "lo", "bs", "poly", "strata")

  # do we have a "log()" pattern here? if yes, get capture region
  # which matches the "cleaned" variable name
  purrr::map_chr(1:length(x), function(i) {
    for (j in 1:length(pattern)) {
      p <- paste0("^", pattern[j], "\\(([^,)]*).*")
      x[i] <- unique(sub(p, "\\1", x[i]))
    }
    x[i]
  })
}
