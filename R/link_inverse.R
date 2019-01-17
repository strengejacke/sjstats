#' @rdname pred_vars
#' @importFrom stats family binomial gaussian make.link
#' @export
link_inverse <- function(x, multi.resp = FALSE, mv = FALSE) {

  if (!missing(multi.resp)) mv <- multi.resp

  # handle glmmTMB models
  if (inherits(x, "glmmTMB")) {
    ff <- stats::family(x)

    if ("linkinv" %in% names(ff))
      return(ff$linkinv)
    else if ("link" %in% names(ff) && is.character(ff$link))
      return(stats::make.link(ff$link)$linkinv)
    else
      return(match.fun("exp"))
  }


  # for gam-components from gamm4, add class attributes, so family
  # function works correctly

  if (inherits(x, "gam") && !inherits(x, c("glm", "lm")))
    class(x) <- c(class(x), "glm", "lm")


  # do we have glm? if so, get link family. make exceptions
  # for specific models that don't have family function

  if (inherits(x, c("truncreg", "coxph", "coxme"))) {
    il <- NULL
  } else if (inherits(x, c("zeroinfl", "hurdle", "zerotrunc"))) {
    il <- stats::make.link("log")$linkinv
  } else if (inherits(x, c("glmmPQL", "MixMod"))) {
    il <- x$family$linkinv
  } else if (inherits(x, c("lme", "plm", "lm_robust", "felm", "gls", "lm", "lmRob")) && !inherits(x, "glm")) {
    il <- stats::gaussian(link = "identity")$linkinv
  } else if (inherits(x, "betareg")) {
    il <- x$link$mean$linkinv
  } else if (inherits(x, c("vgam", "vglm"))) {
    il <- x@family@linkinv
  } else if (inherits(x, "stanmvreg")) {
    fam <- stats::family(x)
    if (mv) {
      il <- purrr::map(fam, ~ .x$linkinv)
    } else {
      fam <- fam[[1]]
      il <- fam$linkinv
    }
  } else if (inherits(x, "brmsfit")) {
    fam <- stats::family(x)
    if (!is.null(stats::formula(x)$response)) {
      if (mv) {
        il <- purrr::map(fam, ~ brms_link_inverse(.x))
      } else {
        fam <- fam[[1]]
        il <- brms_link_inverse(fam)
      }
    } else {
      il <- brms_link_inverse(fam)
    }
  } else if (inherits(x, "polr")) {
    link <- x$method
    if (link == "logistic") link <- "logit"
    il <- stats::make.link(link)$linkinv
  } else if (inherits(x, c("clm", "clmm"))) {
    il <- stats::make.link(x$link)$linkinv
  } else if (inherits(x, "clm2")) {
    il <- switch(
      x$link,
      logistic = ,
      probit = stats::make.link("logit")$linkinv,
      cloglog = ,
      loglog = stats::make.link("log")$linkinv,
      stats::make.link("logit")$linkinv
    )
  } else if (inherits(x, c("lrm", "logistf", "multinom", "Zelig-relogit"))) {
    # "lrm"-object from pkg "rms" have no family method
    # so we construct a logistic-regression-family-object
    il <- stats::make.link(link = "logit")$linkinv
  } else {
    # get family info
    il <- stats::family(x)$linkinv
  }

  il
}

brms_link_inverse <- function(fam) {
  # do we have custom families?
  if (!is.null(fam$family) && (is.character(fam$family) && fam$family == "custom")) {
    il <- stats::make.link(fam$link)$linkinv
  } else {
    if ("linkinv" %in% names(fam)) {
      il <- fam$linkinv
    } else if ("link" %in% names(fam) && is.character(fam$link)) {
      il <- stats::make.link(fam$link)$linkinv
    } else {
      ff <- get(fam$family, asNamespace("stats"))
      il <- ff(fam$link)$linkinv
    }
  }
  il
}

