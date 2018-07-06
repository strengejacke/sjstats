#' @title Access information from model objects
#' @name pred_vars
#'
#' @description Several functions to retrieve information from model objects,
#'    like variable names, link-inverse function, model frame,
#'    model_family etc., in a tidy and consistent way.
#'
#' @param x A fitted model; for \code{var_names()}, \code{x} may also be a
#'    character vector.
#' @param fe.only Logical, if \code{TRUE} (default) and \code{x} is a mixed effects
#'    model, returns the model frame for fixed effects only.
#' @param multi.resp Logical, if \code{TRUE} and model is a multivariate response
#'    model from a \code{brmsfit} object, then a list of family-information
#'    for each regression is returned.
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
#'    removed. \code{model_family()} returns a list with information about the
#'    model family (see 'Details').
#'
#' @details \code{model_family()} returns a list with information about the
#'    model family for many different model objects. Following information
#'    is returned, where all values starting with \code{is_} are logicals.
#'    \itemize{
#'      \item \code{is_bin}: family is binomial (but not negative binomial)
#'      \item \code{is_pois}: family is either poisson or negative binomial
#'      \item \code{is_negbin}: family is negative binomial
#'      \item \code{is_logit}: model has logit link
#'      \item \code{is_linear}: family is gaussian
#'      \item \code{is_ordinal}: family is ordinal or cumulative link
#'      \item \code{is_zeroinf}: model has zero-inflation component
#'      \item \code{is_multivariate}: model is a multivariate response model (currently only works for \emph{brmsfit} objects)
#'      \item \code{link.fun}: the link-function
#'      \item \code{family}: the family-object
#'    }
#'    \code{model_frame()} slighty differs from \code{model.frame()}, especially
#'    for spline terms. Where \code{model.frame()} returns a matrix for splines,
#'    \code{model_frame()} returns the data of the original variable and uses
#'    the same column name as in the \code{data}-argument from the model-function.
#'    This makes it easier, for instance, to get data that should be used as new
#'    data in \code{predict()}. See 'Examples'.
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
#' outcome <- as.numeric(outcome)
#' m <- glm(counts ~ log(outcome) + as.factor(treatment), family = poisson())
#' var_names(m)
#'
#' # model.frame and model_frame behave slightly different
#' library(splines)
#' m <- lm(neg_c_7 ~ e42dep + ns(c160age), data = efc)
#' head(model.frame(m))
#' head(model_frame(m))
#'
#' @importFrom purrr flatten_chr map
#' @importFrom stats formula terms
#' @export
pred_vars <- function(x) {

  fm <- stats::formula(x)

  if (inherits(x, "brmsfit")) {
    if (!is.null(fm$responses)) {
      av <- fm$forms %>%
        purrr::map(~ all.vars(stats::formula(.x)[[3L]])) %>%
        purrr::flatten_chr() %>%
        unique()
    } else
      av <- all.vars(fm$formula[[3L]])
  } else
    av <- all.vars(fm[[3L]])

  if (length(av) == 1 && av == ".")
    av <- all.vars(stats::terms(x))

  av
}

#' @importFrom stats formula
#' @rdname pred_vars
#' @export
resp_var <- function(x) {
  if (inherits(x, "brmsfit")) {
    if (is.null(stats::formula(x)$responses)) {
      deparse(stats::formula(x)$formula[[2L]])
    } else {
      purrr::map_chr(
        x$formula$forms,
        ~ stats::formula(.x)[[2L]] %>% all.vars()
      )
      # stats::formula(x)$responses
    }
  } else
    deparse(stats::formula(x)[[2L]])
}


#' @rdname pred_vars
#' @importFrom nlme getResponse
#' @export
resp_val <- function(x) {
  if (inherits(x, c("lme", "gls")))
    as.vector(nlme::getResponse(x))
  else if (inherits(x, "brmsfit") && !is.null(stats::formula(x)$responses))
    as.vector(model_frame(x)[, var_names(resp_var(x))])
  else
    as.vector(model_frame(x)[[var_names(resp_var(x))]])
}


#' @rdname pred_vars
#' @importFrom stats family binomial gaussian make.link
#' @export
link_inverse <- function(x) {

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

  if (inherits(x, c("truncreg", "coxph"))) {
    il <- NULL
  } else if (inherits(x, c("hurdle", "zeroinfl"))) {
    il <- x$linkinv
  } else if (inherits(x, c("lme", "plm", "gls", "lm", "lmRob")) && !inherits(x, "glm")) {
    il <- stats::gaussian(link = "identity")$linkinv
  } else if (inherits(x, "betareg")) {
    il <- x$link$mean$linkinv
  } else if (inherits(x, c("vgam", "vglm"))) {
    il <- x@family@linkinv
  } else if (inherits(x, "stanmvreg")) {
    ## TODO save different family types for stan multivariate reponse models
    fam <- stats::family(x)
    fam <- fam[[1]]
  } else if (inherits(x, "brmsfit")) {
    fam <- stats::family(x)

    ## TODO save different family types for brms multivariate reponse models

    # in case of multivariate response models for brms, we just take the
    # information from the first model
    if (!is.null(stats::formula(x)$response))
      fam <- fam[[1]]

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
  } else if (inherits(x, "polr")) {
    link <- x$method
    if (link == "logistic") link <- "logit"
    il <- stats::make.link(link)$linkinv
  } else if (inherits(x, c("clm", "clmm"))) {
    il <- stats::make.link(x$link)$linkinv
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


#' @rdname pred_vars
#' @importFrom stats model.frame formula getCall na.omit
#' @importFrom prediction find_data
#' @importFrom purrr map_lgl map
#' @importFrom dplyr select bind_cols
#' @importFrom tibble as_tibble
#' @export
model_frame <- function(x, fe.only = TRUE) {
  # we may store model weights here later
  mw <- NULL

  if (inherits(x, c("merMod", "lmerMod", "glmerMod", "nlmerMod", "merModLmerTest")))
    fitfram <- stats::model.frame(x, fixed.only = fe.only)
  else if (inherits(x, "lme"))
    fitfram <- x$data
  else if (inherits(x, c("vgam", "gee", "gls")))
    fitfram <- prediction::find_data(x)
  else if (inherits(x, "Zelig-relogit"))
    fitfram <- get_zelig_relogit_frame(x)
  else if (inherits(x, "vglm")) {
    if (!length(x@model)) {
      env <- environment(x@terms$terms)
      if (is.null(env)) env <- parent.frame()
      fcall <- x@call
      fcall$method <- "model.frame"
      fcall$smart <- FALSE
      fitfram <- eval(fcall, env, parent.frame())
    } else {
      fitfram <- x@model
    }
  } else
    fitfram <- stats::model.frame(x)


  # clean 1-dimensional matrices

  fitfram <- purrr::modify_if(fitfram, is.matrix, function(x) {
    if (dim(x)[2] == 1)
      as.vector(x)
    else
      x
  })


  # check if we have any matrix columns, e.g. from splines

  mc <- purrr::map_lgl(fitfram, is.matrix)


  # don't change response value, if it's a matrix
  # bound with cbind()

  if (mc[1] && resp_var(x) == colnames(fitfram)[1]) mc[1] <- FALSE


  # if we have any matrix columns, we remove them from original
  # model frame and convert them to regular data frames, give
  # proper column names and bind them back to the original model frame

  if (any(mc)) {
    # try to get model data from environment
    md <- eval(stats::getCall(x)$data, environment(stats::formula(x)))

    # if data not found in environment, reduce matrix variables into regular vectors
    if (is.null(md)) {
      fitfram <- dplyr::bind_cols(purrr::map(fitfram, ~ tibble::as_tibble(.x)))
    } else {

      # get "matrix" terms and "normal" predictors, but exclude
      # response variable(s)

      fitfram_matrix <- dplyr::select(fitfram, -which(mc))
      spline.term <- get_vn_helper(names(which(mc)))
      other.terms <- get_vn_helper(colnames(fitfram_matrix))[-1]

      # now we have all variable names that we need from the original
      # data set

      needed.vars <- c(other.terms, spline.term)

      # if response is a matrix vector (e.g. multivariate response),
      # we need to include all response names as well, because else
      # rows may not match due to additional missings in the response variables

      if (is.matrix(fitfram[[1]])) {
        needed.vars <- c(dimnames(fitfram[[1]])[[2]], needed.vars)
      } else {
        needed.vars <- c(colnames(fitfram)[1], needed.vars)
      }

      # check model weights

      if ("(weights)" %in% needed.vars && !tibble::has_name(md, "(weights)")) {
        needed.vars <- needed.vars[-which(needed.vars == "(weights)")]
        mw <- fitfram[["(weights)"]]
      }


      fitfram <- stats::na.omit(dplyr::select(md, !! needed.vars))

      # add back model weights, if any
      if (!is.null(mw)) fitfram$`(weights)` <- mw
    }

  }

  # clean variable names
  cvn <- get_vn_helper(colnames(fitfram))

  # do we have duplicated names?
  dupes <- which(duplicated(cvn))
  if (!sjmisc::is_empty(dupes)) cvn[dupes] <- sprintf("%s.%s", cvn[dupes], 1:length(dupes))

  colnames(fitfram) <- cvn
  fitfram
}


#' @rdname pred_vars
#' @importFrom sjmisc str_contains is_empty
#' @importFrom stats family formula
#' @importFrom tidyselect starts_with
#' @export
model_family <- function(x, multi.resp = FALSE) {
  zero.inf <- FALSE
  multi.var <- FALSE

  # for gam-components from gamm4, add class attributes, so family
  # function works correctly
  if (inherits(x, "gam") && !inherits(x, c("glm", "lm")))
    class(x) <- c(class(x), "glm", "lm")

  # do we have glm? if so, get link family. make exceptions
  # for specific models that don't have family function
  if (inherits(x, c("lme", "plm", "gls", "truncreg", "lmRob"))) {
    fitfam <- "gaussian"
    logit.link <- FALSE
    link.fun <- "identity"
  } else if (inherits(x, c("vgam", "vglm"))) {
    faminfo <- x@family
    fitfam <- faminfo@vfamily[1]
    logit.link <- sjmisc::str_contains(faminfo@blurb, "logit")
    link.fun <- faminfo@blurb[3]
    if (!sjmisc::is_empty(tidyselect::starts_with("logit(", vars = link.fun)))
      link.fun <- "logit"
  } else if (inherits(x, c("zeroinfl", "hurdle"))) {
    fitfam <- "negative binomial"
    logit.link <- FALSE
    link.fun <- NULL
    zero.inf <- TRUE
  } else if (inherits(x, "betareg")) {
    fitfam <- "beta"
    logit.link <- x$link$mean$name == "logit"
    link.fun <- x$link$mean$linkfun
  } else if (inherits(x, "coxph")) {
    fitfam <- "survival"
    logit.link <- TRUE
    link.fun <- NULL
  } else {
    # "lrm"-object from pkg "rms" have no family method
    # so we construct a logistic-regression-family-object
    if (inherits(x, c("lrm", "polr", "logistf", "clmm", "clm", "multinom", "Zelig-relogit")))
      faminfo <- stats::binomial(link = "logit")
    else
      # get family info
      faminfo <- stats::family(x)

    # in case of multivariate response models for brms or rstanarm,
    # we just take the information from the first model
    if (inherits(x, "brmsfit") && !is.null(stats::formula(x)$response)) {
      multi.var <- TRUE
      if (!multi.resp) faminfo <- faminfo[[1]]
    }

    if (inherits(x, "stanmvreg")) {
      multi.var <- TRUE
      if (!multi.resp) faminfo <- faminfo[[1]]
    }


    if (multi.resp && multi.var) {
      return(purrr::map(faminfo, ~ make_family(
        x,
        .x$family,
        zero.inf,
        .x$link == "logit",
        TRUE,
        .x$link
      )))
    }


    fitfam <- faminfo$family
    logit.link <- faminfo$link == "logit"
    link.fun <- faminfo$link
  }

  make_family(x, fitfam, zero.inf, logit.link, multi.var, link.fun)
}


make_family <- function(x, fitfam, zero.inf, logit.link, multi.var, link.fun) {
  # create logical for family
  binom_fam <-
    fitfam %in% c("bernoulli", "binomial", "quasibinomial", "binomialff") |
    sjmisc::str_contains(fitfam, "binomial", ignore.case = TRUE)

  poisson_fam <-
    fitfam %in% c("poisson", "quasipoisson") |
    sjmisc::str_contains(fitfam, "poisson", ignore.case = TRUE)

  neg_bin_fam <-
    sjmisc::str_contains(fitfam, "negative binomial", ignore.case = T) |
    sjmisc::str_contains(fitfam, "nbinom", ignore.case = TRUE) |
    sjmisc::str_contains(fitfam, "negbinomial", ignore.case = TRUE) |
    sjmisc::str_contains(fitfam, "neg_binomial", ignore.case = TRUE)

  linear_model <- !binom_fam & !poisson_fam & !neg_bin_fam & !logit.link

  zero.inf <- zero.inf | sjmisc::str_contains(fitfam, "zero_inflated", ignore.case = T)

  is.ordinal <-
    inherits(x, c("polr", "clm", "clmm", "multinom")) |
    fitfam %in% c("cumulative", "cratio", "sratio", "acat")

  list(
    is_bin = binom_fam & !neg_bin_fam,
    is_pois = poisson_fam | neg_bin_fam,
    is_negbin = neg_bin_fam,
    is_logit = logit.link,
    is_linear = linear_model,
    is_zeroinf = zero.inf,
    is_ordinal = is.ordinal,
    is_multivariate = multi.var,
    link.fun = link.fun,
    family = fitfam
  )
}


#' @importFrom dplyr select
get_zelig_relogit_frame <- function(x) {
  vars <- c(resp_var(x), pred_vars(x))
  dplyr::select(x$data, !! vars)
}

#' @rdname pred_vars
#' @importFrom purrr map_chr
#' @export
var_names <- function(x) {
  if (is.character(x))
    get_vn_helper(x)
  else
    colnames(model_frame(x))
}


#' @importFrom sjmisc is_empty trim
#' @importFrom purrr map_chr
get_vn_helper <- function(x) {

  # return if x is empty
  if (sjmisc::is_empty(x)) return("")

  # for gam-smoothers/loess, remove s()- and lo()-function in column name
  # for survival, remove strata(), and so on...
  pattern <- c(
    "as.factor", "offset", "log", "lag", "diff", "lo", "bs", "ns", "t2", "te",
    "ti", "mi", "pspline", "poly", "strata", "scale", "s"
  )

  # do we have a "log()" pattern here? if yes, get capture region
  # which matches the "cleaned" variable name
  purrr::map_chr(1:length(x), function(i) {
    for (j in 1:length(pattern)) {
      if (pattern[j] == "offset") {
        x[i] <- sjmisc::trim(unique(sub("^offset\\(([^-+ )]*).*", "\\1", x[i])))
      } else {
        p <- paste0("^", pattern[j], "\\(([^,)]*).*")
        x[i] <- unique(sub(p, "\\1", x[i]))
      }
    }
    x[i]
  })
}
