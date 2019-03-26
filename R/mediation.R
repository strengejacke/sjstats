#' @rdname hdi
#' @export
mediation <- function(x, ...) {
  UseMethod("mediation")
}


#' @rdname hdi
#' @importFrom purrr map
#' @importFrom stats formula
#' @importFrom dplyr pull bind_cols
#' @importFrom sjmisc typical_value
#' @export
mediation.brmsfit <- function(x, treatment, mediator, prob = .9, typical = "median", ...) {
  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  # only one HDI interval
  if (length(prob) > 1) prob <- prob[1]

  # check for binary response. In this case, user should rescale variables
  fitinfo <- model_family(x, mv = TRUE)
  if (any(purrr::map_lgl(fitinfo, ~ .x$is_bin))) {
    message("One of moderator or outcome is binary, so direct and indirect effects may be on different scales. Consider rescaling model predictors, e.g. with `sjmisc::std()`.")
  }


  dv <- insight::find_response(x, combine = TRUE)
  fixm <- FALSE

  if (missing(mediator)) {
    pv <- insight::find_predictors(x)
    mediator <- pv[pv %in% dv]
    fixm <- TRUE
  }

  if (missing(treatment)) {
    pvs <- purrr::map(
      x$formula$forms,
      ~ stats::formula(.x)[[3L]] %>% all.vars()
    )

    treatment <- pvs[[1]][pvs[[1]] %in% pvs[[2]]][1]
    treatment <- fix_factor_name(x, treatment)
  }


  mediator.model <- which(dv == mediator)
  treatment.model <- which(dv != mediator)

  if (fixm) mediator <- fix_factor_name(x, mediator)

  # brms removes underscores from variable names when naming estimates
  # so we need to fix variable names here

  dv <- names(dv)


  # Direct effect: coef(treatment) from model_y_treatment
  coef_treatment <- sprintf("b_%s_%s", dv[treatment.model], treatment)
  eff.direct <- x %>%
    brms::posterior_samples(pars = coef_treatment, exact_match = TRUE) %>%
    dplyr::pull(1)

  # Mediator effect: coef(mediator) from model_y_treatment
  coef_mediator <- sprintf("b_%s_%s", dv[treatment.model], mediator)
  eff.mediator <- x %>%
    brms::posterior_samples(pars = coef_mediator, exact_match = TRUE) %>%
    dplyr::pull(1)

  # Indirect effect: coef(treament) from model_m_mediator * coef(mediator) from model_y_treatment
  coef_indirect <- sprintf("b_%s_%s", dv[mediator.model], treatment)
  tmp.indirect <- brms::posterior_samples(x, pars = c(coef_indirect, coef_mediator), exact_match = TRUE)
  eff.indirect <- tmp.indirect[[coef_indirect]] * tmp.indirect[[coef_mediator]]

  # Total effect
  eff.total <- eff.indirect + eff.direct

  # proportion mediated: indirect effect / total effect
  prop.mediated <- sjmisc::typical_value(eff.indirect, fun = typical) / sjmisc::typical_value(eff.total, fun = typical)
  prop.se <- diff(hdi(eff.indirect / eff.total, prob = prob) / 2)
  prop.hdi <- prop.mediated + c(-1, 1) * prop.se

  res <- data_frame(
    effect = c("direct", "indirect", "mediator", "total", "proportion mediated"),
    value = c(
      sjmisc::typical_value(eff.direct, fun = typical),
      sjmisc::typical_value(eff.indirect, fun = typical),
      sjmisc::typical_value(eff.mediator, fun = typical),
      sjmisc::typical_value(eff.total, fun = typical),
      prop.mediated
    )
  ) %>% dplyr::bind_cols(
    as.data.frame(rbind(
      hdi(eff.direct, prob = prob),
      hdi(eff.indirect, prob = prob),
      hdi(eff.mediator, prob = prob),
      hdi(eff.total, prob = prob),
      prop.hdi
    ))
  )

  colnames(res) <- c("effect", "value", "hdi.low", "hdi.high")

  attr(res, "prob") <- prob
  attr(res, "treatment") <- treatment
  attr(res, "mediator") <- mediator
  attr(res, "response") <- dv[treatment.model]
  attr(res, "formulas") <- lapply(x$formula$forms, function(x) as.character(x[1]))

  class(res) <- c("sj_mediation", class(res))
  res
}


fix_factor_name <- function(model, variable) {
  # check for categorical. if user has not specified a treatment variable
  # and this variable is categorical, the posterior samples contain the
  # samples from each category of the treatment variable - so we need to
  # fix the variable name

  mf <- model_frame(model)
  if (obj_has_name(mf, variable)) {
    check_fac <- mf[[variable]]
    if (is.factor(check_fac)) {
      variable <- sprintf("%s%s", variable, levels(check_fac)[nlevels(check_fac)])
    } else if (is.logical(check_fac)) {
      variable <- sprintf("%sTRUE", variable)
    }
  }

  variable
}
