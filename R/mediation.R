#' @rdname hdi
#' @export
mediation <- function(x, ...) {
  UseMethod("mediation")
}


#' @rdname hdi
#' @importFrom tibble has_name
#' @importFrom purrr map
#' @importFrom stats formula
#' @importFrom dplyr pull bind_cols
#' @export
mediation.brmsfit <- function(x, treatment, mediator, prob = .9, typical = "median", ...) {
  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  # only one HDI interval
  if (length(prob) > 1) prob <- prob[1]


  dv <- resp_var(x)

  if (missing(mediator)) {
    pv <- pred_vars(x)
    mediator <- pv[pv %in% dv]
  }

  if (missing(treatment)) {
    pvs <- purrr::map(
      x$formula$forms,
      ~ stats::formula(.x)[[3L]] %>% all.vars()
    )

    treatment <- pvs[[1]][pvs[[1]] %in% pvs[[2]]][1]

    # check for categorical. if user has not specified a treatment variable
    # and this variable is categorical, the posterior samples contain the
    # samples from each category of the treatment variable - so we need to
    # fix the variable name

    mf <- model_frame(x)
    if (tibble::has_name(mf, treatment)) {
      check_fac <- mf[[treatment]]
      if (is.factor(check_fac)) {
        treatment <- sprintf("%s%s", treatment, levels(check_fac)[nlevels(check_fac)])
      }
    }
  }


  mediator.model <- which(dv == mediator)
  treatment.model <- which(dv != mediator)

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
  prop.mediated <- typical_value(eff.indirect, fun = typical) / typical_value(eff.total, fun = typical)
  prop.se <- diff(hdi(eff.indirect / eff.total, prob = prob) / 2)
  prop.hdi <- prop.mediated + c(-1, 1) * prop.se

  res <- data.frame(
    effect = c("direct", "indirect", "mediator", "total", "proportion mediated"),
    value = c(
      typical_value(eff.direct, fun = typical),
      typical_value(eff.indirect, fun = typical),
      typical_value(eff.mediator, fun = typical),
      typical_value(eff.total, fun = typical),
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
