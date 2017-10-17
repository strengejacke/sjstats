#' @title Standard Error for variables or coefficients
#' @name se
#' @description Compute standard error for a variable, for all variables
#'                of a data frame, for joint random and fixed effects
#'                coefficients of (non-/linear) mixed models, the adjusted
#'                standard errors for generalized linear (mixed) models, or
#'                for intraclass correlation coefficients (ICC).
#'
#' @param x (Numeric) vector, a data frame, an \code{lm} or \code{glm}-object,
#'          a \code{merMod}-object as returned by the functions from the
#'          \pkg{lme4}-package, an ICC object (as obtained by the
#'          \code{\link{icc}}-function) or a list with estimate and p-value.
#'          For the latter case, the list must contain elements named
#'          \code{estimate} and \code{p.value} (see 'Examples' and 'Details').
#' @param nsim Numeric, the number of simulations for calculating the
#'          standard error for intraclass correlation coefficients, as
#'          obtained by the \code{\link{icc}}-function.
#' @param type Type of standard errors for generalized linear mixed models.
#'          \code{type = "fe"} returns the standard errors for fixed effects,
#'          based on the delta-method-approximation. \code{type = "re"} returns
#'          the standard errors for joint random and fixed effects, which are
#'          on the scale of the link function. See 'Details'.
#'
#' @return The standard error of \code{x}.
#'
#' @note Computation of standard errors for coefficients of mixed models
#'         is based \href{http://stackoverflow.com/questions/26198958/extracting-coefficients-and-their-standard-error-from-lme}{on this code}.
#'         Standard errors for generalized linear (mixed) models, if
#'         \code{type = "re"}, are approximations based on the delta
#'         method (Oehlert 1992).
#'         \cr \cr
#'         A remark on standard errors:
#'         \dQuote{Standard error represents variation in the point estimate, but
#'         confidence interval has usual Bayesian interpretation only with flat prior.}
#'         (Gelman 2017)
#'
#' @details For linear mixed models, and generalized linear mixed models \strong{with
#'            \code{type = "re"}}, this function computes the standard errors
#'            for joint (sums of) random and fixed effects coefficients (unlike
#'            \code{\link[arm]{se.coef}}, which returns the standard error
#'            for fixed and random effects separately). Hence, \code{se()}
#'            returns the appropriate standard errors for \code{\link[lme4]{coef.merMod}}.
#'            \cr \cr
#'            For generalized linear models or generalized linear mixed models,
#'            approximated standard errors, using the delta method for transformed
#'            regression parameters are returned (Oehlert 1992). For generalized
#'            linear mixed models, by default, the standard errors refer to the
#'            fixed effects only. Use \code{type = "re"} to compute standard errors
#'            for joint random and fixed effects coefficients. However,
#'            computation for the latter \emph{is not} based on the delta method,
#'            so standard errors from \code{type = "re"} are on the scale of the
#'            link-function (and not back transformed).
#'            \cr \cr
#'            The standard error for the \code{\link{icc}} is based on bootstrapping,
#'            thus, the \code{nsim}-argument is required. See 'Examples'.
#'            \cr \cr
#'            \code{se()} also returns the standard error of an estimate (regression
#'            coefficient) and p-value, assuming a normal distribution to compute
#'            the z-score from the p-value (formula in short: \code{b / qnorm(p / 2)}).
#'            See 'Examples'.
#'
#' @references Oehlert GW. 1992. A note on the delta method. American Statistician 46(1).
#'             \cr \cr
#'             Gelman A 2017. How to interpret confidence intervals? \url{http://andrewgelman.com/2017/03/04/interpret-confidence-intervals/}
#'
#' @examples
#' # compute standard error for vector
#' se(rnorm(n = 100, mean = 3))
#'
#' # compute standard error for each variable in a data frame
#' data(efc)
#' se(efc[, 1:3])
#'
#' # compute standard error for merMod-coefficients
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' se(fit)
#'
#' # compute odds-ratio adjusted standard errors, based on delta method
#' # with first-order Taylor approximation.
#' data(efc)
#' efc$services <- sjmisc::dicho(efc$tot_sc_e, dich.by = 0)
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep,
#'            data = efc, family = binomial(link = "logit"))
#' se(fit)
#'
#' # compute odds-ratio adjusted standard errors for generalized
#' # linear mixed model, also based on delta method
#' library(lme4)
#' library(sjmisc)
#' # create binary response
#' sleepstudy$Reaction.dicho <- dicho(sleepstudy$Reaction, dich.by = "median")
#' fit <- glmer(Reaction.dicho ~ Days + (Days | Subject),
#'              data = sleepstudy, family = binomial("logit"))
#' se(fit)
#'
#' # compute standard error from regression coefficient and p-value
#' se(list(estimate = .3, p.value = .002))
#'
#' \dontrun{
#' # compute standard error of ICC for the linear mixed model
#' icc(fit)
#' se(icc(fit))
#'
#' # the standard error for the ICC can be computed manually in this way,
#' # taking the fitted model example from above
#' library(tidyverse)
#' dummy <- sleepstudy %>%
#'   # generate 100 bootstrap replicates of dataset
#'   bootstrap(100) %>%
#'   # run mixed effects regression on each bootstrap replicate
#'   # and compute ICC for each "bootstrapped" regression
#'   mutate(
#'     models = map(strap, ~lmer(Reaction ~ Days + (Days | Subject), data = .x)),
#'     icc = map_dbl(models, ~icc(.x))
#'   )
#'
#' # now compute SE and p-values for the bootstrapped ICC, values
#' # may differ from above example due to random seed
#' boot_se(dummy, icc)
#' boot_p(dummy, icc)}
#'
#'
#' @importFrom stats qnorm vcov
#' @importFrom broom tidy
#' @importFrom dplyr mutate select
#' @importFrom rlang .data
#' @importFrom purrr map_dbl
#' @export
se <- function(x, nsim = 100, type = c("fe", "re")) {
  # match arguments
  type <- match.arg(type)

  if (inherits(x, c("stanreg", "stanfit"))) {
    se_result <- x %>%
      broom::tidy() %>%
      dplyr::select(.data$term, .data$estimate, .data$std.error)
  } else if (inherits(x, c("lmerMod", "nlmerMod", "merModLmerTest"))) {
    # return standard error for (linear) mixed models
    se_result <- std_merMod(x)
  } else if (inherits(x, "icc.lme4")) {
    # we have a ICC object, so do bootstrapping and compute SE for ICC
    se_result <- std_e_icc(x, nsim)
  } else if (inherits(x, c("svyglm.nb", "svymle"))) {
    se_result <- x %>%
      tidy_svyglm.nb() %>%
      dplyr::select(.data$term, .data$estimate, .data$std.error)
  } else if (inherits(x, c("glm", "glmerMod"))) {
    # check type of se
    if (type == "fe") {
      # for glm, we want to exponentiate coefficients to get odds ratios, however
      # 'exponentiate'-argument currently not works for lme4-tidiers
      # so we need to do this manually for glmer's
      tm <- broom::tidy(x, effects = "fixed")
      tm$estimate <- exp(tm$estimate)

      # # for poisson family, we need a different delta method approach
      # if (get_glm_family(x)$is_pois) {
      #   # standard errors scaled using square root of Pearson
      #   # chi-squared dispersion
      #   pr <- sum(stats::residuals(x, type = "pearson") ^ 2)
      #   dispersion <- pr / x$df.residual
      #   sse <- sqrt(diag(as.matrix(stats::vcov(x)))) * sqrt(dispersion)
      #
      #   return(
      #     tm %>%
      #       # vcov for merMod returns a dpoMatrix-object, so we need
      #       # to coerce to regular matrix here.
      #       dplyr::mutate(std.error = sse) %>%
      #       dplyr::select_("term", "estimate", "std.error")
      #   )
      # }

      se_result <-
        tm %>%
          # vcov for merMod returns a dpoMatrix-object, so we need
          # to coerce to regular matrix here.
          dplyr::mutate(std.error = sqrt(.data$estimate ^ 2 * diag(as.matrix(stats::vcov(x))))) %>%
          dplyr::select(.data$term, .data$estimate, .data$std.error)
    } else {
      # return standard error for mixed models,
      # joint random and fixed effects
      se_result <- std_merMod(x)
    }
  } else if (inherits(x, "lm")) {
    # for convenience reasons, also return se for simple linear models
    se_result <- x %>%
      broom::tidy(effects = "fixed") %>%
      dplyr::select(.data$term, .data$estimate, .data$std.error)
  } else if (is.matrix(x) || is.data.frame(x)) {
    # se for each column
    se_result <- purrr::map_dbl(x, ~ std_e_helper(.x))
    # set names to return vector
    names(se_result) <- colnames(x)
  } else if (is.list(x)) {
    # compute standard error from regression coefficient and p-value
    se_result <- x$estimate / abs(stats::qnorm(x$p.value / 2))
  } else {
    # standard error for a variable
    se_result <- std_e_helper(x)
  }

  se_result
}

std_e_helper <- function(x) sqrt(var(x, na.rm = TRUE) / length(stats::na.omit(x)))


#' @importFrom stats coef setNames vcov
#' @importFrom lme4 ranef
std_merMod <- function(fit) {
  se.merMod <- list()

  # get coefficients
  cc <- stats::coef(fit)

  # get names of intercepts
  inames <- names(cc)

  # variances of fixed effects
  fixed.vars <- diag(as.matrix(stats::vcov(fit)))

  # extract variances of conditional modes
  r1 <- lme4::ranef(fit, condVar = TRUE)

  # we may have multiple random intercepts, iterate all
  for (i in seq_len(length(cc))) {
    cmode.vars <- t(apply(attr(r1[[i]], "postVar"), 3, diag))
    seVals <- sqrt(sweep(cmode.vars, 2, fixed.vars, "+"))

    # add results to return list
    se.merMod[[length(se.merMod) + 1]] <-
      stats::setNames(as.vector(seVals[1, ]), c("intercept_se", "slope_se"))
  }

  # set names of list
  names(se.merMod) <- inames

  se.merMod
}


std_e_icc <- function(x, nsim) {
  # check whether model is still in environment?
  obj.name <- attr(x, ".obj.name", exact = T)
  if (!exists(obj.name, envir = globalenv()))
    stop(sprintf("Can't find merMod-object `%s` (that was used to compute the ICC) in the environment.", obj.name), call. = F)

  # get object, see whether formulas match
  fitted.model <- globalenv()[[obj.name]]
  model.formula <- attr(x, "formula", exact = T)

  if (!identical(model.formula, formula(fitted.model)))
    stop(sprintf("merMod-object `%s` was fitted with a different formula than ICC-model", obj.name), call. = F)

  # get model family, we may have glmer
  model.family <- attr(x, "family", exact = T)

  # check for all required arguments
  if (missing(nsim) || is.null(nsim)) nsim <- 100

  # get ICC, and compute bootstrapped SE, than return both
  bstr <-
    bootstr_icc_se(model_frame(fitted.model),
                   nsim,
                   model.formula,
                   model.family)

  # now compute SE and p-values for the bootstrapped ICC
  res <- data.frame(
    model = obj.name,
    icc = as.vector(x),
    std.err = boot_se(bstr)[["std.err"]],
    p.value = boot_p(bstr)[["p.value"]]
  )

  structure(class = "se.icc.lme4", list(result = res, bootstrap_data = bstr))
}



#' @importFrom dplyr mutate
#' @importFrom lme4 lmer glmer
#' @importFrom utils txtProgressBar
#' @importFrom purrr map map_dbl
#' @importFrom rlang .data
bootstr_icc_se <- function(dd, nsim, formula, model.family) {
  # create progress bar
  pb <- utils::txtProgressBar(min = 1, max = nsim, style = 3)

  # generate bootstraps
  dummy <- dd %>%
    bootstrap(nsim) %>%
    dplyr::mutate(
      models = purrr::map(.data$strap, function(x) {
        # update progress bar
        utils::setTxtProgressBar(pb, x$resample.id)
        # check model family, then compute mixed model
        if (model.family == "gaussian")
          lme4::lmer(formula, data = x)
        else
          lme4::glmer(formula, data = x, family = model.family)
      }),
      # compute ICC(s) for each "bootstrapped" regression
      icc = purrr::map(.data$models, icc)
    )

  # we may have more than one random term in the model, so we might have
  # multiple ICC values. In this case, we need to split the multiple icc-values
  # into multiple columns, i.e. one column per ICC value
  icc_data <-
    tibble::as_tibble(matrix(unlist(purrr::map(dummy$icc, ~ .x)), nrow = nrow(dummy)))

  # close progresss bar
  close(pb)

  icc_data
}
