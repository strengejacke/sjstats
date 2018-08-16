#' @title Standard Error for variables or coefficients
#' @name se
#' @description Compute standard error for a variable, for all variables
#'     of a data frame, for joint random and fixed effects
#'     coefficients of (non-/linear) mixed models, the adjusted
#'     standard errors for generalized linear (mixed) models, or
#'     for intraclass correlation coefficients (ICC).
#'
#' @param x (Numeric) vector, a data frame, an \code{lm}, \code{glm},
#'    \code{merMod} (\pkg{lme4}), or \code{stanreg} model object, an ICC object
#'    (as obtained by the \code{\link{icc}}-function), a \code{table} or
#'    \code{xtabs} object, or a list with estimate and p-value. For the latter
#'    case, the list must contain elements named \code{estimate} and
#'    \code{p.value} (see 'Examples' and 'Details').
#' @param nsim Numeric, the number of simulations for calculating the
#'          standard error for intraclass correlation coefficients, as
#'          obtained by the \code{\link{icc}}-function.
#' @param ... Currently not used.
#'
#' @return The standard error of \code{x}.
#'
#' @note Computation of standard errors for coefficients of mixed models
#'    is based \href{http://stackoverflow.com/questions/26198958/extracting-coefficients-and-their-standard-error-from-lme}{on this code}.
#'    Standard errors for generalized linear (mixed) models, if
#'    \code{type = "re"}, are approximations based on the delta
#'    method (Oehlert 1992).
#'    \cr \cr
#'    A remark on standard errors:
#'    \dQuote{Standard error represents variation in the point estimate, but
#'    confidence interval has usual Bayesian interpretation only with flat prior.}
#'    (Gelman 2017)
#'
#' @details \strong{Standard error for variables}
#'   \cr \cr
#'   For variables and data frames, the standard error is the square root of the
#'   variance divided by the number of observations (length of vector).
#'   \cr \cr
#'   \strong{Standard error for mixed models}
#'   \cr \cr
#'   For linear mixed models, and generalized linear mixed models, this
#'   function computes the standard errors for joint (sums of) random and fixed
#'   effects coefficients (unlike \code{\link[arm]{se.coef}}, which returns the
#'   standard error for fixed and random effects separately). Hence, \code{se()}
#'   returns the appropriate standard errors for \code{\link[lme4]{coef.merMod}}.
#'   \cr \cr
#'   \strong{Standard error for generalized linear models}
#'   \cr \cr
#'   For generalized linear models, approximated standard errors, using the delta
#'   method for transformed regression parameters are returned (Oehlert 1992).
#'   \cr \cr
#'   \strong{Standard error for Intraclass Correlation Coefficient (ICC)}
#'   \cr \cr
#'   The standard error for the \code{\link{icc}} is based on bootstrapping,
#'   thus, the \code{nsim}-argument is required. See 'Examples'.
#'   \cr \cr
#'   \strong{Standard error for proportions and mean value}
#'   \cr \cr
#'   To compute the standard error for relative frequencies (i.e. proportions, or
#'   mean value if \code{x} has only two categories), this vector must be supplied
#'   as table, e.g. \code{se(table(iris$Species))}. \code{se()} than computes the
#'   relative frequencies (proportions) for each value and the related standard
#'   error for each value. This might be useful to add standard errors or confidence
#'   intervals to descriptive statistics. If standard errors for weighted variables
#'   are required, use \code{xtabs()}, e.g. \code{se(xtabs(weights ~ variable))}.
#'   \cr \cr
#'   \strong{Standard error for regression coefficient and p-value}
#'   \cr \cr
#'   \code{se()} also returns the standard error of an estimate (regression
#'   coefficient) and p-value, assuming a normal distribution to compute
#'   the z-score from the p-value (formula in short: \code{b / qnorm(p / 2)}).
#'   See 'Examples'.
#'
#' @references Oehlert GW. 1992. A note on the delta method. American Statistician 46(1).
#'             \cr \cr
#'             Gelman A 2017. How to interpret confidence intervals? \url{http://andrewgelman.com/2017/03/04/interpret-confidence-intervals/}
#'
#' @examples
#' library(lme4)
#' library(sjmisc)
#'
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
#' fit <- glm(
#'   services ~ neg_c_7 + c161sex + e42dep,
#'   data = efc,
#'   family = binomial(link = "logit")
#' )
#' se(fit)
#'
#' # compute odds-ratio adjusted standard errors for generalized
#' # linear mixed model, also based on delta method
#'
#' # create binary response
#' sleepstudy$Reaction.dicho <- dicho(sleepstudy$Reaction, dich.by = "median")
#' fit <- glmer(
#'   Reaction.dicho ~ Days + (Days | Subject),
#'   data = sleepstudy,
#'   family = binomial("logit")
#' )
#' se(fit)
#'
#' # compute standard error for proportions
#' efc$e42dep <- to_label(efc$e42dep)
#' se(table(efc$e42dep))
#'
#' # including weights
#' efc$weights <- rnorm(nrow(efc), 1, .25)
#' se(xtabs(efc$weights ~ efc$e42dep))
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
#' library(dplyr)
#' library(purrr)
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
#' @export
se <- function(x, ...) {
  UseMethod("se")
}


#' @export
se.default <- function(x, ...) {
  std_e_helper(x)
}


#' @importFrom stats qnorm
#' @export
se.list <- function(x, ...) {
  # compute standard error from regression coefficient and p-value
  x$estimate / abs(stats::qnorm(x$p.value / 2))
}


#' @export
se.table <- function(x, ...) {
  se_tab(x)
}

#' @export
se.xtabs <- function(x, ...) {
  se_tab(x)
}

se_tab <- function(x, ...) {
  # compute standard error of proportions
  if (length(dim(x)) == 1) {
    total.n <- as.vector(sum(x))
    rel.frq <- as.vector(x) / total.n
    sev <- data_frame(
      value = names(x),
      proportion = rel.frq,
      std.error = suppressWarnings(sqrt(rel.frq * (1 - rel.frq) / total.n))
    )
  } else {
    sev <- NA
  }

  sev
}


#' @importFrom purrr map_dbl map_lgl
#' @export
se.data.frame <- function(x, ...) {
  # se for each column
  cols <- purrr::map_lgl(x, is.numeric)
  se_result <- purrr::map_dbl(x[, cols], ~ std_e_helper(.x))
  # set names to return vector
  names(se_result) <- colnames(x)[cols]
  se_result
}


#' @importFrom broom tidy
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
se.lm <- function(x, ...) {
  x %>%
    broom::tidy(effects = "fixed") %>%
    dplyr::select(.data$term, .data$estimate, .data$std.error)
}


#' @importFrom stats qnorm vcov
#' @importFrom broom tidy
#' @importFrom dplyr mutate select
#' @importFrom rlang .data
#' @export
se.glm <- function(x, ...) {
  # for glm, we want to exponentiate coefficients to get odds ratios, however
  # 'exponentiate'-argument currently not works for lme4-tidiers
  # so we need to do this manually for glmer's

  tm <- broom::tidy(x, effects = "fixed")
  tm$estimate <- exp(tm$estimate)

  tm %>%
    dplyr::mutate(std.error = sqrt(.data$estimate^2 * diag(as.matrix(stats::vcov(x))))) %>%
    dplyr::select(.data$term, .data$estimate, .data$std.error)
}


#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
se.svymle <- function(x, ...) {
  x %>%
    tidy_svyglm.nb() %>%
    dplyr::select(.data$term, .data$estimate, .data$std.error)
}


#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
se.svyglm.nb <- function(x, ...) {
  x %>%
    tidy_svyglm.nb() %>%
    dplyr::select(.data$term, .data$estimate, .data$std.error)
}


#' @rdname se
#' @export
se.icc.lme4 <- function(x, nsim = 100, ...) {
  std_e_icc(x, nsim)
}


#' @export
se.glmerMod <- function(x, ...) {
  std_merMod(x)
}


#' @export
se.lmerMod <- function(x, ...) {
  std_merMod(x)
}


#' @export
se.nlmerMod <- function(x, ...) {
  std_merMod(x)
}


#' @export
se.merModLmerTest <- function(x, ...) {
  std_merMod(x)
}


#' @importFrom broom tidy
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
se.stanreg <- function(x, ...) {
  x %>%
    broom::tidy() %>%
    dplyr::select(.data$term, .data$estimate, .data$std.error)
}


#' @importFrom broom tidy
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
se.stanfit <- function(x, ...) {
  x %>%
    broom::tidy() %>%
    dplyr::select(.data$term, .data$estimate, .data$std.error)
}


#' @importFrom stats var na.omit
std_e_helper <- function(x) sqrt(stats::var(x, na.rm = TRUE) / length(stats::na.omit(x)))


#' @importFrom purrr map
#' @importFrom stats coef setNames vcov
#' @importFrom lme4 ranef
std_merMod <- function(fit) {

  # see https://stackoverflow.com/a/26206495/2094622

  # get coefficients
  cc <- stats::coef(fit)

  # get names of intercepts
  inames <- names(cc)

  # variances of fixed effects
  fixed.vars <- diag(as.matrix(stats::vcov(fit)))

  # extract variances of conditional modes
  r1 <- lme4::ranef(fit, condVar = TRUE)

  # we may have multiple random intercepts, iterate all
  se.merMod <- purrr::map(1:length(cc), function(i) {
    cmode.vars <- t(apply(attr(r1[[i]], "postVar"), 3, diag))
    seVals <- sqrt(sweep(cmode.vars, 2, fixed.vars[names(r1[[i]])], "+", check.margin = F))

    if (length(r1[[i]]) == 1) {
      seVals <- as.data.frame(t(seVals))
      stats::setNames(seVals, names(r1[[i]]))
    } else {
      seVals <- seVals[, 1:2]
      stats::setNames(as.data.frame(seVals), names(r1[[i]]))
    }
  })

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
    bootstr_icc_se(model_frame(fitted.model, fe.only = FALSE),
                   nsim,
                   model.formula,
                   model.family)

  # now compute SE and p-values for the bootstrapped ICC
  res <- data_frame(
    model = obj.name,
    icc = as.vector(x),
    std.err = boot_se(bstr)[["std.err"]],
    p.value = boot_p(bstr)[["p.value"]]
  )

  structure(class = "sj_se_icc", list(result = res, bootstrap_data = bstr))
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
    as.data.frame(matrix(unlist(purrr::map(dummy$icc, ~ .x)), nrow = nrow(dummy)))

  # close progresss bar
  close(pb)

  icc_data
}

## TODO se for poisson etc

# # for poisson family, we need a different delta method approach
# if (get_glm_family(x)$is_pois) {
#   # standard errors scaled using square root of Pearson
#   # chi-squared dispersion
#   pr <- sum(stats::residuals(x, type = "pearson")^2)
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
