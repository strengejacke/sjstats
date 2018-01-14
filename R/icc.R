#' @title Intraclass-Correlation Coefficient
#' @name icc
#' @description This function calculates the intraclass-correlation
#'                (icc) - sometimes also called \emph{variance partition coefficient}
#'                (vpc) - for random intercepts of mixed effects models.
#'                Currently, \code{\link[lme4]{merMod}}, \code{\link[glmmTMB]{glmmTMB}}
#'                \code{stanreg} and \code{\link[brms]{brmsfit}}
#'                objects are supported.
#'
#' @param x Fitted mixed effects model (of class \code{merMod}, \code{glmmTMB},
#'          \code{stanreg} or \code{brmsfit}).
#' @param ... More fitted model objects, to compute multiple intraclass-correlation
#'              coefficients at once.
#'
#' @return A numeric vector with all random intercept intraclass-correlation-coefficients,
#'           or a list of numeric vectors, when more than one model were used
#'           as arguments. Furthermore, between- and within-group variances as well
#'           as random-slope variance are returned as attributes.
#'
#' @references \itemize{
#'               \item Aguinis H, Gottfredson RK, Culpepper SA. 2013. Best-Practice Recommendations for Estimating Cross-Level Interaction Effects Using Multilevel Modeling. Journal of Management 39(6): 1490–1528 (\doi{10.1177/0149206313478188})
#'               \item Aly SS, Zhao J, Li B, Jiang J. 2014. Reliability of environmental sampling culture results using the negative binomial intraclass correlation coefficient. Springerplus 14(3) (\doi{10.1186/2193-1801-3-40})
#'               \item Grace-Martion K. The Intraclass Correlation Coefficient in Mixed Models, \href{http://www.theanalysisfactor.com/the-intraclass-correlation-coefficient-in-mixed-models/}{web}
#'               \item Hox J. 2002. Multilevel analysis: techniques and applications. Mahwah, NJ: Erlbaum
#'               \item Rabe-Hesketh S, Skrondal A. 2012. Multilevel and longitudinal modeling using Stata. 3rd ed. College Station, Tex: Stata Press Publication
#'               \item Raudenbush SW, Bryk AS. 2002. Hierarchical linear models: applications and data analysis methods. 2nd ed. Thousand Oaks: Sage Publications
#'               \item Stryhn H, Sanchez J, Morley P, Booker C, Dohoo IR. 2006. Interpretation of variance parameters in multilevel Poisson regression models. Proceedings of the 11th International Symposium on Veterinary Epidemiology and Economics, 2006 Available at \url{http://www.sciquest.org.nz/node/64294}
#'               \item Wu S, Crespi CM, Wong WK. 2012. Comparison of methods for estimating the intraclass correlation coefficient for binary responses in cancer prevention cluster randomized trials. Contempory Clinical Trials 33: 869-880 (\doi{10.1016/j.cct.2012.05.004})
#'             }
#'             Further helpful online-ressources:
#'             \itemize{
#'               \item \href{http://stats.stackexchange.com/questions/18088/intraclass-correlation-icc-for-an-interaction/28100#28100}{CrossValidated (2012) \emph{Intraclass correlation (ICC) for an interaction?}}
#'               \item \href{http://stats.stackexchange.com/questions/113577/interpreting-the-random-effect-in-a-mixed-effect-model/113825#113825}{CrossValidated (2014) \emph{Interpreting the random effect in a mixed-effect model}}
#'               \item \href{http://stats.stackexchange.com/questions/67247/how-to-partition-the-variance-explained-at-group-level-and-individual-level/67356#67356}{CrossValidated (2014) \emph{how to partition the variance explained at group level and individual level}}
#'             }
#'
#'
#' @note Some notes on why the ICC is useful, based on \cite{Grace-Martin}:
#'       \itemize{
#'        \item It can help you determine whether or not a linear mixed model is even necessary. If you find that the correlation is zero, that means the observations within clusters are no more similar than observations from different clusters. Go ahead and use a simpler analysis technique.
#'        \item It can be theoretically meaningful to understand how much of the overall variation in the response is explained simply by clustering. For example, in a repeated measures psychological study you can tell to what extent mood is a trait (varies among people, but not within a person on different occasions) or state (varies little on average among people, but varies a lot across occasions).
#'        \item It can also be meaningful to see how the ICC (as well as the between and within cluster variances) changes as variable are added to the model.
#'       }
#'       In short, the ICC can be interpreted as \dQuote{the proportion of the variance
#'       explained by the grouping structure in the population} \cite{(Hox 2002: 15)}.
#'       \cr \cr
#'       Usually, the ICC is calculated for the null model ("unconditional model").
#'       However, according to \cite{Raudenbush and Bryk (2002)} or
#'       \cite{Rabe-Hesketh and Skrondal (2012)} it is also feasible to compute the ICC
#'       for full models with covariates ("conditional models") and compare how
#'       much a level-2 variable explains the portion of variation in the grouping
#'       structure (random intercept).
#'       \cr \cr
#'       \strong{Caution:} For three-level-models, depending on the nested structure
#'       of the model, the ICC only reports the proportion of variance explained
#'       for each grouping level. However, the proportion of variance for specific
#'       levels related to each other (e.g., similarity of level-1-units within
#'       level-2-units or level-2-units within level-3-units) must be computed
#'       manually. Use \code{\link{get_re_var}} to get the between-group-variances
#'       and residual variance of the model, and calculate the ICC for the various level
#'       correlations.
#'       \cr \cr
#'       For example, for the ICC between level 1 and 2: \cr
#'       \code{sum(get_re_var(fit)) / (sum(get_re_var(fit)) + get_re_var(fit, "sigma_2"))}
#'       \cr \cr
#'       or for the ICC between level 2 and 3: \cr
#'       \code{get_re_var(fit)[2] / sum(get_re_var(fit))}
#'
#' @details The ICC is calculated by dividing the between-group-variance (random
#'          intercept variance) by the total variance (i.e. sum of between-group-variance
#'          and within-group (residual) variance). \cr \cr
#'       The calculation of the ICC for generalized linear mixed models with binary outcome is based on
#'       \cite{Wu et al. (2012)}. For Poisson multilevel models, please refere to \cite{Stryhn et al. (2006)}.
#'       \cite{Aly et al. (2014)} describe computation of ICC for negative binomial models.
#'       \cr \cr
#'       There is a \code{print}-method that prints the variance parameters using
#'       the \code{comp}-argument set to \code{"var"}: \code{print(x, comp = "var")}
#'       (see 'Examples'). The \code{\link{re_var}}-function is a convenient wrapper.
#'       \cr \cr
#'       The random effect variances indicate the between- and within-group
#'         variances as well as random-slope variance and random-slope-intercept
#'         correlation. The components are denoted as following:
#'         \itemize{
#'          \item Within-group (residual) variance: sigma_2
#'          \item Between-group-variance: tau.00 (variation between individual intercepts and average intercept)
#'          \item Random-slope-variance: tau.11 (variation between individual slopes and average slope)
#'          \item Random-Intercept-Slope-covariance: tau.01
#'          \item Random-Intercept-Slope-correlation: rho.01
#'         }
#'
#' @seealso \code{\link{re_var}}
#'
#' @examples
#' library(lme4)
#' fit0 <- lmer(Reaction ~ 1 + (1 | Subject), sleepstudy)
#' icc(fit0)
#'
#' fit1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' icc(fit1)
#'
#' sleepstudy$mygrp <- sample(1:45, size = 180, replace = TRUE)
#' fit2 <- lmer(Reaction ~ Days + (1 | mygrp) + (Days | Subject), sleepstudy)
#' icc(fit2)
#'
#' # return icc for all models at once
#' icc(fit0, fit1, fit2)
#'
#' icc1 <- icc(fit1)
#' icc2 <- icc(fit2)
#'
#' print(icc1, comp = "var")
#' print(icc2, comp = "var")
#'
#'
#' @importFrom purrr map2
#' @export
icc <- function(x, ...) {
  # return value
  icc_ <- icc.lme4(x, deparse(substitute(x)))

  # check if we have multiple parameters
  if (nargs() > 1) {
    # evaluate dots
    dots <- match.call(expand.dots = FALSE)$`...`
    # get paramater names
    dot.names <- dot_names(dots)

    # get input list
    params_ <- list(...)
    icc_ <- list(icc_)

    for (i in seq_len(length(params_))) {
      icc_[[length(icc_) + 1]] <- icc.lme4(params_[[i]], dot.names[i])
    }

    names(icc_) <- NULL
  }

  icc_
}

# icc <- function(...) {
#   # evaluate dots
#   dots <- match.call(expand.dots = FALSE)$`...`
#   # get paramater names
#   dot.names <- dot_names(dots)
#
#   icc_ <- purrr::map2(list(...), dot.names, ~ icc.lme4(.x, .y))
#   names(icc_) <- NULL
#
#   if (length(icc_) == 1)
#     icc_[[1]]
#   else
#     icc_
# }


#' @importFrom lme4 VarCorr fixef getME
#' @importFrom glmmTMB VarCorr fixef getME
#' @importFrom stats family formula
#' @importFrom purrr map map_dbl map_lgl
#' @importFrom sjmisc str_contains
icc.lme4 <- function(fit, obj.name) {
  # check object class
  if (is_merMod(fit) || inherits(fit, c("glmmTMB", "brmsfit"))) {

    if (inherits(fit, "brmsfit") && !requireNamespace("brms", quietly = TRUE))
      stop("Please install and load package `brms` first.", call. = F)

    # get family
    fitfam <- stats::family(fit)$family


    # is neg. binomial?

    is_negbin <-
      sjmisc::str_contains(
        fitfam,
        c("Negative Binomial", "nbinom"),
        ignore.case = TRUE,
        logic = "OR"
      )


    # is logistic?

    is_logistic <-
      inherits(fit, c("glmerMod", "glmmTMB", "brmsfit")) &&
      fitfam %in% c("bernoulli", "binomial")


    # random effects variances
    # for details on tau and sigma, see
    # Aguinis H, Gottfredson RK, Culpepper SA2013. Best-Practice Recommendations
    # for Estimating Cross-Level Interaction Effects Using Multilevel Modeling.
    # Journal of Management 39(6): 1490–1528. doi:10.1177/0149206313478188.

    if (inherits(fit, "glmmTMB")) {
      reva <- glmmTMB::VarCorr(fit)[[1]]
    } else if (inherits(fit, "brmsfit")) {
      reva <- brms::VarCorr(fit, old = TRUE)
    } else
      reva <- lme4::VarCorr(fit)


    # for brmsfit-objects, remove "RESIDUAL" element from list
    # and save in separate object
    if (inherits(fit, "brmsfit")) {
      reva.resid <- reva[names(reva) == "RESIDUAL"]
      reva <- reva[!(names(reva) == "RESIDUAL")]
    }


    # retrieve only intercepts

    if (inherits(fit, "brmsfit"))
      vars <- purrr::map(reva, ~ .x$cov$mean[1])
    else
      vars <- purrr::map(reva, ~ .x[1])


    # random intercept-variances, i.e.
    # between-subject-variance (tau 00)

    tau.00 <- purrr::map_dbl(vars, ~ .x)


    # random slope-variances (tau 11)
    if (inherits(fit, "brmsfit"))
      tau.11 <- unlist(lapply(reva, function(x) diag(x$cov$mean)[-1]))
    else
      tau.11 <- unlist(lapply(reva, function(x) diag(x)[-1]))


    # get residual standard deviation sigma
    if (inherits(fit, "brmsfit"))
      sig <- reva.resid$RESIDUAL$sd[1]
    else
      sig <- attr(reva, "sc")


    # set default, if no residual variance is available

    if (is.null(sig)) {
      if (is_logistic)
        sig <- sqrt((pi ^ 2) / 3)
      else
        sig <- 1
    }


    # residual variances, i.e.
    # within-cluster-variance (sigma^2)

    if (is_logistic) {
      # for logistic models, we use pi / 3
      resid_var <- (pi ^ 2) / 3
    } else if (inherits(fit, "glmerMod") && is_negbin) {
      # for negative binomial models, we use 0
      resid_var <- 1
    } else {
      # for linear and poisson models, we have a clear residual variance
      resid_var <- sig ^ 2
    }


    # total variance, sum of random intercept and residual variances
    total_var <- sum(purrr::map_dbl(vars, ~ sum(.x)), resid_var)


    # check whether we have negative binomial

    if (is_negbin) {
      if (is_merMod(fit)) {
        # for negative binomial models, we also need the intercept...
        beta <- as.numeric(lme4::fixef(fit)["(Intercept)"])
        # ... and the theta value to compute the ICC
        r <- lme4::getME(fit, "glmer.nb.theta")
      } else {
        # for negative binomial models, we also need the intercept...
        beta <- as.numeric(glmmTMB::fixef(fit)[[1]]["(Intercept)"])
        # ... and the theta value to compute the ICC
        r <- sig
      }


      # make formula more readable

      numerator <- (exp(tau.00) - 1)
      denominator <- ((exp(total_var) - 1) + (exp(total_var) / r) + exp(-beta - (total_var / 2)))

      ri.icc <- numerator / denominator
    } else {
      # random intercept icc
      ri.icc <- tau.00 / total_var
    }


    # get random slope random intercept correlations
    # do we have any rnd slopes?

    if (inherits(fit, "brmsfit"))
      has_rnd_slope <- purrr::map_lgl(reva, ~ dim(.x$cor$mean)[1] > 1)
    else
      has_rnd_slope <- purrr::map_lgl(reva, ~ dim(attr(.x, "correlation"))[1] > 1)

    tau.01 <- rho.01 <- NULL


    # get rnd slopes

    if (any(has_rnd_slope)) {

      rnd_slope <- reva[has_rnd_slope]

      if (inherits(fit, "brmsfit")) {
        # get slope-intercept-correlations
        rho.01 <- purrr::map_dbl(rnd_slope, ~ .x$cor$mean[1, 2])
        # get standard deviations, multiplied
        std_ <- purrr::map_dbl(rnd_slope, ~ prod(.x$sd))
      } else {
        # get slope-intercept-correlations
        rho.01 <- purrr::map_dbl(rnd_slope, ~ attr(.x, "correlation")[1, 2])
        # get standard deviations, multiplied
        std_ <- purrr::map_dbl(rnd_slope, ~ prod(attr(.x, "stddev")))
      }

      # bind to matrix
      tau.01 <- apply(cbind(rho.01, std_), MARGIN = 1, FUN = prod)
    }

    # name values
    names(ri.icc) <- names(reva)


    if (inherits(fit, c("glmerMod", "glmmTMB")))
      mt <- "Generalized linear mixed model"
    else if (inherits(fit, "brmsfit"))
      mt <- "Bayesian mixed model"
    else
      mt <- "Linear mixed model"

    # add attributes, for print method
    class(ri.icc) <- c("icc.lme4", class(ri.icc))
    attr(ri.icc, "family") <- stats::family(fit)$family
    attr(ri.icc, "link") <- stats::family(fit)$link
    attr(ri.icc, "formula") <- stats::formula(fit)
    attr(ri.icc, "model") <- mt
    attr(ri.icc, "tau.00") <- tau.00
    attr(ri.icc, "tau.01") <- tau.01
    attr(ri.icc, "rho.01") <- rho.01
    attr(ri.icc, "tau.11") <- tau.11
    attr(ri.icc, "sigma_2") <- resid_var


    # finally, save name of fitted model object. May be needed for
    # the 'se()' function, which accesses the global environment

    attr(ri.icc, ".obj.name") <- obj.name

    # return results
    return(ri.icc)
  } else {
    warning("Function `icc` currently only supports `merMod` (package `lme4`), `glmmTMB` (package `glmmTMB`) or `brmsfit` (package brms) objects.", call. = TRUE)
  }
}


#' @title Random effect variances
#' @name re_var
#' @description These functions extracts random effect variances as well as
#'                random-intercept-slope-correlation of mixed effects models.
#'                Currently, \code{\link[lme4]{merMod}}, \code{\link[glmmTMB]{glmmTMB}}
#'                \code{stanreg} and \code{\link[brms]{brmsfit}}
#'                objects are supported.
#'
#' @param x Fitted mixed effects model (of class \code{merMod}, \code{glmmTMB},
#'          \code{stanreg} or \code{brmsfit}). \code{get_re_var()} also accepts
#'           an object of class \code{icc.lme4}, as returned by the
#'           \code{\link{icc}} function.
#' @param comp Name of the variance component to be returned. See 'Details'.
#'
#' @return \code{get_re_var()} returns the value of the requested variance component,
#'           \code{re_var()} returns all random effects variances.
#'
#' @references Aguinis H, Gottfredson RK, Culpepper SA. 2013. Best-Practice Recommendations for Estimating Cross-Level Interaction Effects Using Multilevel Modeling. Journal of Management 39(6): 1490–1528 (\doi{10.1177/0149206313478188})
#'
#' @details The random effect variances indicate the between- and within-group
#'         variances as well as random-slope variance and random-slope-intercept
#'         correlation. Use following values for \code{comp} to get the particular
#'         variance component:
#'         \describe{
#'          \item{\code{"sigma_2"}}{Within-group (residual) variance}
#'          \item{\code{"tau.00"}}{Between-group-variance (variation between individual intercepts and average intercept)}
#'          \item{\code{"tau.11"}}{Random-slope-variance (variation between individual slopes and average slope)}
#'          \item{\code{"tau.01"}}{Random-Intercept-Slope-covariance}
#'          \item{\code{"rho.01"}}{Random-Intercept-Slope-correlation}
#'         }
#'         The within-group-variance is affected by factors at level one, i.e.
#'         by the lower-level direct effects. Level two factors (i.e. cross-level
#'         direct effects) affect the between-group-variance. Cross-level
#'         interaction effects are group-level factors that explain the
#'         variance in random slopes (Aguinis et al. 2013).
#'
#' @seealso \code{\link{icc}}
#'
#' @examples
#' library(lme4)
#' fit1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#'
#' # all random effect variance components
#' re_var(fit1)
#'
#' # just the rand. slope-intercept covariance
#' get_re_var(fit1, "tau.01")
#'
#' sleepstudy$mygrp <- sample(1:45, size = 180, replace = TRUE)
#' fit2 <- lmer(Reaction ~ Days + (1 | mygrp) + (Days | Subject), sleepstudy)
#' re_var(fit2)
#'
#' @importFrom stats family
#' @importFrom purrr map map2 flatten_dbl flatten_chr
#' @importFrom sjmisc trim
#' @export
re_var <- function(x) {
  # iterate all attributes and return them as vector
  rv <- c("sigma_2", "tau.00", "tau.11", "tau.01", "rho.01")

  rv_ <- purrr::map(rv, ~ attr(icc(x), .x, exact = TRUE))
  rn <- purrr::map2(1:length(rv_), rv, ~ sjmisc::trim(paste(names(rv_[[.x]]), .y, sep = "_")))
  rv_ <- purrr::flatten_dbl(rv_)

  names(rv_) <- purrr::flatten_chr(rn)[1:length(rv_)]

  class(rv_) <- c("sj_revar", class(rv_))

  rv_
}


#' @rdname re_var
#' @export
get_re_var <- function(x, comp = c("tau.00", "tau.01", "tau.11", "rho.01", "sigma_2")) {
  # check if we have a valid object
  if (!inherits(x, "icc.lme4") && !is_merMod(x) && !inherits(x, c("glmmTMB", "brmsfit"))) {
    stop("`x` must either be an object returned by the `icc` function, or a merMod-, glmmTMB- or brmsfit-object.", call. = F)
  }

  # check arguments
  comp <- match.arg(comp)

  # do we have a merMod object? If yes, get ICC and var components
  if (is_merMod(x) || inherits(x, c("glmmTMB", "brmsfit"))) x <- icc(x)

  # return results
  attr(x, comp, exact = TRUE)
}
