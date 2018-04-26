#' @title Intraclass-Correlation Coefficient
#' @name icc
#' @description This function calculates the intraclass-correlation
#'    (icc) - sometimes also called \emph{variance partition coefficient}
#'    (vpc) - for random intercepts of mixed effects models. Currently,
#'    \code{\link[lme4]{merMod}}, \code{\link[glmmTMB]{glmmTMB}},
#'    \code{stanreg} and \code{\link[brms]{brmsfit}} objects are supported.
#'
#' @param x Fitted mixed effects model (of class \code{merMod}, \code{glmmTMB},
#'    \code{stanreg} or \code{brmsfit}).
#' @param ... More fitted model objects, to compute multiple intraclass-correlation
#'    coefficients at once.
#' @param posterior Logical, if \code{TRUE} and \code{x} is a \code{brmsfit}
#'    object, ICC values are computed for each sample of the posterior
#'    distribution. In this case, a data frame is returned with the same
#'    number of rows as samples in \code{x}, with one column per random
#'    effect ICC.
#'
#' @return If \code{posterior = FALSE} (the default), a numeric vector with all
#'    random intercept intraclass-correlation-coefficients, or a list of
#'    numeric vectors, when more than one model were used as arguments.
#'    Furthermore, between- and within-group variances as well as random-slope
#'    variance are returned as attributes.
#'    \cr \cr
#'    If \code{posterior = TRUE}, \code{icc()} returns a data frame with ICC
#'    and variance components for each sample of the posterior distribution.
#'
#' @references \itemize{
#'               \item Aguinis H, Gottfredson RK, Culpepper SA. 2013. Best-Practice Recommendations for Estimating Cross-Level Interaction Effects Using Multilevel Modeling. Journal of Management 39(6): 1490–1528 (\doi{10.1177/0149206313478188})
#'               \item Aly SS, Zhao J, Li B, Jiang J. 2014. Reliability of environmental sampling culture results using the negative binomial intraclass correlation coefficient. Springerplus 14(3) (\doi{10.1186/2193-1801-3-40})
#'               \item Goldstein H, Browne W, Rasbash J. 2010. Partitioning Variation in Multilevel Models. Understanding Statistics, 1:4, 223-231 (\doi{10.1207/S15328031US0104_02})
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
#'       \cite{Wu et al. (2012)}. For Poisson multilevel models, please refer to \cite{Stryhn et al. (2006)}.
#'       \cite{Aly et al. (2014)} describe computation of ICC for negative binomial models.
#'       \cr \cr
#'       \strong{Caution:} For models with random slopes and random intercepts,
#'       the ICC would differ at each unit of the predictors. Hence, the ICC for these
#'       kind of models cannot be understood simply as proportion of variance
#'       (see \cite{Goldstein et al. 2010}). For convenience reasons, as the
#'       \code{icc()} function also extracts the different random effects
#'       variances, the ICC for random-slope-intercept-models is reported
#'       nonetheless, but it is usually no meaningful summary of the
#'       proportion of variances.
#'       \cr \cr
#'       If \code{posterior = FALSE}, there is a \code{print()}-method that prints
#'       the variance parameters using the \code{comp}-argument set to \code{"var"}:
#'       \code{print(x, comp = "var")} (see 'Examples'). The
#'       \code{\link{re_var}}-function is a convenient wrapper. If
#'       \code{posterior = TRUE}, the \code{print()}-method accepts the arguments
#'       \code{prob} and \code{digits}, which indicate the probability of the
#'       uncertainty interval for the ICC and variance components, and the digits
#'       in the output (see also 'Examples').
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
#' # note: ICC for random-slope-intercept model usually not
#' # meaningful - see 'Note'.
#' fit1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' icc(fit1)
#'
#' sleepstudy$mygrp <- sample(1:45, size = 180, replace = TRUE)
#' fit2 <- lmer(Reaction ~ Days + (1 | mygrp) + (1 | Subject), sleepstudy)
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
#' \dontrun{
#' # compute ICC for Bayesian mixed model, with an ICC for each
#' # sample of the posterior. The print()-method then shows
#' # the median ICC as well as 89% HDI for the ICC.
#' # Change interval with print-method:
#' # print(icc(m, posterior = TRUE), prob = .5)
#'
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   library(dplyr)
#'   sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
#'   sleepstudy <- sleepstudy %>%
#'     group_by(mygrp) %>%
#'     mutate(mysubgrp = sample(1:30, size = n(), replace = TRUE))
#'   m <- brms::brm(
#'     Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
#'     data = sleepstudy
#'   )
#'
#'   # by default, 89% interval
#'   icc(m, posterior = TRUE)
#'
#'   # show 50% interval
#'   print(icc(m, posterior = TRUE), prob = .5, digits = 3)
#' }}
#'
#' @importFrom purrr map2
#' @export
icc <- function(x, ..., posterior = FALSE) {

  if (isTRUE(posterior) && !inherits(x, "brmsfit")) {
    warning("ICC from posterior samples only possible for `brmsfit`-objects.")
    posterior <- FALSE
  }

  # return value
  if (posterior) {
    return(icc.posterior(x, deparse(substitute(x))))
  }

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
#' @importFrom tibble has_name
icc.lme4 <- function(fit, obj.name) {
  # check object class
  if (is_merMod(fit) || inherits(fit, c("glmmTMB", "brmsfit"))) {

    if (inherits(fit, "brmsfit") && !requireNamespace("brms", quietly = TRUE))
      stop("Please install and load package `brms` first.", call. = F)

    # get family
    fitfam <- model_family(fit)


    # random effects variances
    # for details on tau and sigma, see
    # Aguinis H, Gottfredson RK, Culpepper SA2013. Best-Practice Recommendations
    # for Estimating Cross-Level Interaction Effects Using Multilevel Modeling.
    # Journal of Management 39(6): 1490–1528. doi:10.1177/0149206313478188.

    if (inherits(fit, "glmmTMB")) {
      reva <- glmmTMB::VarCorr(fit)[[1]]
    } else if (inherits(fit, "brmsfit")) {
      reva <- brms::VarCorr(fit)
    } else
      reva <- lme4::VarCorr(fit)


    # for brmsfit-objects, remove "residual__" element from list
    # and save in separate object
    if (inherits(fit, "brmsfit")) {
      reva.resid <- reva[names(reva) == "residual__"]
      reva <- reva[!(names(reva) == "residual__")]
    }


    # retrieve only intercepts

    if (inherits(fit, "brmsfit"))
      vars <- purrr::map(reva, ~ .x$sd[1] ^ 2)
    else
      vars <- purrr::map(reva, ~ .x[1])


    # random intercept-variances, i.e.
    # between-subject-variance (tau 00)

    tau.00 <- purrr::map_dbl(vars, ~ .x)


    # random slope-variances (tau 11)
    if (inherits(fit, "brmsfit"))
      tau.11 <- unlist(lapply(reva, function(x) diag(x$cov[, 1, ])[-1]))
    else
      tau.11 <- unlist(lapply(reva, function(x) diag(x)[-1]))


    # get residual standard deviation sigma
    if (inherits(fit, "brmsfit"))
      sig <- reva.resid[["residual__"]]$sd[1]
    else
      sig <- attr(reva, "sc")


    # set default, if no residual variance is available

    if (is.null(sig)) {
      if (fitfam$is_bin)
        sig <- sqrt((pi ^ 2) / 3)
      else
        sig <- 1
    }


    # residual variances, i.e.
    # within-cluster-variance (sigma^2)

    if (fitfam$is_bin) {
      # for logistic models, we use pi / 3
      resid_var <- (pi ^ 2) / 3
    } else if (inherits(fit, "glmerMod") && fitfam$is_negbin) {
      # for negative binomial models, we use 1
      resid_var <- 1
    } else {
      # for linear and poisson models, we have a clear residual variance
      resid_var <- sig ^ 2
    }


    # total variance, sum of random intercept and residual variances
    total_var <- sum(purrr::map_dbl(vars, ~ sum(.x)), resid_var)


    # check whether we have negative binomial

    if (fitfam$is_negbin) {
      if (is_merMod(fit)) {
        # for negative binomial models, we also need the intercept...
        beta <- as.numeric(lme4::fixef(fit)["(Intercept)"])
        # ... and the theta value to compute the ICC
        r <- lme4::getME(fit, "glmer.nb.theta")
      } else if (inherits(fit, "brms")) {
        # for negative binomial models, we also need the intercept...
        beta <- as.numeric(brms::fixef(fit)[[1]])
        # ... and the theta value to compute the ICC
        r <- sig
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
      has_rnd_slope <- purrr::map_lgl(reva, ~ tibble::has_name(.x, "cor"))
    else
      has_rnd_slope <- purrr::map_lgl(reva, ~ dim(attr(.x, "correlation"))[1] > 1)

    tau.01 <- rho.01 <- NULL


    # get rnd slopes

    if (any(has_rnd_slope)) {

      rnd_slope <- reva[has_rnd_slope]

      if (inherits(fit, "brmsfit")) {
        # get slope-intercept-correlations
        rho.01 <- purrr::map_dbl(rnd_slope, ~ .x$cor[1, 1, 2])
        # get standard deviations, multiplied
        std_ <- purrr::map_dbl(rnd_slope, ~ prod(.x$sd[, 1]))
      } else {
        # get slope-intercept-correlations
        rho.01 <- purrr::map_dbl(rnd_slope, ~ attr(.x, "correlation")[1, 2])
        # get standard deviations, multiplied
        std_ <- purrr::map_dbl(rnd_slope, ~ prod(attr(.x, "stddev")))
      }

      # bind to matrix
      tau.01 <- apply(cbind(rho.01, std_), MARGIN = 1, FUN = prod)

      message("Caution! ICC for random-slope-intercept models usually not meaningful. See 'Note' in `?icc`.")
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
    attr(ri.icc, "rnd.slope.model") <- any(has_rnd_slope)


    # finally, save name of fitted model object. May be needed for
    # the 'se()' function, which accesses the global environment

    attr(ri.icc, ".obj.name") <- obj.name

    # return results
    return(ri.icc)
  } else {
    warning("`icc()` does not support this model-object.", call. = TRUE)
  }
}


#' @importFrom purrr map_df map_if
icc.posterior <- function(fit, obj.name) {
  if (inherits(fit, "brmsfit") && !requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.", call. = F)

  # get family
  fitfam <- model_family(fit)

  # get random effect variances for each sample of posterior
  reva <- brms::VarCorr(fit, summary = FALSE)

  # remove "residual__" element from list
  # and save in separate object
  reva.resid <- reva[names(reva) == "residual__"]
  reva <- reva[!(names(reva) == "residual__")]


  # retrieve only intercepts
  vars <- purrr::map(reva, ~ .x$sd[, 1] ^ 2)

  # random intercept-variances, i.e.
  # between-subject-variance (tau 00)
  tau.00 <- purrr::map(vars, ~ .x)

  # random slope-variances (tau 11)
  tau.11 <- purrr::map(reva, ~ .x$cov[, 2, 2])

  # get residual standard deviation sigma
  sig <- reva.resid[["residual__"]]$sd[, 1]

  # set default, if no residual variance is available
  if (is.null(sig)) {
    if (fitfam$is_bin)
      sig <- sqrt((pi ^ 2) / 3)
    else
      sig <- 1
  }


  # residual variances, i.e.
  # within-cluster-variance (sigma^2)

  if (fitfam$is_bin) {
    # for logistic models, we use pi / 3
    resid_var <- (pi ^ 2) / 3
  } else {
    # for linear and poisson models, we have a clear residual variance
    resid_var <- sig ^ 2
  }


  # total variance, sum of random intercept and residual variances
  total_var <- apply(as.data.frame(vars), MARGIN = 1, FUN = sum) + resid_var

  # make sure residual variance has same length as other components
  # if not, just repeat the current value to match number of samples
  if (length(resid_var) == 1) resid_var <- rep(resid_var, length(total_var))

  # check whether we have negative binomial

  if (fitfam$is_negbin) {

    # for negative binomial models, we also need the intercept...
    beta <- as.numeric(brms::fixef(fit)[[1]])
    # ... and the theta value to compute the ICC
    r <- sig

    # make formula more readable

    numerator <- purrr::map(tau.00, ~ exp(.x) - 1)
    denominator <- ((exp(total_var) - 1) + (exp(total_var) / r) + exp(-beta - (total_var / 2)))

    ri.icc <- purrr::map(numerator, ~ .x / denominator)
  } else {
    # random intercept icc
    ri.icc <- purrr::map(tau.00, ~ .x / total_var)
  }

  tau.11 <- purrr::map_if(tau.11, is.null, ~ rep(NA, length(resid_var)))

  names(ri.icc) <- sprintf("icc_%s", names(ri.icc))
  names(tau.00) <- sprintf("tau.00_%s", names(tau.00))
  names(tau.11) <- sprintf("tau.11_%s", names(tau.11))

  icc_ <- dplyr::bind_cols(ri.icc, tau.00, tau.11, data.frame(resid_var = resid_var))

  attr(icc_, "family") <- stats::family(fit)$family
  attr(icc_, "link") <- stats::family(fit)$link
  attr(icc_, "formula") <- stats::formula(fit)
  attr(icc_, "model") <- "Bayesian mixed model"
  attr(icc_, "tau.00") <- tau.00
  attr(icc_, "tau.11") <- tau.11
  attr(icc_, "sigma_2") <- resid_var

  class(icc_) <- c("icc.posterior", class(icc_))

  # return results
  icc_
}


#' @title Random effect variances
#' @name re_var
#' @description These functions extracts random effect variances as well as
#'                random-intercept-slope-correlation of mixed effects models.
#'                Currently, \code{\link[lme4]{merMod}}, \code{\link[glmmTMB]{glmmTMB}},
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

  # compute icc
  icc_ <- suppressMessages(icc(x))

  rv_ <- purrr::map(rv, ~ attr(icc_, .x, exact = TRUE))
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
  if (is_merMod(x) || inherits(x, c("glmmTMB", "brmsfit"))) x <- suppressMessages(icc(x))

  # return results
  attr(x, comp, exact = TRUE)
}
