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
#' @param ... Currently not used.
#' @param re.form Formula containing group-level effects to be considered in
#'   the prediction. If \code{NULL} (default), include all group-level effects.
#'   Else, for instance for nested models, name a specific group-level effect
#'   to calculate the ICC for this group-level. Only applies if \code{ppd = TRUE}.
#' @param typical Character vector, naming the function that will be used as
#'   measure of central tendency for the ICC. The default is "mean". See
#'   \link{typical_value} for options.
#' @param ppd Logical, if \code{TRUE}, variance decomposition is based on the
#'   posterior predictive distribution, which is the correct way for Bayesian
#'   non-Gaussian models. By default, \code{ppd} is set to \code{TRUE} for
#'   non-Gaussian models.If \code{adjusted = TRUE} and \code{ppd = FALSE},
#'   variance decomposition is approximated following the suggestion by
#'   \cite{Nakagawa et al. 2017} (see 'Details'), however, this is currently
#'   only implemented for Gaussian models.
#' @param adjusted Logical, if \code{TRUE}, the adjusted (and
#'   conditional) ICC is calculated, which reflects the uncertainty of all
#'   random effects (see 'Details'). For Bayesian models, if \code{ppd = TRUE},
#'   \code{adjusted} will be ignored.
#'
#' @inheritParams hdi
#'
#' @return A numeric vector with all random intercept intraclass-correlation-coefficients.
#'    Furthermore, if \code{adjusted = FALSE}, between- and within-group variances
#'    as well as random-slope variance are returned as attributes.
#'    \cr \cr
#'    For \code{stanreg} or \code{brmsfit} objects, the HDI for each statistic
#'    is also  included as attribute.
#'
#' @references \itemize{
#'    \item Aguinis H, Gottfredson RK, Culpepper SA. 2013. Best-Practice Recommendations for Estimating Cross-Level Interaction Effects Using Multilevel Modeling. Journal of Management 39(6): 1490–1528 (\doi{10.1177/0149206313478188})
#'    \item Goldstein H, Browne W, Rasbash J. 2010. Partitioning Variation in Multilevel Models. Understanding Statistics, 1:4, 223-231 (\doi{10.1207/S15328031US0104_02})
#'    \item Grace-Martion K. The Intraclass Correlation Coefficient in Mixed Models, \href{http://www.theanalysisfactor.com/the-intraclass-correlation-coefficient-in-mixed-models/}{web}
#'    \item Hox J. 2002. Multilevel analysis: techniques and applications. Mahwah, NJ: Erlbaum
#'    \item Johnson PC, O'Hara RB. 2014. Extension of Nakagawa & Schielzeth's R2GLMM to random slopes models. Methods Ecol Evol, 5: 944-946. (\doi{10.1111/2041-210X.12225})
#'    \item Nakagawa S, Johnson P, Schielzeth H (2017) The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisted and expanded. J. R. Soc. Interface 14. \doi{10.1098/rsif.2017.0213}
#'    \item Rabe-Hesketh S, Skrondal A. 2012. Multilevel and longitudinal modeling using Stata. 3rd ed. College Station, Tex: Stata Press Publication
#'    \item Raudenbush SW, Bryk AS. 2002. Hierarchical linear models: applications and data analysis methods. 2nd ed. Thousand Oaks: Sage Publications
#'    \item Wu S, Crespi CM, Wong WK. 2012. Comparison of methods for estimating the intraclass correlation coefficient for binary responses in cancer prevention cluster randomized trials. Contempory Clinical Trials 33: 869-880 (\doi{10.1016/j.cct.2012.05.004})
#'  }
#'  Further helpful online-ressources:
#'  \itemize{
#'    \item \href{http://stats.stackexchange.com/questions/18088/intraclass-correlation-icc-for-an-interaction/28100#28100}{CrossValidated (2012) \emph{Intraclass correlation (ICC) for an interaction?}}
#'    \item \href{http://stats.stackexchange.com/questions/113577/interpreting-the-random-effect-in-a-mixed-effect-model/113825#113825}{CrossValidated (2014) \emph{Interpreting the random effect in a mixed-effect model}}
#'    \item \href{http://stats.stackexchange.com/questions/67247/how-to-partition-the-variance-explained-at-group-level-and-individual-level/67356#67356}{CrossValidated (2014) \emph{how to partition the variance explained at group level and individual level}}
#'  }
#'
#'
#' @note Some notes on why the ICC is useful, based on \cite{Grace-Martin}:
#'    \itemize{
#'      \item It can help you determine whether or not a linear mixed model is even necessary. If you find that the correlation is zero, that means the observations within clusters are no more similar than observations from different clusters. Go ahead and use a simpler analysis technique.
#'      \item It can be theoretically meaningful to understand how much of the overall variation in the response is explained simply by clustering. For example, in a repeated measures psychological study you can tell to what extent mood is a trait (varies among people, but not within a person on different occasions) or state (varies little on average among people, but varies a lot across occasions).
#'      \item It can also be meaningful to see how the ICC (as well as the between and within cluster variances) changes as variable are added to the model.
#'    }
#'    In short, the ICC can be interpreted as \dQuote{the proportion of the variance
#'    explained by the grouping structure in the population} \cite{(Hox 2002: 15)}.
#'   \cr \cr
#'   The random effect variances indicate the between- and within-group
#'   variances as well as random-slope variance and random-slope-intercept
#'   correlation. The components are denoted as following:
#'   \itemize{
#'     \item Within-group (residual) variance: sigma_2
#'     \item Between-group-variance: tau.00 (variation between individual intercepts and average intercept)
#'     \item Random-slope-variance: tau.11 (variation between individual slopes and average slope)
#'     \item Random-Intercept-Slope-covariance: tau.01
#'     \item Random-Intercept-Slope-correlation: rho.01
#'   }
#'
#' @details The "simple" ICC (with both \code{ppd} and \code{adjusted} set to
#'    \code{FALSE}) is calculated by dividing the between-group-variance (random
#'    intercept variance) by the total variance (i.e. sum of between-group-variance
#'    and within-group (residual) variance). \cr \cr
#'    The calculation of the ICC for generalized linear mixed models with binary outcome is based on
#'    \cite{Wu et al. (2012)}. For other distributions (negative binomial, poisson, ...),
#'    calculation is based on \cite{Nakagawa et al. 2017}, \strong{however}, for
#'    non-Gaussian models it is recommended to compute the adjusted ICC (with
#'    \code{adjusted = TRUE}, see below).
#'    \cr \cr
#'    \strong{ICC for unconditional and conditional models}
#'    \cr \cr
#'    Usually, the ICC is calculated for the null model ("unconditional model").
#'    However, according to \cite{Raudenbush and Bryk (2002)} or
#'    \cite{Rabe-Hesketh and Skrondal (2012)} it is also feasible to compute the ICC
#'    for full models with covariates ("conditional models") and compare how
#'    much a level-2 variable explains the portion of variation in the grouping
#'    structure (random intercept).
#'    \cr \cr
#'    \strong{ICC for random-slope models}
#'    \cr \cr
#'    \strong{Caution:} For models with random slopes and random intercepts,
#'    the ICC would differ at each unit of the predictors. Hence, the ICC for these
#'    kind of models cannot be understood simply as proportion of variance
#'    (see \cite{Goldstein et al. 2010}). For convenience reasons, as the
#'    \code{icc()} function also extracts the different random effects
#'    variances, the ICC for random-slope-intercept-models is reported
#'    nonetheless, but it is usually no meaningful summary of the
#'    proportion of variances.
#'    \cr \cr
#'    To get a meaningful ICC also for models with random slopes, use \code{adjusted = TRUE}.
#'    The adjusted ICC uses the mean random effect variance, which is based
#'    on the random effect variances for each value of the random slope
#'    (see \cite{Johnson et al. 2014}).
#'    \cr \cr
#'    \strong{ICC for models with multiple or nested random effects}
#'    \cr \cr
#'    \strong{Caution:} By default, for three-level-models, depending on the
#'    nested structure of the model, or for models with multiple random effects,
#'    \code{icc()} only reports the proportion of variance explained for each
#'    grouping level. Use \code{adjusted = TRUE} to calculate the adjusted and
#'    conditional ICC, which condition on \emph{all random effects}.
#'    \cr \cr
#'    \strong{Adjusted and conditional ICC}
#'    \cr \cr
#'    If \code{adjusted = TRUE}, an adjusted and conditional ICC are calculated,
#'    which take all sources of uncertainty (of \emph{all random effects})
#'    into account to report an "adjusted" ICC, as well as the conditional ICC.
#'    The latter also takes the fixed effects variances into account (see
#'    \cite{Nakagawa et al. 2017}). If random effects are not nested and not
#'    cross-classified, the adjusted (\code{adjusted = TRUE}) and unadjusted
#'    (\code{adjusted = FALSE}) ICC are identical. \code{adjust = TRUE} returns
#'    a meaningful ICC for models with random slopes. Furthermore, the adjusted
#'    ICC is recommended for models with other distributions than Gaussian.
#'    \cr \cr
#'    \strong{ICC for specific group-levels}
#'    \cr \cr
#'    To calculate the proportion of variance for specific levels related to each
#'    other (e.g., similarity of level-1-units within
#'    level-2-units or level-2-units within level-3-units) must be computed
#'    manually. Use \code{\link{get_re_var}} to get the between-group-variances
#'    and residual variance of the model, and calculate the ICC for the various level
#'    correlations.
#'    \cr \cr
#'    For example, for the ICC between level 1 and 2: \cr
#'    \code{sum(get_re_var(fit)) / (sum(get_re_var(fit)) + get_re_var(fit, "sigma_2"))}
#'    \cr \cr
#'    or for the ICC between level 2 and 3: \cr
#'    \code{get_re_var(fit)[2] / sum(get_re_var(fit))}
#'    \cr \cr
#'    \strong{ICC for Bayesian models}
#'    \cr \cr
#'    If \code{ppd = TRUE}, \code{icc()} calculates a variance decomposition based on
#'    the posterior predictive distribution. In this case, first, the draws from
#'    the posterior predictive distribution \emph{not conditioned} on group-level
#'    terms (\code{posterior_predict(..., re.form = NA)}) are calculated as well
#'    as draws from this distribution \emph{conditioned} on \emph{all random effects}
#'    (by default, unless specified else in \code{re.form}) are taken. Then, second,
#'    the variances for each of these draws are calculated. The "ICC" is then the
#'    ratio between these two variances. This is the recommended way to
#'    analyse random-effect-variances for non-Gaussian models. It is then possible
#'    to compare variances accross models, also by specifying different group-level
#'    terms via the \code{re.form}-argument.
#'    \cr \cr
#'    Sometimes, when the variance of the posterior predictive distribution is
#'    very large, the variance ratio in the output makes no sense, e.g. because
#'    it is negative. In such cases, it might help to use a more robust measure
#'    to calculate the central tendency of the variances. For example, use
#'    \code{typical = "median"}.
#'
#' @seealso \code{\link{re_var}}
#'
#' @examples
#' library(lme4)
#' fit0 <- lmer(Reaction ~ 1 + (1 | Subject), sleepstudy)
#' icc(fit0)
#'
#' # note: ICC for random-slope-intercept model usually not
#' # meaningful, unless you use "adjusted = TRUE" - see 'Note'.
#' fit1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' icc(fit1)
#' icc(fit1, adjusted = TRUE)
#'
#' sleepstudy$mygrp <- sample(1:45, size = 180, replace = TRUE)
#' fit2 <- lmer(Reaction ~ Days + (1 | mygrp) + (1 | Subject), sleepstudy)
#' icc(fit2)
#' icc(fit2, adjusted = TRUE)
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
#'   icc(m)
#'
#'   # show 50% interval
#'   icc(m, prob = .5)
#'
#'   # adjusted ICC, 89% interval
#'   icc(m, adjusted = TRUE)
#'
#'   # variances based on posterior predictive distribution
#'   icc(m, ppd = TRUE)
#' }}
#'
#' @importFrom purrr map2
#' @export
icc <- function(x, ...) {
  UseMethod("icc")
}


#' @importFrom lme4 VarCorr fixef getME
#' @importFrom stats formula
#' @importFrom purrr map map_dbl map_lgl map2
#' @rdname icc
#' @export
icc.merMod <- function(x, adjusted = FALSE, ...) {

  add.args <- lapply(match.call(expand.dots = F)$`...`, function(x) x)

  if (obj_has_name(add.args, "type"))
    type <- add.args[["type"]]
  else
    type <- "icc"

  # compute adjusted and conditional ICC
  if (adjusted) return(r2_mixedmodel(x, type = type, obj.name = deparse(substitute(x))))

  # get family
  fitfam <- model_family(x)


  # random effects variances
  # for details on tau and sigma, see
  # Aguinis H, Gottfredson RK, Culpepper SA2013. Best-Practice Recommendations
  # for Estimating Cross-Level Interaction Effects Using Multilevel Modeling.
  # Journal of Management 39(6): 1490–1528. doi:10.1177/0149206313478188.
  reva <- lme4::VarCorr(x)

  # retrieve only intercepts
  vars <- purrr::map(reva, ~ .x[1])

  # random intercept-variances, i.e.
  # between-subject-variance (tau 00)
  tau.00 <- purrr::map_dbl(vars, ~ .x)

  # random slope-variances (tau 11)
  tau.11 <- unlist(lapply(reva, function(x) diag(x)[-1]))

  # residual variances, i.e. within-cluster-variance
  resid.var <- get_residual_variance(x, var.cor = reva, fitfam, type = "ICC")

  # total variance, sum of random intercept and residual variances
  total_var <- sum(purrr::map_dbl(vars, ~ sum(.x)), resid.var)

  # random intercept icc
  ri.icc <- tau.00 / total_var


  # get random slope random intercept correlations
  # do we have any rnd slopes?

  has_rnd_slope <- purrr::map_lgl(reva, ~ dim(attr(.x, "correlation"))[1] > 1)
  tau.01 <- rho.01 <- NULL


  # get rnd slopes

  if (any(has_rnd_slope)) {

    rnd_slope <- reva[has_rnd_slope]

    # get slope-intercept-correlations
    rho.01 <- purrr::map(rnd_slope, ~ attr(.x, "correlation")[-1, 1])

    # get standard deviations, multiplied
    std_ <- purrr::map(rnd_slope, function(.x) {
      .sd <- attr(.x, "stddev")
      .sd[1] * .sd[2:length(.sd)]
    })

    # bind to matrix
    std.matrix <- purrr::map2(std_, rho.01, ~ cbind(.y, .x))
    tau.01 <- purrr::map(std.matrix, ~ apply(.x, MARGIN = 1, FUN = prod))

    rho.01 <- unlist(rho.01)
    tau.01 <- unlist(tau.01)

    message("Caution! ICC for random-slope-intercept models usually not meaningful. Use `adjusted = TRUE` to use the mean random effect variance to calculate the ICC. See 'Note' in `?icc`.")
  }

  # name values
  names(ri.icc) <- names(reva)


  if (inherits(x, "glmerMod"))
    mt <- "Generalized linear mixed model"
  else
    mt <- "Linear mixed model"

  # add attributes, for print method
  class(ri.icc) <- c("sj_icc_merMod", class(ri.icc))
  attr(ri.icc, "family") <- fitfam$family
  attr(ri.icc, "link") <- fitfam$link.fun
  attr(ri.icc, "formula") <- stats::formula(x)
  attr(ri.icc, "model") <- mt
  attr(ri.icc, "tau.00") <- tau.00
  attr(ri.icc, "tau.01") <- tau.01
  attr(ri.icc, "rho.01") <- rho.01
  attr(ri.icc, "tau.11") <- tau.11
  attr(ri.icc, "sigma_2") <- resid.var
  attr(ri.icc, "rnd.slope.model") <- any(has_rnd_slope)


  # finally, save name of fitted model object. May be needed for
  # the 'se()' function, which accesses the global environment

  attr(ri.icc, ".obj.name") <- deparse(substitute(x))

  # return results
  ri.icc
}


#' @importFrom lme4 VarCorr fixef getME
#' @importFrom glmmTMB VarCorr fixef getME
#' @importFrom stats family formula
#' @importFrom purrr map map_dbl map_lgl
#' @rdname icc
#' @export
icc.glmmTMB <- function(x, adjusted = FALSE, ...) {

  add.args <- lapply(match.call(expand.dots = F)$`...`, function(x) x)

  if (obj_has_name(add.args, "type"))
    type <- add.args[["type"]]
  else
    type <- "icc"

  # compute adjusted and conditional ICC
  if (adjusted) return(r2_mixedmodel(x, type = type, obj.name = deparse(substitute(x))))

  # get family
  fitfam <- model_family(x)


  # random effects variances
  # for details on tau and sigma, see
  # Aguinis H, Gottfredson RK, Culpepper SA2013. Best-Practice Recommendations
  # for Estimating Cross-Level Interaction Effects Using Multilevel Modeling.
  # Journal of Management 39(6): 1490–1528. doi:10.1177/0149206313478188.
  reva <- glmmTMB::VarCorr(x)[[1]]

  # retrieve only intercepts
  vars <- purrr::map(reva, ~ .x[1])

  # random intercept-variances, i.e.
  # between-subject-variance (tau 00)
  tau.00 <- purrr::map_dbl(vars, ~ .x)

  # random slope-variances (tau 11)
  tau.11 <- unlist(lapply(reva, function(x) diag(x)[-1]))

  # residual variances, i.e. within-cluster-variance
  resid.var <- get_residual_variance(x, var.cor = reva, fitfam, type = "ICC")

  # total variance, sum of random intercept and residual variances
  total_var <- sum(purrr::map_dbl(vars, ~ sum(.x)), resid.var)

  # random intercept icc
  ri.icc <- tau.00 / total_var


  # get random slope random intercept correlations
  # do we have any rnd slopes?

  has_rnd_slope <- purrr::map_lgl(reva, ~ dim(attr(.x, "correlation"))[1] > 1)
  tau.01 <- rho.01 <- NULL


  # get rnd slopes

  if (any(has_rnd_slope)) {

    rnd_slope <- reva[has_rnd_slope]

    # get slope-intercept-correlations
    rho.01 <- purrr::map(rnd_slope, ~ attr(.x, "correlation")[-1, 1])

    # get standard deviations, multiplied
    std_ <- purrr::map(rnd_slope, function(.x) {
      .sd <- attr(.x, "stddev")
      .sd[1] * .sd[2:length(.sd)]
    })

    # bind to matrix
    std.matrix <- purrr::map2(std_, rho.01, ~ cbind(.y, .x))
    tau.01 <- purrr::map(std.matrix, ~ apply(.x, MARGIN = 1, FUN = prod))

    rho.01 <- unlist(rho.01)
    tau.01 <- unlist(tau.01)

    message("Caution! ICC for random-slope-intercept models usually not meaningful. Use `adjusted = TRUE` to use the mean random effect variance to calculate the ICC. See 'Note' in `?icc`.")
  }

  # name values
  names(ri.icc) <- names(reva)


  mt <- "Generalized linear mixed model"

  # add attributes, for print method
  class(ri.icc) <- c("sj_icc_merMod", class(ri.icc))
  attr(ri.icc, "family") <- fitfam$family
  attr(ri.icc, "link") <- fitfam$link.fun
  attr(ri.icc, "formula") <- stats::formula(x)
  attr(ri.icc, "model") <- mt
  attr(ri.icc, "tau.00") <- tau.00
  attr(ri.icc, "tau.01") <- tau.01
  attr(ri.icc, "rho.01") <- rho.01
  attr(ri.icc, "tau.11") <- tau.11
  attr(ri.icc, "sigma_2") <- resid.var
  attr(ri.icc, "rnd.slope.model") <- any(has_rnd_slope)


  # finally, save name of fitted model object. May be needed for
  # the 'se()' function, which accesses the global environment

  attr(ri.icc, ".obj.name") <- deparse(substitute(x))

  # return results
  ri.icc
}


#' @importFrom stats formula
#' @importFrom purrr map map_dbl map_lgl
#' @importFrom sjmisc row_sums is_empty
#' @rdname icc
#' @export
icc.stanreg <- function(x, re.form = NULL, typical = "mean", prob = .89, ppd = FALSE, adjusted = FALSE, ...) {

  if (!requireNamespace("rstanarm", quietly = TRUE))
    stop("Please install and load package `rstanarm` first.", call. = F)

  # get family
  fitfam <- model_family(x)
  xdat <- as.data.frame(x)

  if (missing(ppd) && missing(adjusted) && !fitfam$is_linear) {
    #message("Variance decomposition is based on the posterior predictive distribution. Set `ppd = FALSE` to calculate \"classical\" ICC, and `adjusted = TRUE` for adjusted ICC.")
    message("Variance decomposition for non-Gaussian models should be based on the posterior predictive distribution. To do this, set `ppd = TRUE`.")
    ## TODO set ppd to FALSE by default later
    # ppd <- TRUE
  }

  if (ppd) {

    ## TODO automatically calculate for multiple levels / nested models

    PPD <- rstanarm::posterior_predict(x, re.form = re.form)
    total_var <- apply(PPD, MARGIN = 1, FUN = stats::var)

    PPD_0 <- rstanarm::posterior_predict(x, re.form = NA)
    tau.00 <- apply(PPD_0, MARGIN = 1, FUN = stats::var)

    ri.icc <- tau.00 / total_var
    resid.var <- total_var - tau.00

    icc_ <- c(
      1 - typical_value(ri.icc, fun = typical),
      typical_value(tau.00, fun = typical),
      typical_value(resid.var, fun = typical),
      typical_value(total_var, fun = typical)
    )

    attr(icc_, "hdi.icc") <- rev(1 - hdi(ri.icc, prob = prob))
    attr(icc_, "hdi.tau.00") <- hdi(tau.00, prob = prob)
    attr(icc_, "hdi.resid") <- hdi(resid.var, prob = prob)
    attr(icc_, "hdi.total") <- hdi(total_var, prob = prob)
    attr(icc_, "re.form") <- re.form
    attr(icc_, "ranef") <- x$ranef$group[1]

    has_rnd_slope <- FALSE
    names(icc_) <- c("icc", "tau.00", "resid.var", "total.var")
    class(icc_) <- c("icc_ppd", class(icc_))

  } else if (adjusted) {
    # compute adjusted and conditional ICC
    return(r2_mixedmodel(x, type = "ICC", obj.name = deparse(substitute(x))))
  } else {

    # random intercept-variances, i.e.
    # between-subject-variance (tau 00)
    tau00s <- grep("^Sigma\\[(.*):\\(Intercept\\),\\(Intercept\\)", colnames(xdat))
    tau.00 <- xdat[, tau00s, drop = FALSE]

    names(tau.00) <- gsub(
      "^Sigma\\[(.*):\\(Intercept\\),\\(Intercept\\)\\]",
      "\\1",
      colnames(xdat))[tau00s]


    # random slope-variances (tau 11)
    tau11s <- grep("^Sigma\\[(.*):[^\\(\\)](.*),[^\\(\\)](.*)\\]", colnames(xdat))

    if (!sjmisc::is_empty(tau11s)) {
      tau.11 <- xdat[, tau11s, drop = FALSE]
      names(tau.11) <- gsub(
        "^Sigma\\[(.*):[^\\(\\)](.*),[^\\(\\)](.*)\\]",
        "\\1",
        colnames(xdat))[tau11s]
    } else {
      tau.11 <- NULL
    }


    # get residual standard deviation sigma
    sig <- xdat[["sigma"]]

    # set default, if no residual variance is available
    if (is.null(sig)) {
      if (fitfam$is_bin)
        sig <- sqrt((pi^2) / 3)
      else
        sig <- 1
    }

    # residual variance
    resid.var <- sig^2

    # total variance, sum of random intercept and residual variances
    total_var <- sjmisc::row_sums(
      cbind(tau.00, data.frame(resid.var)),
      n = 1,
      var = "total_var",
      append = FALSE
    )

    # random intercept icc
    ri.icc <- tau.00 / total_var$total_var


    # get random slope random intercept correlations
    # do we have any rnd slopes?

    has_rnd_slope <- !sjmisc::is_empty(tau11s)
    tau.01 <- rho.01 <- NULL


    # get rnd slopes

    if (any(has_rnd_slope)) {

      # get slope-intercept-covariance
      tau01s <- grep("^Sigma\\[(.*):[^\\(\\)](.*),\\(Intercept\\)\\]", colnames(xdat))
      tau.01 <- xdat[, tau01s, drop = FALSE]

      tau.00.sums <- sjmisc::row_sums(tau.00, var = "t0sums", n = 1)$t0sums
      tau.11.sums <- sjmisc::row_sums(tau.11, var = "t1sums", n = 1)$t1sums

      # get slope-intercept-correlations
      rho.01 <- tau.01 / sqrt(tau.00.sums * tau.11.sums)

      message("Caution! ICC for random-slope-intercept models usually not meaningful. See 'Note' in `?icc`.")

    }

    if (inherits(x, "glmerMod"))
      mt <- "Generalized linear mixed model"
    else
      mt <- "Linear mixed model"


    icc_ <- purrr::map_dbl(ri.icc, ~ typical_value(.x, fun = typical))
    attr(icc_, "hdi.icc") <- purrr::map(ri.icc, ~ hdi(.x, prob = prob))

    attr(icc_, "hdi.tau.00") <- purrr::map(tau.00, ~ hdi(.x, prob = prob))
    tau.00 <- purrr::map_dbl(tau.00, ~ typical_value(.x, fun = typical))

    if (length(resid.var) > 10)
      attr(icc_, "hdi.sigma_2") <- hdi(resid.var, prob = prob)
    resid.var <- typical_value(resid.var, fun = typical)

    if (!is.null(tau.11)) {
      attr(icc_, "hdi.tau.11") <- purrr::map(tau.11, ~ hdi(.x, prob = prob))
      tau.11 <- purrr::map_dbl(tau.11, ~ typical_value(.x, fun = typical))
    }

    if (!is.null(rho.01)) {
      attr(icc_, "hdi.rho.01") <- purrr::map(rho.01, ~ hdi(.x, prob = prob))
      rho.01 <- purrr::map_dbl(rho.01, ~ typical_value(.x, fun = typical))
    }

    if (!is.null(tau.01)) {
      attr(icc_, "hdi.tau.01") <- purrr::map(tau.01, ~ hdi(.x, prob = prob))
      tau.01 <- purrr::map_dbl(tau.01, ~ typical_value(.x, fun = typical))
    }

    attr(icc_, "tau.00") <- tau.00
    attr(icc_, "tau.01") <- tau.01
    attr(icc_, "rho.01") <- rho.01
    attr(icc_, "tau.11") <- tau.11
    attr(icc_, "sigma_2") <- resid.var
    attr(icc_, "rnd.slope.model") <- any(has_rnd_slope)
    attr(icc_, "model") <- mt

    class(icc_) <- c("sj_icc_stanreg", class(icc_))
  }

  # add attributes, for print method
  attr(icc_, "family") <- fitfam$family
  attr(icc_, "link") <- fitfam$link.fun
  attr(icc_, "formula") <- stats::formula(x)
  attr(icc_, "prob") <- prob


  # return results
  icc_
}


#' @importFrom purrr map_df map_if map_lgl map_dbl
#' @importFrom dplyr bind_cols
#' @importFrom sjmisc all_na
#' @rdname icc
#' @export
icc.brmsfit <- function(x, re.form = NULL, typical = "mean", prob = .89, ppd = FALSE, adjusted = FALSE, ...) {

  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.", call. = F)

  # get family
  fitfam <- model_family(x)

  if (missing(ppd) && missing(adjusted) && !fitfam$is_linear) {
    #message("Variance decomposition is based on the posterior predictive distribution. Set `ppd = FALSE` to calculate \"classical\" ICC, and `adjusted = TRUE` for adjusted ICC.")
    message("Variance decomposition for non-Gaussian models should be based on the posterior predictive distribution. To do this, set `ppd = TRUE`.")
    ## TODO set ppd to FALSE by default later
    # ppd <- TRUE
  }


  if (ppd) {

    ## TODO automatically calculate for multiple levels / nested models

    PPD <- brms::posterior_predict(x, re.form = re.form, summary = FALSE)
    total_var <- apply(PPD, MARGIN = 1, FUN = stats::var)

    PPD_0 <- brms::posterior_predict(x, re.form = NA, summary = FALSE)
    tau.00 <- apply(PPD_0, MARGIN = 1, FUN = stats::var)

    ri.icc <- tau.00 / total_var
    resid.var <- total_var - tau.00

    icc_ <- c(
      1 - typical_value(ri.icc, fun = typical),
      typical_value(tau.00, fun = typical),
      typical_value(resid.var, fun = typical),
      typical_value(total_var, fun = typical)
    )

    attr(icc_, "hdi.icc") <- rev(1 - hdi(ri.icc, prob = prob))
    attr(icc_, "hdi.tau.00") <- hdi(tau.00, prob = prob)
    attr(icc_, "hdi.resid") <- hdi(resid.var, prob = prob)
    attr(icc_, "hdi.total") <- hdi(total_var, prob = prob)
    attr(icc_, "prob") <- prob
    attr(icc_, "re.form") <- re.form
    attr(icc_, "ranef") <- x$ranef$group[1]

    has_rnd_slope <- FALSE
    names(icc_) <- c("icc", "tau.00", "resid.var", "total.var")
    class(icc_) <- c("icc_ppd", class(icc_))

  } else if (adjusted) {
    # compute adjusted and conditional ICC
    return(r2_mixedmodel(x, type = "ICC", obj.name = deparse(substitute(x))))
  } else {

    # get random effect variances for each sample of posterior
    reva <- brms::VarCorr(x, summary = FALSE)

    # remove "residual__" element from list
    # and save in separate object
    reva.resid <- reva[names(reva) == "residual__"]
    reva <- reva[!(names(reva) == "residual__")]


    # retrieve only intercepts
    vars <- purrr::map(reva, ~ .x$sd[, 1]^2)

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
        sig <- sqrt((pi^2) / 3)
      else
        sig <- 1
    }


    # residual variances, i.e.
    # within-cluster-variance (sigma^2)

    resid.var <- sig^2


    # total variance, sum of random intercept and residual variances
    total_var <- apply(as.data.frame(vars), MARGIN = 1, FUN = sum) + resid.var

    # make sure residual variance has same length as other components
    # if not, just repeat the current value to match number of samples
    if (length(resid.var) == 1) resid.var <- rep(resid.var, length(total_var))

    # random intercept icc
    ri.icc <- purrr::map(tau.00, ~ .x / total_var)

    # random slope-variances (tau 11)
    tau.11 <- purrr::map_if(tau.11, is.null, ~ rep(NA, length(resid.var)))

    names(ri.icc) <- sprintf("icc_%s", names(ri.icc))
    names(tau.00) <- sprintf("tau.00_%s", names(tau.00))
    names(tau.11) <- sprintf("tau.11_%s", names(tau.11))

    icc_ <- purrr::map_dbl(ri.icc, ~ typical_value(.x, fun = typical))

    attr(icc_, "tau.00") <- purrr::map_dbl(tau.00, ~ typical_value(.x, fun = typical))
    attr(icc_, "hdi.icc") <- purrr::map(ri.icc, ~ hdi(.x, prob = prob))
    attr(icc_, "hdi.tau.00") <- purrr::map(tau.00, ~ hdi(.x, prob = prob))

    attr(icc_, "sigma_2") <- typical_value(resid.var, fun = typical)
    attr(icc_, "hdi.sigma_2") <- hdi(resid.var, prob = prob)

    attr(icc_, "prob") <- prob

    check_tau <- purrr::map_lgl(tau.11, ~ sjmisc::all_na(.x))
    if (any(!check_tau)) {
      tau.11 <- tau.11[!check_tau]
      attr(icc_, "tau.11") <- purrr::map_dbl(tau.11, ~ typical_value(.x, fun = typical))
      attr(icc_, "hdi.tau.11") <- purrr::map(tau.11, ~ hdi(.x, prob = prob))
    }

    has_rnd_slope <- any(isTRUE(purrr::map_lgl(brms::ranef(x), ~ dim(.x)[3] > 1)))

    if (has_rnd_slope)
      message("Caution! ICC for random-slope-intercept models usually not meaningful. See 'Note' in `?icc`.")

    class(icc_) <- c("sj_icc_brms", class(icc_))
  }


  attr(icc_, "family") <- fitfam$family
  attr(icc_, "link") <- fitfam$link.fun
  attr(icc_, "formula") <- stats::formula(x)
  attr(icc_, "model") <- "Bayesian mixed model"
  attr(ri.icc, "rnd.slope.model") <- any(has_rnd_slope)

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
#'   \code{stanreg} or \code{brmsfit}). \code{get_re_var()} also accepts
#'   an object returned by the \code{\link{icc}} function.
#' @param comp Name of the variance component to be returned. See 'Details'.
#' @param adjusted Logical, if \code{TRUE}, returns the variance of the fixed
#'   and random effects as well as of the additive dispersion and
#'   distribution-specific variance, which are used to calculate the
#'   adjusted and conditional \code{\link{r2}} and \code{\link{icc}}.
#'
#' @return \code{get_re_var()} returns the value of the requested variance component,
#'           \code{re_var()} returns all random effects variances.
#'
#' @references \itemize{
#'    \item Aguinis H, Gottfredson RK, Culpepper SA. 2013. Best-Practice Recommendations for Estimating Cross-Level Interaction Effects Using Multilevel Modeling. Journal of Management 39(6): 1490–1528 (\doi{10.1177/0149206313478188})
#'    \item Johnson PC, O'Hara RB. 2014. Extension of Nakagawa & Schielzeth's R2GLMM to random slopes models. Methods Ecol Evol, 5: 944-946. (\doi{10.1111/2041-210X.12225})
#'    \item Nakagawa S, Johnson P, Schielzeth H (2017) The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisted and expanded. J. R. Soc. Interface 14. \doi{10.1098/rsif.2017.0213}
#'  }
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
#'         \cr \cr
#'         If \code{adjusted = TRUE}, the variance of the fixed and random
#'         effects as well as of the additive dispersion and
#'         distribution-specific variance are returned (see \cite{Johnson et al. 2014}
#'         and \cite{Nakagawa et al. 2017}):
#'         \describe{
#'          \item{\code{"fixed"}}{variance attributable to the fixed effects}
#'          \item{\code{"random"}}{(mean) variance of random effects}
#'          \item{\code{"dispersion"}}{variance due to additive dispersion}
#'          \item{\code{"distribution"}}{distribution-specific variance}
#'          \item{\code{"residual"}}{sum of dispersion and distribution}
#'         }
#'
#' @seealso \code{\link{icc}}
#'
#' @examples
#' library(lme4)
#' fit1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#'
#' # all random effect variance components
#' re_var(fit1)
#' re_var(fit1, adjusted = TRUE)
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
re_var <- function(x, adjusted = FALSE) {

  if (adjusted) {

    rv <- r2_mixedmodel(x)

    rv_ <- list(
      var.fixef = attr(rv, "var.fixef", exact = TRUE),
      var.ranef = attr(rv, "var.ranef", exact = TRUE),
      var.disp = attr(rv, "var.disp", exact = TRUE),
      var.dist = attr(rv, "var.dist", exact = TRUE),
      var.resid = attr(rv, "var.resid", exact = TRUE),
      formula = attr(rv, "formula", exact = TRUE),
      family = attr(rv, "family", exact = TRUE),
      link = attr(rv, "link", exact = TRUE)
    )

    class(rv_) <- c("sj_revar_adjust", class(rv_))

  } else {
    # iterate all attributes and return them as vector
    rv <- c("sigma_2", "tau.00", "tau.11", "tau.01", "rho.01")

    # compute icc
    icc_ <- suppressMessages(icc(x, ppd = FALSE))

    rv_ <- purrr::map(rv, ~ attr(icc_, .x, exact = TRUE))
    rn <- purrr::map2(1:length(rv_), rv, ~ sjmisc::trim(paste(names(rv_[[.x]]), .y, sep = "_")))
    rv_ <- purrr::flatten_dbl(rv_)

    names(rv_) <- purrr::flatten_chr(rn)[1:length(rv_)]

    class(rv_) <- c("sj_revar", class(rv_))
  }

  rv_
}


#' @rdname re_var
#' @export
get_re_var <- function(x, comp = c("tau.00", "tau.01", "tau.11", "rho.01", "sigma_2")) {
  # check if we have a valid object
  if (!inherits(x, c("sj_icc_merMod", "sj_icc_stanreg", "sj_icc_brms")) && !is_merMod(x) && !inherits(x, c("glmmTMB", "brmsfit"))) {
    stop("`x` must either be an object returned by the `icc()` function, or a merMod-, glmmTMB- or brmsfit-object.", call. = F)
  }

  # check arguments
  comp <- match.arg(comp)

  # do we have a merMod object? If yes, get ICC and var components
  if (is_merMod(x) || inherits(x, c("glmmTMB", "brmsfit"))) x <- suppressMessages(icc(x, ppd = FALSE))

  # return results
  attr(x, comp, exact = TRUE)
}
