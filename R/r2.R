#' @title Goodness-of-fit measures for regression models
#' @name cod
#'
#' @description Compute Goodness-of-fit measures for various regression models,
#'   including mixed and Bayesian regression models.
#'
#' @param x Fitted model of class \code{lm}, \code{glm}, \code{merMod},
#'    \code{glmmTMB}, \code{lme}, \code{plm}, \code{stanreg} or \code{brmsfit}.
#'    For method \code{cod()}, only a \code{glm} with binrary response.
#' @param n Optional, an \code{lme} object, representing the fitted null-model
#'    (unconditional model) to \code{x}. If \code{n} is given, the pseudo-r-squared
#'    for random intercept and random slope variances are computed
#'    (\cite{Kwok et al. 2008}) as well as the Omega squared value
#'    (\cite{Xu 2003}). See 'Examples' and 'Details'.
#' @param loo Logical, if \code{TRUE} and \code{x} is a \code{stanreg} or
#'    \code{brmsfit} object, a LOO-adjusted r-squared is calculated. Else,
#'    a rather "unadjusted" r-squared will be returned by calling
#'    \code{rstantools::bayes_R2()}.
#' @param ... Currently not used.
#'
#' @return For \code{r2()}, depending on the model, returns:
#'         \itemize{
#'           \item For linear models, the r-squared and adjusted r-squared values.
#'           \item For mixed models, the marginal and conditional r-squared values.
#'           \item For \code{glm} objects, Cox & Snell's and Nagelkerke's pseudo r-squared values.
#'           \item For \code{brmsfit} or \code{stanreg} objects, the Bayesian version of r-squared is computed, calling \code{rstantools::bayes_R2()}.
#'           \item If \code{loo = TRUE}, for \code{brmsfit} or \code{stanreg} objects a LOO-adjusted version of r-squared is returned.
#'           \item Models that are not currently supported return \code{NULL}.
#'         }
#'         For \code{cod()}, returns the \code{D} Coefficient of Discrimination,
#'         also known as Tjur's R-squared value.
#'
#' @references \itemize{
#'               \item \href{http://glmm.wikidot.com/faq}{DRAFT r-sig-mixed-models FAQ}
#'               \item Bolker B et al. (2017): \href{http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html}{GLMM FAQ}
#'               \item Byrnes, J. 2008. Re: Coefficient of determination (R^2) when using lme() (\url{https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q2/000713.html})
#'               \item Kwok OM, Underhill AT, Berry JW, Luo W, Elliott TR, Yoon M. 2008. Analyzing Longitudinal Data with Multilevel Models: An Example with Individuals Living with Lower Extremity Intra-Articular Fractures. Rehabilitation Psychology 53(3): 370-86. \doi{10.1037/a0012765}
#'               \item Nakagawa S, Schielzeth H. 2013. A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution, 4(2):133-142. \doi{10.1111/j.2041-210x.2012.00261.x}
#'               \item Nakagawa S, Johnson P, Schielzeth H (2017) The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisted and expanded. J. R. Soc. Interface 14. \doi{10.1098/rsif.2017.0213}
#'               \item Rabe-Hesketh S, Skrondal A. 2012. Multilevel and longitudinal modeling using Stata. 3rd ed. College Station, Tex: Stata Press Publication
#'               \item Raudenbush SW, Bryk AS. 2002. Hierarchical linear models: applications and data analysis methods. 2nd ed. Thousand Oaks: Sage Publications
#'               \item Snijders TAB, Bosker RJ. 2012. Multilevel analysis: an introduction to basic and advanced multilevel modeling. 2nd ed. Los Angeles: Sage
#'               \item Xu, R. 2003. Measuring explained variation in linear mixed effects models. Statist. Med. 22:3527-3541. \doi{10.1002/sim.1572}
#'               \item Tjur T. 2009. Coefficients of determination in logistic regression models - a new proposal: The coefficient of discrimination. The American Statistician, 63(4): 366-372
#'             }
#'
#' @details For linear models, the r-squared and adjusted r-squared value is returned,
#'          as provided by the \code{summary}-function.
#'          \cr \cr
#'          For mixed models (from \pkg{lme4} or \pkg{glmmTMB}) marginal and
#'          conditional r-squared values are calculated, based on
#'          \cite{Nakagawa et al. 2017}. The distributional variance
#'          (or observation-level variance) is based on lognormal approximation,
#'          \code{log(1+var(x)/mu^2)}.
#'          \cr \cr
#'          For \code{lme}-models, an r-squared approximation by computing the
#'          correlation between the fitted and observed values, as suggested by
#'          \cite{Byrnes (2008)}, is returned as well as a simplified version of
#'          the Omega-squared value (1 - (residual variance / response variance),
#'          \cite{Xu (2003)}, \cite{Nakagawa, Schielzeth 2013}), unless \code{n}
#'          is specified.
#'          \cr \cr
#'          If \code{n} is given, for \code{lme}-models pseudo r-squared measures based
#'          on the variances of random intercept (tau 00, between-group-variance)
#'          and random slope (tau 11, random-slope-variance), as well as the
#'          r-squared statistics as proposed by \cite{Snijders and Bosker 2012} and
#'          the Omega-squared value (1 - (residual variance full model / residual
#'          variance null model)) as suggested by \cite{Xu (2003)} are returned.
#'          \cr \cr
#'          For generalized linear models, Cox & Snell's and Nagelkerke's
#'          pseudo r-squared values are returned.
#'          \cr \cr
#'          The ("unadjusted") r-squared value and its standard error for
#'          \code{brmsfit} or \code{stanreg} objects are robust measures, i.e.
#'          the median is used to compute r-squared, and the median absolute
#'          deviation as the measure of variability. If \code{loo = TRUE},
#'          a LOO-adjusted r-squared is calculated, which comes conceptionally
#'          closer to an adjusted r-squared measure.
#'
#' @note \describe{
#'         \item{\strong{cod()}}{
#'          This method calculates the Coefficient of Discrimination \code{D}
#'          for generalized linear (mixed) models for binary data. It is
#'          an alternative to other Pseudo-R-squared values like Nakelkerke's
#'          R2 or Cox-Snell R2. The Coefficient of Discrimination \code{D}
#'          can be read like any other (Pseudo-)R-squared value.
#'         }
#'         \item{\strong{r2()}}{
#'          For mixed models, the marginal r-squared considers only the variance
#'          of the fixed effects, while the conditional r-squared takes both
#'          the fixed and random effects into account.
#'          \cr \cr
#'          For \code{lme}-objects, if \code{n} is given, the Pseudo-R2 statistic
#'          is the proportion of explained variance in the random effect after
#'          adding co-variates or predictors to the model, or in short: the
#'          proportion of the explained variance in the random effect of the
#'          full (conditional) model \code{x} compared to the null (unconditional)
#'          model \code{n}.
#'          \cr \cr
#'          The Omega-squared statistics, if \code{n} is given, is 1 - the proportion
#'          of the residual variance of the full model compared to the null model's
#'          residual variance, or in short: the the proportion of the residual
#'          variation explained by the covariates.
#'          \cr \cr
#'          Alternative ways to assess the "goodness-of-fit" is to compare the ICC
#'          of the null model with the ICC of the full model (see \code{\link{icc}}).
#'         }
#'       }
#'
#'
#' @examples
#' data(efc)
#'
#' # Tjur's R-squared value
#' efc$services <- ifelse(efc$tot_sc_e > 0, 1, 0)
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep,
#'            data = efc, family = binomial(link = "logit"))
#' cod(fit)
#'
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' r2(fit)
#'
#' fit <- lm(barthtot ~ c160age + c12hour, data = efc)
#' r2(fit)
#'
#' # Pseudo-R-squared values
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep,
#'            data = efc, family = binomial(link = "logit"))
#' r2(fit)
#'
#' @importFrom stats predict predict.glm residuals
#' @export
cod <- function(x) {
  # check for valid object class
  if (!inherits(x, c("glmerMod", "glm"))) {
    stop("`x` must be an object of class `glm` or `glmerMod`.", call. = F)
  }

  # mixed models (lme4)
  if (inherits(x, "glmerMod")) {
    # check for package availability
    y <- lme4::getME(x, "y")
    pred <- stats::predict(x, type = "response", re.form = NULL)
  } else {
    y <- x$y
    pred <- stats::predict.glm(x, type = "response")
  }

  # delete pred for cases with missing residuals
  if (anyNA(stats::residuals(x))) pred <- pred[!is.na(stats::residuals(x))]

  categories <- unique(y)
  m1 <- mean(pred[which(y == categories[1])], na.rm = T)
  m2 <- mean(pred[which(y == categories[2])], na.rm = T)

  cod = abs(m2 - m1)
  names(cod) <- "Tjur's D"

  structure(class = "sj_r2", list(cod = cod))
}


#' @rdname cod
#' @importFrom stats model.response fitted var residuals median mad
#' @importFrom sjmisc is_empty
#' @export
r2 <- function(x, ...) {
  UseMethod("r2")
}


#' @export
r2.glmmTMB <- function(x, ...) {
  if (!requireNamespace("glmmTMB", quietly = TRUE))
    stop("Package `glmmTMB` needed for this function to work. Please install it.", call. = FALSE)

  r2_mixedmodel(x, type = "r2", obj.name = deparse(substitute(x)))
}


#' @rdname cod
#' @export
r2.lme <- function(x, n = NULL, ...) {
  r2linmix(x, n)
}


#' @rdname cod
#' @export
r2.stanreg <- function(x, loo = FALSE, ...) {
  if (!requireNamespace("rstanarm", quietly = TRUE))
    stop("Package `rstanarm` needed for this function to work. Please install it.", call. = FALSE)

  if (inherits(x, "stanmvreg"))
    return(r2.stanmvreg(x = x, loo = loo, ...))

  if (isTRUE(loo)) {
    rsq <- looR2(x)
    names(rsq) <- "LOO-adjusted R2"
    structure(class = "sj_r2", list(r2 = rsq))
  } else {
    brs <- rstanarm::bayes_R2(x)
    rsq <- stats::median(brs)
    rsq.se <- stats::mad(brs)

    names(rsq) <- "Bayes R2"
    names(rsq.se) <- "Standard Error"

    structure(class = "sj_r2", list(r2 = rsq, se = rsq.se))
  }
}


#' @export
r2.stanmvreg <- function(x, loo = FALSE, ...) {
  NULL
}


#' @rdname cod
#' @export
r2.brmsfit <- function(x, loo = FALSE, ...) {
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Package `brms` needed for this function to work. Please install it.", call. = FALSE)

  if (isTRUE(loo)) {
    rsq <- looR2(x)
    names(rsq) <- "LOO-adjusted R2"
    structure(class = "sj_r2", list(r2 = rsq))
  } else {
    brs <- brms::bayes_R2(x, summary = TRUE, robust = TRUE)
    rsq <- brs[1]
    rsq.se <- brs[2]

    names(rsq) <- "Bayes R2"
    names(rsq.se) <- "Standard Error"

    structure(class = "sj_r2", list(r2 = rsq, se = rsq.se))
  }
}


#' @importFrom stats nobs
#' @export
r2.glm <- function(x, ...) {
  n <- stats::nobs(x)
  CoxSnell <- (1 - exp((x$dev - x$null) / n))
  Nagelkerke <- CoxSnell / (1 - exp(-x$null / n))

  names(CoxSnell) <- "Cox & Snell's R-squared"
  names(Nagelkerke) <- "Nagelkerke's R-squared"

  structure(class = "sj_r2", list(CoxSnell = CoxSnell, Nagelkerke = Nagelkerke))
}


#' @export
r2.mlogit <- function(x, ...) {
  McFadden <- as.vector(summary(x)$mfR2)
  names(McFadden) <- "McFadden R-squared"
  structure(class = "sj_r2", list(McFadden = McFadden))
}


#' @export
r2.merMod <- function(x, ...) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package `lme4` needed for this function to work. Please install it.", call. = FALSE)

  r2_mixedmodel(x, type = "r2", obj.name = deparse(substitute(x)))
}


#' @export
r2.lm <- function(x, ...) {
  # do we have a simple linear model?
  rsq <- summary(x)$r.squared
  adjr2 <- summary(x)$adj.r.squared

  # name vectors
  names(rsq) <- "R-squared"
  names(adjr2) <- "adjusted R-squared"

  # return results
  structure(class = "sj_r2", list(r2 = rsq, adjr2 = adjr2))
}


#' @export
r2.default <- function(x, ...) {
  tryCatch(
    {
      # do we have a simple linear model?
      rsq <- summary(x)$r.squared
      adjr2 <- summary(x)$adj.r.squared

      # name vectors
      names(rsq) <- "R-squared"
      names(adjr2) <- "adjusted R-squared"

      # return results
      structure(class = "sj_r2", list(r2 = rsq, adjr2 = adjr2))
    },
    error = function(x) { NULL }
  )
}


#' @importFrom stats logLik update
#' @export
r2.polr <- function(x, ...) {
  L.base <- stats::logLik(stats::update(x, ~ 1))
  r2glm(x, L.base)
}


#' @importFrom stats logLik update
#' @export
r2.clm2 <- function(x, ...) {
  L.base <- stats::logLik(stats::update(x, location = ~ 1, scale = ~ 1))
  r2glm(x, L.base)
}


#' @importFrom stats logLik update
#' @export
r2.clm <- function(x, ...) {
  L.base <- stats::logLik(stats::update(x, ~ 1))
  r2glm(x, L.base)
}


#' @importFrom stats logLik update
#' @export
r2.vglm <- function(x, ...) {
  if (!(is.null(x@call$summ) && !identical(x@call$summ, 0)))
    stop("Can't get log-likelihood when `summ` is not zero.", call. = FALSE)

  L.base <- stats::logLik(stats::update(x, ~ 1))
  r2glm(x, L.base)
}


#' @importFrom stats logLik update
#' @export
r2.multinom <- function(x, ...) {
  L.base <- stats::logLik(stats::update(x, ~ 1, trace = FALSE))
  r2glm(x, L.base)
}


#' @export
r2.plm <- function(x, ...) {
  # else do we have a mixed model?
  rsq <- summary(x)$r.squared[1]
  adjr2 <- summary(x)$r.squared[2]

  # name vectors
  names(rsq) <- "R-squared"
  names(adjr2) <- "adjusted R-squared"

  # return results
  structure(class = "sj_r2", list(r2 = rsq, adjr2 = adjr2))
}


#' @importFrom stats var
looR2 <- function(fit) {

  if (!requireNamespace("rstantools", quietly = TRUE))
    stop("Package `rstantools` required. Please install.", call. = FALSE)

  if (!requireNamespace("loo", quietly = TRUE))
    stop("Package `loo` required. Please install.", call. = FALSE)

  y <- resp_val(fit)
  ypred <- rstantools::posterior_linpred(fit)


  # for some weird models, not all response values can be
  # predicted, resulting in different lengths between y and ypred

  if (length(y) > ncol(ypred)) {
    tryCatch(
      {
        y <- y[as.numeric(attr(ypred, "dimnames")[[2]])]
      },
      error = function(x) { NULL }
    )
  }

  ll <- rstantools::log_lik(fit)

  r_eff <- loo::relative_eff(
    exp(ll),
    chain_id = rep(1:n_of_chains(fit), each = n_of_samples(fit) / n_of_chains(fit))
  )

  psis_object <- loo::psis(log_ratios = -ll, r_eff = r_eff)
  ypredloo <- loo::E_loo(ypred, psis_object, log_ratios = -ll)$value
  eloo <- ypredloo - y

  1 - stats::var(eloo) / stats::var(y)
}


r2linmix <- function(x, n) {
  # do we have null model?
  if (!is.null(n)) {
    # compute tau for both models
    tau_full <- suppressMessages(icc(x))
    tau_null <- suppressMessages(icc(n))

    # get taus. tau.00 is the random intercept variance, i.e. for growth models,
    # the difference in the outcome's mean at first time point
    rsq0 <- (attr(tau_null, "tau.00") - attr(tau_full, "tau.00")) / attr(tau_null, "tau.00")

    # tau.11 is the variance of the random slopes, i.e. how model predictors
    # affect the trajectory of subjects over time (for growth models)
    rsq1 <- (attr(tau_null, "tau.11") - attr(tau_full, "tau.11")) / attr(tau_null, "tau.11")

    # get r2
    rsq <- ((attr(tau_null, "tau.00") + attr(tau_null, "sigma_2")) -
              (attr(tau_full, "tau.00") + attr(tau_full, "sigma_2"))) /
      (attr(tau_null, "tau.00") + attr(tau_null, "sigma_2"))

    # get omega-squared
    osq <- 1 - ((attr(tau_full, "sigma_2") / attr(tau_null, "sigma_2")))

    # if model has no random slope, we need to set this value to NA
    if (is.null(rsq1) || sjmisc::is_empty(rsq1)) rsq1 <- NA

    # name vectors
    names(rsq0) <- "R-squared (tau-00)"
    names(rsq1) <- "R-squared (tau-11)"
    names(rsq) <- "R-squared"
    names(osq) <- "Omega-squared"

    # return results
    structure(class = "sj_r2", list(
      r2_tau00 = rsq0,
      r2_tau11 = rsq1,
      r2 = rsq,
      o2 = osq
    ))
  } else {
    # compute "correlation"
    lmfit <-  lm(resp_val(x) ~ stats::fitted(x))
    # get r-squared
    rsq <- summary(lmfit)$r.squared
    # get omega squared
    osq <- 1 - stats::var(stats::residuals(x)) / stats::var(resp_val(x))

    # name vectors
    names(rsq) <- "R-squared"
    names(osq) <- "Omega-squared"

    # return results
    structure(class = "sj_r2", list(r2 = rsq, o2 = osq))
  }
}


#' @importFrom insight find_formula model_info
r2_mixedmodel <- function(x, type = NULL, obj.name = NULL) {

  if (is.null(type) || type == "r2") {
    ws <- "r2()"
    ws2 <- "R2"
  } else {
    ws <- "icc()"
    ws2 <- "ICC"
  }

  faminfo <- insight::model_info(x)
  vars <- .compute_variances(x, name_fun = ws, name_full = ws2, faminfo = faminfo)

  if (length(vars) == 1 && is.na(vars)) {
    return(NA)
  }

  # Calculate R2 values

  rsq.marginal <- vars$var.fixef / (vars$var.fixef + vars$var.ranef + vars$var.resid)
  rsq.conditional <- (vars$var.fixef + vars$var.ranef) / (vars$var.fixef + vars$var.ranef + vars$var.resid)

  names(rsq.marginal) <- "Marginal R2"
  names(rsq.conditional) <- "Conditional R2"


  # Calculate ICC values

  icc.adjusted <- vars$var.ranef / (vars$var.ranef + vars$var.resid)
  icc.conditional <- vars$var.ranef / (vars$var.fixef + vars$var.ranef + vars$var.resid)

  names(icc.adjusted) <-    "Adjusted ICC"
  names(icc.conditional) <- "Conditional ICC"


  if (is.null(type) || type == "r2") {
    var.measure <- structure(
      class = "sj_r2",
      list(rsq.marginal = rsq.marginal, rsq.conditional = rsq.conditional)
    )
  } else if (type == "all") {
    var.measure <- structure(
      class = "sj_iccr2",
      list(
        list(rsq.marginal = rsq.marginal, rsq.conditional = rsq.conditional),
        list(icc.adjusted = icc.adjusted, icc.conditional = icc.conditional)
      )
    )
  } else {
    var.measure <- structure(
      class = "sj_icc",
      list(icc.adjusted = icc.adjusted, icc.conditional = icc.conditional)
    )
  }


  # save variance information

  attr(var.measure, "var.fixef") <- vars$var.fixef
  attr(var.measure, "var.ranef") <- vars$var.ranef
  attr(var.measure, "var.disp") <- vars$var.disp
  attr(var.measure, "var.dist") <- vars$var.dist
  attr(var.measure, "var.resid") <- vars$var.resid

  attr(var.measure, "family") <- faminfo$family
  attr(var.measure, "link") <- faminfo$link_function
  attr(var.measure, "formula") <- insight::find_formula(x)

  # finally, save name of fitted model object. May be needed for
  # the 'se()' function, which accesses the global environment

  attr(var.measure, ".obj.name") <- obj.name

  var.measure
}


#' @importFrom stats nobs logLik
r2glm <- function(x, L.base) {
  L.full <- stats::logLik(x)
  D.full <- -2 * L.full

  D.base <- -2 * L.base
  G2 <- -2 * (L.base - L.full)

  if (inherits(x, c("vglm", "clm2")))
    n <- stats::nobs(x)
  else
    n <- attr(L.full, "nobs")

  Nagelkerke <- (1 - exp((D.full - D.base) / n)) / (1 - exp(-D.base / n))
  CoxSnell <- 1 - exp(-G2 / n)

  names(CoxSnell) <- "Cox & Snell's R-squared"
  names(Nagelkerke) <- "Nagelkerke's R-squared"

  # return results
  structure(class = "sj_r2", list(CoxSnell = CoxSnell, Nagelkerke = Nagelkerke))
}
