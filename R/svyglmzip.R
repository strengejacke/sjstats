utils::globalVariables("scaled.weights")

#' @title Survey-weighted zero-inflated Poisson model
#' @name svyglm.zip
#' @description \code{svyglm.zip()} is an extension to the \CRANpkg{survey}-package
#'                to fit survey-weighted zero-inflated Poisson models. It uses
#'                \code{\link[survey]{svymle}} to fit sampling-weighted
#'                maximum likelihood estimates, based on starting values provided
#'                by \code{\link[pscl]{zeroinfl}}.
#'
#'
#' @param formula An object of class \code{formula}, i.e. a symbolic description
#'          of the model to be fitted. See 'Details' in \code{\link[pscl]{zeroinfl}}.
#' @param design An object of class \code{\link[survey]{svydesign}}, providing
#'          a specification of the survey design.
#' @param ... Other arguments passed down to \code{\link[pscl]{zeroinfl}}.
#'
#' @return An object of class \code{\link[survey]{svymle}} and \code{svyglm.zip},
#'           with some additional information about the model.
#'
#' @details Code modified from https://notstatschat.rbind.io/2015/05/26/zero-inflated-poisson-from-complex-samples/.
#'
#' @examples
#' library(survey)
#' data(nhanes_sample)
#' set.seed(123)
#' nhanes_sample$malepartners <- rpois(nrow(nhanes_sample), 2)
#' nhanes_sample$malepartners[sample(1:2992, 400)] <- 0
#'
#' # create survey design
#' des <- svydesign(
#'   id = ~SDMVPSU,
#'   strat = ~SDMVSTRA,
#'   weights = ~WTINT2YR,
#'   nest = TRUE,
#'   data = nhanes_sample
#' )
#'
#' # fit negative binomial regression
#' fit <- svyglm.zip(
#'   malepartners ~ age + factor(RIDRETH1) | age + factor(RIDRETH1),
#'   des
#' )
#'
#' # print coefficients and standard errors
#' fit
#'
#' @importFrom insight find_formula
#' @importFrom stats weights update model.frame coef as.formula family
#' @export
svyglm.zip <- function(formula, design, ...) {
  # check if pkg survey is available
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package `survey` needed to for this function to work. Please install it.", call. = FALSE)
  }

  if (!requireNamespace("pscl", quietly = TRUE)) {
    stop("Package `pscl` needed to for this function to work. Please install it.", call. = FALSE)
  }


  # get design weights. we need to scale these weights for the glm.nb() function
  dw <- stats::weights(design)

  # update design with scaled weights
  design <- stats::update(design, scaled.weights = dw / mean(dw, na.rm = TRUE))

  # fit ZIP model, with scaled design weights
  mod <- pscl::zeroinfl(formula, data = stats::model.frame(design), weights = scaled.weights, ...)
  ff <- insight::find_formula(mod)

  # fit survey model, using maximum likelihood estimation
  svyfit <-
    survey::svymle(
      loglike = sjstats_loglik_zip,
      grad = sjstats_score_zip,
      design = design,
      formulas = list(eta = ff$conditional, logitp = ff$zero_inflated),
      start = stats::coef(mod),
      na.action = "na.omit"
    )

  # add additoinal information
  class(svyfit) <- c("svyglm.zip", class(svyfit))
  attr(svyfit, "zip.terms") <- all.vars(formula)
  attr(svyfit, "zip.formula") <- formula

  svyfit$deviance <- mod$deviance
  svyfit$df.residuals <- mod$df.residuals
  svyfit$df <- length(stats::coef(mod)) + 1
  svyfit$aic <- mod$aic

  svyfit
}


#' @importFrom stats dpois
# log-likelihood function used in "svymle()"
sjstats_loglik_zip <- function(y, eta, logitp) {
  mu <- exp(eta)
  p <- exp(logitp) / (1 + exp(logitp))
  log(p * (y == 0) + (1 - p) * stats::dpois(y, mu))
}


sjstats_dlogitp = function(y, eta, logitp) {
  mu <- exp(eta)
  p <- exp(logitp) / (1 + exp(logitp))
  dexpit <- p / (1 + p) ^ 2
  num <- dexpit * (y == 0) - dexpit * stats::dpois(y, mu)
  denom <- p * (y == 0) + (1 - p) * stats::dpois(y, mu)
  num / denom
}

# derivative
sjstats_deta_zip <- function(y, eta, logitp) {
  mu <- exp(eta)
  p <-  exp(logitp) / (1 + exp(logitp))
  dmutoy <- 0 * y
  dmutoy[y > 0] = exp(-mu[y > 0]) * mu[y > 0] ^ (y[y > 0] - 1) / factorial(y[y > 0] - 1)
  num = (1 - p) * (-stats::dpois(y, mu) + dmutoy)
  denom = p * (y == 0) + (1 - p) * stats::dpois(y, mu)
  num / denom
}

# score function, combines derivatives
sjstats_score_zip <- function(y, eta, logitp) {
  cbind(sjstats_deta_zip(y, eta, logitp), sjstats_dlogitp(y, eta, logitp))
}
