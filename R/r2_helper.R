# Compute variance associated with a random-effects term
# (Johnson 2014)

get_ranef_variance <- function(terms, x, vals) {
  sum(sapply(
    vals$vc[terms],
    function(Sigma) {
      rn <- rownames(Sigma)

      if (!is.null(rn)) {
        valid <- rownames(Sigma) %in% colnames(vals$X)
        if (!all(valid)) {
          rn <- rn[valid]
          Sigma <- Sigma[valid, valid]
        }
      }

      Z <- vals$X[, rn, drop = FALSE]
      # Z <- vals$X[, rownames(Sigma), drop = FALSE]
      Z.m <- Z %*% Sigma
      sum(diag(crossprod(Z.m, Z))) / stats::nobs(x)
    }))
}


# get fixed effects variance

#' @importFrom stats var
get_fixef_variance <- function(vals) {
  with(vals, stats::var(as.vector(beta %*% t(X))))
}


# Get residual (distribution specific) variance from random effects

get_residual_variance <- function(x, var.cor, faminfo, type) {

  sig <- attr(var.cor, "sc")
  if (is.null(sig)) sig <- 1

  if (faminfo$is_linear) {
    residual.variance <- sig^2
  } else {

    if (faminfo$is_bin) {
      residual.variance <- switch(
        faminfo$link.fun,
        logit = pi^2 / 3,
        probit = 1,
        badlink(faminfo$link.fun, faminfo$family)
      )
    } else if (faminfo$is_pois) {
      residual.variance <- switch(
        faminfo$link.fun,
        log = logVarDist(x, null_model(x), faminfo, sig, type = type),
        sqrt = 0.25,
        badlink(faminfo$link.fun, faminfo$family)
      )
    } else if (faminfo$family == "beta") {
      residual.variance <- switch(
        faminfo$link.fun,
        logit = logVarDist(x, null_model(x), faminfo, sig, type = type),
        badlink(faminfo$link.fun, faminfo$family)
      )
    }
  }

  residual.variance
}


# get dispersion-specific variance

get_disp_variance <- function(x, vals, faminfo, obs.terms) {
  if (faminfo$is_linear) {
    0
  } else {
    if (length(obs.terms) == 0 )
      0
    else
      get_ranef_variance(obs.terms, x = x, vals = vals)
  }
}


# helper-function, telling user if model is supported or not

badlink <- function(link, family) {
  warning(sprintf("Model link '%s' is not yet supported for the %s distribution.", link, family), call. = FALSE)
  return(NA)
}


# glmmTMB returns a list of model information, one for conditional and one
# for zero-inflated part, so here we "unlist" it

collapse_cond <- function(fit) {
  if (is.list(fit) && "cond" %in% names(fit))
    fit[["cond"]]
  else
    fit
}


# Generate null model (intercept and random effects only, no fixed effects)

null_model <- function(x) {

  ## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q4/023013.html
  ## FIXME: deparse is a *little* dangerous
  rterms <- paste0("(", sapply(lme4::findbars(stats::formula(x)), deparse), ")")
  nullform <- stats::reformulate(rterms, response = ".")
  null.model <- stats::update(x, nullform)

  ## Get the fixed effects of the null model
  unname(collapse_cond(lme4::fixef(null.model)))
}


# helper function, compute distributional variance for beta-family

beta_variance <- function(mu, phi) {
  mu * (1 - mu) / (1 + phi)
}


# distributional variance for different models

logVarDist <- function(x, null.fixef, faminfo, sig, type) {
  ## in general want log(1+var(x)/mu^2)
  mu <- exp(null.fixef)
  if (mu < 6)
    warning(sprintf("mu of %0.1f is too close to zero, estimate of %s may be unreliable.\n", mu, type), call. = FALSE)

  vv <- switch(
    faminfo$family,
    poisson = stats::family(x)$variance(mu),
    truncated_poisson = stats::family(x)$variance(sig),
    beta = beta_variance(mu, sig),
    genpois = ,
    nbinom1 = ,
    nbinom2 = stats::family(x)$variance(mu, sig),

    if (inherits(x,"merMod"))
      mu * (1 + mu / lme4::getME(x, "glmer.nb.theta"))
    else
      mu * (1 + mu / x$theta)
  )

  cvsquared <- vv / mu^2
  log1p(cvsquared)
}
