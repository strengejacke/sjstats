#' @title Tidy summary output for stan models
#' @name tidy_stan
#'
#' @description Returns a tidy summary output for stan models.
#'
#' @param x A \code{stanreg}, \code{stanfit} or \code{brmsfit} object.
#' @param trans Name of a function or character vector naming a function, used
#'        to apply transformations on the estimates and uncertainty intervals. The
#'        values for standard errors are \emph{not} transformed! If \code{trans}
#'        is not \code{NULL}, \emph{credible intervals} instead of \emph{HDI}
#'        are computed, due to the possible asymmetry of the HDI.
#' @param digits Amount of digits to round numerical values in the output.
#' @param prob Vector of scalars between 0 and 1, indicating the mass within
#'   the credible interval that is to be estimated.
#' @param typical The typical value that will represent the Bayesian point estimate.
#'   By default, the posterior median is returned. See \code{\link[sjmisc]{typical_value}}
#'   for possible values for this argument.
#' @inheritParams bayestestR::hdi
#'
#' @return A data frame, summarizing \code{x}, with consistent column names.
#'         To distinguish multiple HDI values, column names for the HDI get a suffix
#'         when \code{prob} has more than one element.
#'
#' @details The returned data frame has an additonal class-attribute,
#'    \code{tidy_stan}, to pass the result to its own \code{print()}-method.
#'    The \code{print()}-method creates a cleaner output, especially for multilevel,
#'    zero-inflated or multivariate response models, where - for instance -
#'    the conditional part of a model is printed separately from the zero-inflated
#'    part, or random and fixed effects are printed separately.
#'    \cr \cr
#'    The returned data frame gives information on:
#'    \itemize{
#'      \item{The Bayesian point estimate (column \emph{estimate}, which is by
#'            default the posterior median; other statistics are also possible,
#'            see argument \code{typical}).}
#'      \item{
#'        The standard error (which is actually the \emph{median absolute deviation}).
#'      }
#'      \item{
#'        The HDI. Computation for HDI is based on the
#'        code from Kruschke 2015, pp. 727f.
#'      }
#'      \item{
#'        The Probability of Direction (pd), which is an index for "effect significance"
#'        (see \cite{Makowski et al. 2019}). A value of 95\% or higher indicates a
#'        "significant" (i.e. statistically clear) effect.
#'      }
#'      \item{
#'        The effective numbers of samples, \emph{ESS}.
#'      }
#'      \item{
#'        The Rhat statistics. When Rhat is above 1, it usually indicates that
#'        the chain has not yet converged, indicating that the drawn samples
#'        might not be trustworthy. Drawing more iteration may solve this issue.
#'      }
#'      \item{
#'        The Monte Carlo standard error (see \code{\link{mcse}}). It is defined
#'        as standard deviation of the chains divided by their effective sample
#'        size and \dQuote{provides a quantitative suggestion of how big the
#'        estimation noise is} (\emph{Kruschke 2015, p.187}).
#'      }
#'    }
#'
#' @references Kruschke JK. \emph{Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan} 2nd edition. Academic Press, 2015
#' \cr \cr
#' Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB. \emph{Bayesian data analysis} 3rd ed. Boca Raton: Chapman and Hall/CRC, 2013
#' \cr \cr
#' Gelman A, Rubin DB. \emph{Inference from iterative simulation using multiple sequences} Statistical Science 1992;7: 457-511
#' \cr \cr
#' Makowski D, Ben-Shachar MS, Lüdecke D. bayestestR: Describing Effects and their Uncertainty, Existence and Significance within the Bayesian Framework. Journal of Open Source Software 2019;4:1541. \doi{10.21105/joss.01541}
#' \cr \cr
#' McElreath R. \emph{Statistical Rethinking. A Bayesian Course with Examples in R and Stan} Chapman and Hall, 2015
#'
#' @examples
#' \dontrun{
#' if (require("rstanarm")) {
#'   fit <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#'   tidy_stan(fit)
#'   tidy_stan(fit, prob = c(.89, .5))
#' }}
#'
#' @importFrom purrr map_dbl
#' @importFrom stats mad formula sd
#' @importFrom sjmisc is_empty trim seq_col typical_value
#' @importFrom insight model_info get_parameters  is_multivariate print_parameters
#' @importFrom bayestestR hdi ci effective_sample mcse pd
#' @export
tidy_stan <- function(x, prob = .89, typical = "median", trans = NULL, effects = c("all", "fixed", "random"), component = c("all", "conditional", "zero_inflated", "zi"), digits = 2) {

  # only works for rstanarm- or brms-models
  if (!inherits(x, c("stanreg", "stanfit", "brmsfit")))
    stop("`x` needs to be a stanreg- or brmsfit-object.", call. = F)

  # check arguments
  effects <- match.arg(effects)
  component <- match.arg(component)

  # get transformaton function
  if (!is.null(trans)) trans <- match.fun(trans)

  # family info
  faminfo <- insight::model_info(x)

  if (insight::is_multivariate(x)) {
    faminfo <- faminfo[[1]]
  }

  # get parameters ----
  out.pars <- insight::get_parameters(x, effects = effects, component = component, parameters = "^(?!prior)")

  # compute HDI / ci ----
  if (!is.null(trans)) {
    out.hdi <- bayestestR::ci(x, ci = prob, effects = effects, component = component)
    out.hdi$CI_low <- trans(out.hdi$CI_low)
    out.hdi$CI_high <- trans(out.hdi$CI_high)
    is_hdi <- FALSE
  } else {
    out.hdi <- bayestestR::hdi(x, ci = prob, effects = effects, component = component)
    is_hdi <- TRUE
  }

  out.hdi$CI_low <- sprintf("%.*f", digits, out.hdi$CI_low)
  out.hdi$CI_high <- sprintf("%.*f", digits, out.hdi$CI_high)

  # transform data frame for multiple ci-levels
  if (length(unique(out.hdi$CI)) > 1) {
    hdi_list <- lapply(
      split(out.hdi, out.hdi$CI, drop = FALSE),
      function(i) {
        .rename_ci_column(i, ifelse(is_hdi, "HDI", "CI"))
      })
    hdi_frame <- Reduce(function(x, y) merge(x, y, all.y = TRUE, by = "Parameter"), hdi_list)
    to_remove <- string_starts_with("CI.", colnames(hdi_frame))
    to_remove <- c(to_remove, string_starts_with("Component", colnames(hdi_frame)))
    to_remove <- c(to_remove, string_starts_with("Group", colnames(hdi_frame)))
    to_remove <- c(to_remove, string_starts_with("Effects", colnames(hdi_frame)))
    to_remove <- c(to_remove, string_starts_with("Response", colnames(hdi_frame)))
    to_remove <- c(to_remove, which(colnames(hdi_frame) == "CI"))
    out.hdi <- hdi_frame[, -to_remove, drop = FALSE]
  } else {
    out.hdi <- .rename_ci_column(out.hdi, ifelse(is_hdi, "HDI", "CI"))
  }

  remove <- which(colnames(out.hdi) == "CI" | grepl("^CI_", colnames(out.hdi)))
  out.hdi <- out.hdi[, -remove]

  # compute effective sample size ESS ----
  out.ess <- bayestestR::effective_sample(x, effects = effects, component = component)


  # compute MCSE ----
  out.mcse <- bayestestR::mcse(x, effects = effects, component = component)


  # compute Probability of Direction ----
  out.pd <- bayestestR::pd(x, effects = effects, component = component, method = "direct")


  # compute RHat ----
  out.rhat <- .rhat(x)


  # transform estimate, if requested
  if (!is.null(trans)) {
    all.cols <- sjmisc::seq_col(out.pars)
    simp.pars <- string_starts_with("simo_mo", colnames(out.pars))
    if (!sjmisc::is_empty(simp.pars)) all.cols <- all.cols[-simp.pars]
    for (i in all.cols) out.pars[[i]] <- trans(out.pars[[i]])
  }

  se.fun <- switch(
    typical,
    "median" = stats::mad,
    stats::sd
  )

  out.parameters <- data.frame(
    Parameter = colnames(out.pars),
    Estimate = purrr::map_dbl(out.pars, ~ sjmisc::typical_value(.x, fun = typical)),
    Std.Error = purrr::map_dbl(out.pars, se.fun),
    stringsAsFactors = FALSE
  )

  out <- insight::print_parameters(x, out.parameters, out.hdi, out.pd, out.ess, out.rhat, out.mcse)
  class(out) <- c("tidy_stan", class(out))
  attr(out, "digits") <- digits

  out
}



.rename_ci_column <- function(x, col_name) {
  x$CI_low <- format(x$CI_low, justify = "right")
  x$CI_high <- format(x$CI_high, justify = "right")
  x$.ci <- sprintf("[%s, %s]", x$CI_low, x$CI_high)
  colnames(x)[ncol(x)] <- sprintf("%s(%g%%)", col_name, x$CI[1])
  x
}



.rhat <- function(x) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Package `rstan` is required. Please install it first.", call. = FALSE)

  if (inherits(x, "brmsfit")) x <- x$fit

  if (inherits(x, "stanfit")) {
    s <- rstan::summary(x)$summary
  } else if (inherits(x, "stanreg")) {
    s <- summary(x, pars = NULL, regex_pars = NULL)
  }

  data.frame(
    Parameter = make.names(rownames(s)),
    Rhat = s[, "Rhat"],
    stringsAsFactors = FALSE
  )
}
