#' @title Compute statistics for MCMC samples and Stan models
#' @name hdi
#'
#' @description \code{hdi()} computes the highest density interval for values from
#'   MCMC samples. \code{rope()} calculates the proportion of a posterior
#'   distribution that lies within a region of practical equivalence.
#'   \code{equi_test()} combines these two functions and performs a
#'   "HDI+ROPE decision rule" (Test for Practical Equivalence) (\cite{Kruschke 2018})
#'   to check whether parameter values should be accepted or rejected against
#'   the background of a formulated null hypothesis. \code{n_eff()} calculates
#'   the the number of effective samples (effective sample size). \code{mcse()}
#'   returns the Monte Carlo standard error. \code{mediation()} is a short
#'   summary for multivariate-response mediation-models.
#'
#' @param x A \code{stanreg}, \code{stanfit}, or \code{brmsfit} object. For
#'   \code{hdi()} and \code{rope()}, may also be a data frame or a vector
#'   of values from a probability distribution (e.g., posterior probabilities
#'   from MCMC sampling).
#' @param prob Vector of scalars between 0 and 1, indicating the mass within
#'   the credible interval that is to be estimated. See \code{\link{hdi}}.
#' @param rope Vector of length two, indicating the lower and upper limit of a
#'   range around zero, which indicates the region of practical equivalence.
#'   Values of the posterior distribution within this range are considered as
#'   being "practically equivalent to zero".
#' @param eff_size A scalar indicating the effect size that is used
#'   to calculate the limits of the ROPE for the test of practical equivalence.
#'   If not specified, an effect size of .1 is used, as suggested by
#'   \cite{Kruschke 2018} (see 'Details'). If \code{rope} is specified, this
#'   argument will be ignored.
#' @param trans Name of a function or character vector naming a function, used
#'   to apply transformations on the returned HDI-values resp.
#'   (for \code{rope()}) on the values of the posterior distribution, before
#'   calculating the rope based on the boundaries given in \code{rope}. Note
#'   that the values in \code{rope} are not transformed.
#' @param type For mixed effects models, specify the type of effects that should
#'   be returned. \code{type = "fixed"} returns fixed effects only,
#'   \code{type = "random"} the random effects and \code{type = "all"} returns
#'   both fixed and random effects.
#' @param out Character vector, indicating whether the results should be printed
#'    to console (\code{out = "txt"}) or as HTML-table in the viewer-pane
#'    (\code{out = "viewer"}) or browser (\code{out = "browser"}), of if the
#'    results should be plotted (\code{out = "plot"}, only applies to certain
#'    functions). May be abbreviated.
#' @param treatment Character, name of the treatment variable (or direct effect)
#'   in a (multivariate response) mediator-model. If missing, \code{mediation()}
#'   tries to find the treatment variable automatically, however, this may fail.
#' @param mediator Character, name of the mediator variable in a (multivariate
#'   response) mediator-model. If missing, \code{mediation()} tries to find the
#'   treatment variable automatically, however, this may fail.
#' @param typical The typical value that will represent the Bayesian point estimate.
#'   By default, the posterior median is returned. See \code{\link{typical_value}}
#'   for possible values for this argument.
#' @param ... Further arguments passed down to \code{equi_test()} when
#'   \code{plot = TRUE}:
#'   \itemize{
#'     \item{\code{colors}: }{
#'       Color of the density regions for the 95\% distribution of the posterior
#'       samples.
#'     }
#'     \item{\code{rope.color} and \code{rope.alpha}: }{
#'       Fill color and alpha-value of the ROPE (region of practical equivalence).
#'     }
#'     \item{\code{x.title}: }{
#'       Title for the x-axis of the plot.
#'     }
#'     \item{\code{legend.title}: }{
#'       Title for the plot legend.
#'     }
#'     \item{\code{labels}: }{
#'       Character vector of same length as terms plotted on the y-axis, to
#'       give axis labels user-defined labels.
#'     }
#'   }
#'
#'
#' @return For \code{hdi()}, if \code{x} is a vector, returns a vector of length
#'   two with the lower and upper limit of the HDI; if \code{x} is a
#'   \code{stanreg}, \code{stanfit} or \code{brmsfit} object, returns a
#'   tibble with lower and upper HDI-limits for each predictor. To distinguish
#'   multiple HDI values, column names for the HDI get a suffix when \code{prob}
#'   has more than one element.
#'   \cr \cr
#'   For \code{rope()}, returns a tibble with two columns: the proportion of
#'   values from \code{x} that are within and outside the boundaries of
#'   \code{rope}.
#'   \cr \cr
#'   \code{equi_test()} returns a tibble with a column \code{decision} that
#'   indicates whether or not a parameter value is accepted/rejected;
#'   \code{inside.rope}, which indicates the proportion of the whole posterior
#'   distribution that lies inside the ROPE (\emph{not} just the proportion of
#'   values from the 95\% HDI); and the lower and upper interval from the 95\%-HDI.
#'   \cr \cr
#'   \code{mcse()} and \code{n_eff()} return a tibble with two columns: one
#'   with the term names and one with the related statistic resp. effective
#'   sample size.
#'   \cr \cr
#'   \code{mediation()} returns a data frame with direct, indirect, mediator and
#'   total effect of a multivariate-response mediation-model, as well as the
#'   proportion mediated. The effect sizes are mean values of the posterior
#'   samples.
#'
#' @details
#'   \describe{
#'     \item{\strong{HDI}}{
#'       Computation for HDI is based on the code from Kruschke 2015, pp. 727f.
#'       For default sampling in Stan (4000 samples), the 90\% intervals for HDI are
#'       more stable than, for instance, 95\% intervals. An effective sample size
#'       (see \code{\link[brms]{nsamples}}) of at least 10.000 is recommended if
#'       95\% intervals should be computed (see Kruschke 2015, p. 183ff).
#'     }
#'     \item{\strong{MCSE}}{
#'       The Monte Carlo Standard Error is another useful measure of accuracy of
#'       the chains. It is defined as standard deviation of the chains divided by
#'       their effective sample size (the formula for \code{mcse()} is from
#'       Kruschke 2015, p. 187). The MCSE \dQuote{provides a quantitative suggestion
#'       of how big the estimation noise is}.
#'     }
#'     \item{\strong{Number of Effective Samples}}{
#'       The effective sample size divides the actual sample size by the amount
#'       of autocorrelation. The effective sample size is a measure of \dQuote{how
#'       much independent information there is in autocorrelated chains}, or:
#'       \dQuote{What would be the sample size of a completely non-autocorrelated chain
#'       that yielded the same information?} (\emph{Kruschke 2015, p182-3}).
#'       The ratio of effective number of samples and total number of samples
#'       (provided in \code{tidy_stan()}) ranges from 0 to 1, and should be close
#'       to 1. The closer this ratio comes to zero means that the chains may be
#'       inefficient, but possibly still okay.
#'     }
#'     \item{\strong{ROPE}}{
#'       There are no fixed rules to set the limits for the region of practical
#'       equivalence. However, there are some conventions described by
#'       \cite{Kruschke (2018)} how to specify the limits of the rope. One
#'       convention for linear models is to set the limits about .1 SD of the
#'       dependent variable around zero (i.e. \code{0 +/- .1 * sd(y)}), where
#'       .1 stands for half of a small effect size. Another, more conservative
#'       convention to set the ROPE limits is a range of half a standard
#'       deviation around zero (see \cite{Norman et al. 2003}), which indicates
#'       a clinical relevant effect (i.e. \code{0 +/- .25 * sd(y)} or even
#'       \code{0 +/- .5 * sd(y)})
#'     }
#'     \item{\strong{Test for Practical Equivalence}}{
#'       \code{equi_test()} computes the 95\%-HDI for \code{x} and checks if a
#'       model predictor's HDI lies completely outside, completely inside or
#'       partially inside the ROPE. If the HDI is completely outside the ROPE,
#'       the "null hypothesis" for this parameter is "rejected". If the ROPE
#'       completely covers the HDI, i.e. all most credible values of a parameter
#'       are inside the region of practical equivalence, the null hypothesis
#'       is accepted. Else, it's undecided whether to accept or reject the
#'       null hypothesis. In short, desirable results are low proportions inside
#'       the ROPE (the closer to zero the better) and the H0 should be rejected.
#'       \cr \cr
#'       If neither the \code{rope} nor \code{eff_size} argument are specified,
#'       the effect size will be set to 0.1 and the ROPE is \code{0 +/- .1 * sd(y)}
#'       for linear models. This is the suggested way to specify the ROPE limits
#'       according to \cite{Kruschke (2018)}. For models with binary outcome, there
#'       is no direct way to specify the effect size that defines the ROPE limits.
#'       Two examples from Kruschke suggest that a negligible change is about
#'       .05 on the logit-scale. In these cases, it is recommended to specify
#'       the \code{rope} argument, however, if not specified, the ROPE limits
#'       are calculated in this way: \code{0 +/- .1 * sd(intercept) / 4}. For
#'       all other models, \code{0 +/- .1 * sd(intercept)} is used to determine
#'       the ROPE limits.
#'       \cr \cr
#'       If \code{eff_size} is specified, but \code{rope} is not, then
#'       the same formulas apply, except that \code{.1} is replaced by the
#'       value in \code{eff_size}. If \code{rope} is specified, \code{eff_size}
#'       will be ignored. See also section \emph{ROPE} in 'Details'.
#'       \cr \cr
#'       The advantage of Bayesian testing for practical equivalence over
#'       classical frequentist null hypothesis significance testing is that
#'       discrete decisions are avoided, \dQuote{because such decisions encourage
#'       people to ignore the magnitude of the parameter value and its uncertainty}
#'       (\cite{Kruschke (2018)}).
#'     }
#'     \item{\strong{Mediation Analysis}}{
#'       \code{mediation()} returns a data frame with information on the
#'       \emph{direct effect} (mean value of posterior samples from \code{treatment}
#'       of the outcome model), \emph{mediator effect} (mean value of posterior
#'       samples from \code{mediator} of the outcome model), \emph{indirect effect}
#'       (mean value of the multiplication of the posterior samples from
#'       \code{mediator} of the outcome model and the posterior samples from
#'       \code{treatment} of the mediation model) and the total effect (mean
#'       value of sums of posterior samples used for the direct and indirect
#'       effect). The \emph{proportion mediated} is the indirect effect divided
#'       by the total effect.
#'       \cr \cr
#'       For all values, the 90\% HDIs are calculated by default. Use \code{prob}
#'       to calculate a different interval.
#'       \cr \cr
#'       The arguments \code{treatment} and \code{mediator} do not necessarily
#'       need to be specified. If missing, \code{mediation()} tries to find the
#'       treatment and mediator variable automatically. If this does not work,
#'       specify these variables.
#'     }
#'   }
#'
#' @note Since \code{equi_test()} computes 95\% HDI, a number of 10.000 samples
#'   produces more stable results (see Kruschke 2015, p. 183ff).
#'
#' @references
#'   Kruschke JK. Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. 2nd edition. Academic Press, 2015
#'   \cr \cr
#'   Kruschke JK. Rejecting or Accepting Parameter Values in Bayesian Estimation. Advances in Methods and Practices in Psychological Science. 2018; \doi{10.1177/2515245918771304}
#'   \cr \cr
#'   Norman GR, Sloan JA, Wyrwich KW. Interpretation of Changes in Health-related Quality of Life: The Remarkable Universality of Half a Standard Deviation. Medical Care. 2003;41: 582–592. \doi{10.1097/01.MLR.0000062554.74615.4C}
#'
#' @examples
#' \dontrun{
#' if (require("rstanarm")) {
#'   fit <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#'   hdi(fit)
#'
#'   # return multiple intervals
#'   hdi(fit, prob = c(.5, .7, .9))
#'
#'   # fit logistic regression model
#'   fit <- stan_glm(
#'     vs ~ wt + am,
#'     data = mtcars,
#'     family = binomial("logit"),
#'     chains = 1
#'   )
#'   # compute hdi, transform on "odds ratio scale"
#'   hdi(fit, trans = exp)
#'
#'   # compute rope, on scale of linear predictor. finds proportion
#'   # of posterior distribution values between -1 and 1.
#'   rope(fit, rope = c(-1, 1))
#'
#'   # compute rope, boundaries as "odds ratios". finds proportion of
#'   # posterior distribution values, which - after being exponentiated -
#'   # are between .8 and 1.25 (about -.22 and .22 on linear scale)
#'   rope(fit, rope = c(.8, 1.25), trans = exp)
#'
#'   # Test for Practical Equivalence
#'   equi_test(fit)
#'   equi_test(fit, out = "plot")
#' }}
#'
#' @export
hdi <- function(x, ...) {
  UseMethod("hdi")
}


#' @rdname hdi
#' @export
hdi.stanreg <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  hdi_worker(x = x, prob = prob, trans = trans, type = type)
}


#' @rdname hdi
#' @export
hdi.brmsfit <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all"), ...) {
  # check arguments
  type <- match.arg(type)

  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  hdi_worker(x = x, prob = prob, trans = trans, type = type)
}


#' @export
hdi.stanfit <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  hdi_worker(x = x, prob = prob, trans = trans, type = type)
}


#' @export
hdi.data.frame <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  hdi_worker(x = x, prob = prob, trans = trans, type = type)
}


#' @export
hdi.default <- function(x, prob = .9, trans = NULL, ...) {
  hdi_helper(x, prob, trans)
}


#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom purrr map_df
#' @importFrom sjmisc rotate_df
hdi_worker <- function(x, prob, trans, type) {
  dat <- purrr::map(
    prob,
    function(i) {

      out <- x %>%
        tibble::as_tibble() %>%
        purrr::map_df(~ hdi_helper(.x, i, trans)) %>%
        sjmisc::rotate_df() %>%
        tibble::rownames_to_column()

      colnames(out) <- c("term", "hdi.low", "hdi.high")
      out

    }) %>%
    dplyr::bind_cols() %>%
    dplyr::select(1, tidyselect::starts_with("hdi."))



  # for multiple HDIs, fix column names

  if (length(prob) > 1) {
    suffix <- prob %>%
      purrr::map(~ rep(.x, 2)) %>%
      purrr::flatten_dbl()

    colnames(dat)[2:ncol(dat)] <-
      sprintf(
        "%s_%s",
        rep(c("hdi.low", "hdi.high"), length(prob)),
        as.character(suffix)
      )
  }

  attr(dat, "prob") <- prob

  if (is_stan_model(x)) {
    # check if we need to remove random or fixed effects
    dat <- remove_effects_from_stan(dat, type, is.brms = inherits(x, "brmsfit"))
  }

  class(dat) <- c("sj_hdi", class(dat))
  dat
}


# based on Kruschke 2015, pp727f
#' @importFrom purrr map_dbl map_df
hdi_helper <- function(x, prob, trans) {
  x <- sort(x)
  ci.index <- ceiling(prob * length(x))
  nCIs <- length(x) - ci.index
  ci.width <- purrr::map_dbl(1:nCIs, ~ x[.x + ci.index] - x[.x])
  HDImin <- x[which.min(ci.width)]
  HDImax <- x[which.min(ci.width) + ci.index]

  # check if we have correct function
  if (!is.null(trans)) {
    trans <- match.fun(trans)
    HDImin <- trans(HDImin)
    HDImax <- trans(HDImax)
  }

  c(HDImin, HDImax)
}
