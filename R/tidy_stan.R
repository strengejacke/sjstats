#' @title Tidy summary output for stan models
#' @name tidy_stan
#'
#' @description Returns a tidy summary output for stan models.
#'
#' @param x A \code{stanreg}, \code{stanfit} or \code{brmsfit} object.
#' @param probs Vector of scalars between 0 and 1, indicating the mass within
#'        the credible interval that is to be estimated. See \code{\link{hdi}}.
#' @param typical The typical value that will represent the Bayesian point estimate.
#'        By default, the posterior median is returned. See \code{\link{typical_value}}
#'        for possible values for this argument.
#' @param trans Name of a function or character vector naming a function, used
#'        to apply transformations on the estimate and HDI-values. The
#'        values for standard errors are \emph{not} transformed!
#' @param digits Amount of digits to round numerical values in the output.
#'
#' @inheritParams hdi
#'
#' @return A tidy data frame, summarizing \code{x}, with consistent column names.
#'         To distinguish multiple HDI values, column names for the HDI get a suffix
#'         when \code{probs} has more than one element.
#'
#' @details The returned data frame gives information on the Bayesian point
#'          estimate (column \emph{estimate}, which is by default the posterior
#'          median; other statistics are also possible, see \code{typical}), the
#'          standard error (which are actually \emph{median absolute deviations}),
#'          the HDI, the ratio of effective numbers of samples, \emph{n_eff},
#'          (i.e. effective number of samples divided by total number of samples)
#'          and Rhat statistics.
#'          \cr \cr
#'          The ratio of effective number of samples ranges from 0 to 1,
#'          and should be close to 1. The closer this ratio comes to zero means
#'          that the chains may be inefficient, but possibly still okay.
#'          \cr \cr
#'          When Rhat is above 1, it usually indicates that the chain has not
#'          yet converged, indicating that the drawn samples might not be
#'          trustworthy. Drawing more iteration may solve this issue.
#'          \cr \cr
#'          Computation for HDI is based on the code from Kruschke 2015, pp. 727f.
#'
#' @seealso \code{\link{hdi}}
#'
#' @references Kruschke JK. \emph{Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan.} 2nd edition. Academic Press, 2015
#' \cr \cr
#' Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB. \emph{Bayesian data analysis.} 3rd ed. Boca Raton: Chapman & Hall/CRC, 2013
#' \cr \cr
#' Gelman A, Rubin DB. \emph{Inference from iterative simulation using multiple sequences.} Statistical Science 1992;7: 457–511
#' \cr \cr
#' McElreath R. \emph{Statistical Rethinking. A Bayesian Course with Examples in R and Stan.} Chapman and Hall, 2015
#'
#' @examples
#' \dontrun{
#' if (require("rstanarm")) {
#'   fit <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#'   tidy_stan(fit)
#'   tidy_stan(fit, probs = c(.89, .5))
#' }}
#'
#' @importFrom purrr map flatten_dbl map_dbl modify_if
#' @importFrom dplyr bind_cols select mutate
#' @importFrom tidyselect starts_with
#' @importFrom tibble add_column
#' @importFrom stats mad
#' @importFrom bayesplot rhat neff_ratio
#' @export
tidy_stan <- function(x, probs = .89, typical = "median", trans = NULL, type = c("fixed", "random", "all"), digits = 3) {

  # only works for rstanarm-models

  if (!inherits(x, c("stanreg", "stanfit", "brmsfit")))
    stop("`x` needs to be a stanreg- or brmsfit-object.", call. = F)


  # check arguments
  type <- match.arg(type)


  # get data frame

  mod.dat <- as.data.frame(x)
  brmsfit.removers <- NULL

  # for brmsfit models, we need to remove some columns here to
  # match data rows later

  if (inherits(x, "brmsfit")) {
    re.sd <- tidyselect::starts_with("sd_", vars = colnames(mod.dat))
    re.cor <- tidyselect::starts_with("cor_", vars = colnames(mod.dat))

    brmsfit.removers <- unique(c(re.sd, re.cor))

    if (!sjmisc::is_empty(brmsfit.removers))
      mod.dat <- dplyr::select(mod.dat, !! -brmsfit.removers)
  }


  # compute HDI

  out <- purrr::map(probs, ~ hdi(x, prob = .x, trans = trans, type = "all")) %>%
    dplyr::bind_cols() %>%
    dplyr::select(1, tidyselect::starts_with("hdi."))


  # for multiple HDIs, fix column names

  if (length(probs) > 1) {
    suffix <- probs %>%
      purrr::map(~ rep(.x, length(probs))) %>%
      purrr::flatten_dbl()

    colnames(out)[2:ncol(out)] <-
      sprintf(
        "%s_%s",
        rep(c("hdi.low", "hdi.high"), length(probs)),
        as.character(suffix)
      )
  }


  # compute additional statistics, like point estimate, standard errors etc.

  if (sjmisc::is_empty(brmsfit.removers)) {
    nr <- bayesplot::neff_ratio(x)[1:nrow(out)]
    rh <- bayesplot::rhat(x)[1:nrow(out)]
    se <- dplyr::pull(mcse(x), "mcse")[1:nrow(out)]
  } else {
    nr <- bayesplot::neff_ratio(x)[-brmsfit.removers]
    rh <- bayesplot::rhat(x)[-brmsfit.removers]
    se <- dplyr::pull(mcse(x), "mcse")[1:nrow(out)]
  }


  out <- out %>%
    tibble::add_column(
      estimate = purrr::map_dbl(mod.dat, ~ typical_value(.x, fun = typical)),
      std.error = purrr::map_dbl(mod.dat, stats::mad),
      .after = 1
    ) %>%
    dplyr::mutate(
      n_eff = nr,
      Rhat = rh,
      mcse = se
    )


  # transform estimate, if requested

  if (!is.null(trans)) {
    trans <- match.fun(trans)
    out$estimate <- trans(out$estimate)
  }


  # check if we need to remove random or fixed effects

  out <- remove_effects_from_stan(out, type, is.brms = inherits(x, "brmsfit"))


  # round values

  purrr::modify_if(out, is.numeric, ~ round(.x, digits = digits))
}


#' @importFrom tidyselect starts_with ends_with
#' @importFrom dplyr slice
#' @importFrom tibble as_tibble
#' @importFrom sjmisc is_empty
remove_effects_from_stan <- function(out, type, is.brms) {

  # brmsfit-objects also include sd and cor for mixed
  # effecs models, so remove these here

  if (is.brms) out <- brms_clean(out)

  # if user wants all terms, return data here

  if (type == "all") return(tibble::as_tibble(out))


  # check if we need to remove random or fixed effects
  # therefor, find random effect parts first

  re <- tidyselect::starts_with("b[", vars = out$term)
  re.s <- tidyselect::starts_with("Sigma[", vars = out$term)
  re.i <- intersect(
    tidyselect::starts_with("r_", vars = out$term),
    tidyselect::ends_with(".", vars = out$term)
  )

  removers <- unique(c(re, re.s, re.i))

  if (!sjmisc::is_empty(removers)) {
    if (type == "fixed") {
      # remove all random effects
      out <- dplyr::slice(out, !! -removers)
    } else if (type == "random") {
      # remove all fixed effects
      out <- dplyr::slice(out, !! removers)
    }
  }


  tibble::as_tibble(out)
}


#' @importFrom tidyselect starts_with ends_with
#' @importFrom dplyr slice
#' @importFrom sjmisc is_empty
brms_clean <- function(out) {

  # brmsfit-objects also include sd and cor for mixed
  # effecs models, so remove these here

  re.sd <- tidyselect::starts_with("sd_", vars = out$term)
  re.cor <- tidyselect::starts_with("cor_", vars = out$term)

  removers <- unique(c(re.sd, re.cor))

  if (!sjmisc::is_empty(removers))
    out <- dplyr::slice(out, !! -removers)

  out
}
