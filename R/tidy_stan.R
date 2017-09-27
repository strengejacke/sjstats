#' @title Tidy summary output for stan models
#' @name tidy_stan
#'
#' @description Returns a tidy summary output for stan models.
#'
#' @param x A \code{stanreg}, \code{stanfit} or \code{brmsfit} object.
#' @param probs Vector of scalars between 0 and 1, indicating the mass within
#'        the credible interval that is to be estimated.
#' @param typical The typical value that will represent the Bayesian point estimate.
#'        By default, the posterior median is returned. See \code{\link{typical_value}}
#'        for possible values for this argument.
#' @param trans Name of a function or character vector naming a function, used
#'        to apply transformations on the estimate and HDI-values. The
#'        values for standard errors are \emph{not} transformed!
#' @param digits Amount of digits to round numerical values in the output.
#'
#' @return A tidy data frame, summarizing \code{x}, giving information on the
#'         Bayesian point estimate, standard error (which are effectively
#'         \emph{median absolute deviations}), HDI, ratio of effective
#'         sample size numbers and Rhat statistics.
#'
#' @details Computation for HDI is based on the code from Kruschke 2015, pp. 727f.
#'
#' @seealso \code{\link{hdi}}
#'
#' @references Kruschke JK. Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. 2nd edition. Academic Press, 2015
#'
#' @examples
#' if (require("rstanarm")) {
#'   fit <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#'   tidy_stan(fit)
#'   tidy_stan(fit, probs = c(.89, .5))
#' }
#'
#' @importFrom purrr map flatten_dbl map_dbl modify_if
#' @importFrom dplyr bind_cols select starts_with mutate pull
#' @importFrom tibble add_column
#' @importFrom stats mad
#' @export
tidy_stan <- function(x, probs = .89, typical = "median", trans = NULL, digits = 3) {

  # only works for rstanarm-models

  if (!inherits(x, c("stanreg", "stanfit", "brmsfit")))
    stop("`x` needs to be a stanreg-object.", call. = F)


  # compute HDI

  out <- purrr::map(probs, ~ hdi(x, prob = .x, trans = trans)) %>%
    dplyr::bind_cols() %>%
    dplyr::select(1, dplyr::starts_with("hdi."))


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

  out <- out %>%
    tibble::add_column(
      estimate = purrr::map_dbl(as.data.frame(x), ~ typical_value(.x, fun = typical)),
      .after = 1
    ) %>%
    tibble::add_column(
      std.error = purrr::map_dbl(as.data.frame(x), stats::mad),
      .after = 2
    ) %>%
    dplyr::mutate(
      n_eff = neff_ratio(x)[1:nrow(out)],
      Rhat = rhat(x)[1:nrow(out)]
    )


  # transform estimate, if requested

  if (!is.null(trans)) {
    trans <- match.fun(trans)
    out$estimate <- trans(out$estimate)
  }


  # round values

  purrr::modify_if(out, is.numeric, ~ round(.x, digits = digits))
}
