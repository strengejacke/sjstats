#' @title Effect size statistics for anova
#' @name eta_sq
#' @description Returns the (partial) eta-squared, omega-squared statistic
#'              or Cohen's F for all terms in an anovas. \code{anova_stats()} returns
#'              a tidy summary, including all these statistics and power for each term.
#'
#' @param model A fitted anova-model of class \code{aov} or \code{anova}. Other
#'              models are coerced to \code{\link[stats]{anova}}.
#' @param partial Logical, if \code{TRUE}, the partial eta-squared is returned.
#' @param digits Number of decimal points in the returned data frame.
#'
#' @return A numeric vector with the effect size statistics; for \code{anova_stats()},
#'         a tidy data frame with all these statistics.
#'
#' @references Levine TR, Hullett CR (2002): Eta Squared, Partial Eta Squared, and Misreporting of Effect Size in Communication Research (\href{https://www.msu.edu/~levinet/eta\%20squared\%20hcr.pdf}{pdf})
#'
#' @examples
#' # load sample data
#' data(efc)
#'
#' # fit linear model
#' fit <- aov(
#'   c12hour ~ as.factor(e42dep) + as.factor(c172code) + c160age,
#'   data = efc
#' )
#'
#' eta_sq(fit)
#' omega_sq(fit)
#' eta_sq(fit, partial = TRUE)
#'
#' anova_stats(car::Anova(fit, type = 2))
#'
#' @importFrom broom tidy
#' @importFrom purrr map_dbl
#' @importFrom stats anova
#' @export
eta_sq <- function(model, partial = FALSE) {
  if (partial)
    aov_stat(model, type = "peta")
  else
    aov_stat(model, type = "eta")
}


#' @rdname eta_sq
#' @export
omega_sq <- function(model) {
  aov_stat(model, type = "omega")
}


#' @rdname eta_sq
#' @export
cohens_f <- function(model) {
  aov_stat(model, type = "cohens.f")
}



#' @importFrom tibble tibble add_row add_column
#' @importFrom sjmisc add_columns
#' @importFrom broom tidy
#' @importFrom stats anova
#' @importFrom pwr pwr.f2.test
#' @rdname eta_sq
#' @export
anova_stats <- function(model, digits = 3) {
  # get tidy summary table
  aov.sum <- aov_stat_summary(model)

  # compute all model statstics
  etasq <- aov_stat_core(aov.sum, type = "eta")
  partial.etasq <- aov_stat_core(aov.sum, type = "peta")
  omegasq <- aov_stat_core(aov.sum, type = "omega")

  # compute power for each estimate
  cohens.f <- sqrt(partial.etasq / (1 - partial.etasq))

  # bind as data frame
  as <- tibble::tibble(etasq, partial.etasq, omegasq, cohens.f) %>%
    tibble::add_row(etasq = NA, partial.etasq = NA, omegasq = NA, cohens.f = NA) %>%
    sjmisc::add_columns(aov.sum)

  # get nr of terms
  nt <- nrow(as) - 1

  # finally, compute power
  power <- c(
    pwr::pwr.f2.test(u = as$df[1:nt], v = as$df[nrow(as)], f2 = as$cohens.f[1:nt] ^ 2)[["power"]],
    NA
  )

  tibble::add_column(as, power = power) %>%
    purrr::map_if(is.numeric, ~ round(.x, digits = digits)) %>%
    tibble::as_tibble()
}



#' @importFrom tibble has_name add_column
#' @importFrom dplyr mutate
#' @importFrom rlang .data
aov_stat <- function(model, type) {
  aov.sum <- aov_stat_summary(model)
  aov_stat_core(aov.sum, type)
}



aov_stat_summary <- function(model) {
  # check that model inherits from correct class
  # else, try to coerce to anova table
  if (!inherits(model, c("aov", "anova"))) model <- stats::anova(model)

  # get summary table
  aov.sum <- broom::tidy(model)

  # for car::Anova, the meansq-column might be missing, so add it manually
  if (!tibble::has_name(aov.sum, "meansq"))
    aov.sum <- tibble::add_column(aov.sum, meansq = aov.sum$sumsq / aov.sum$df, .after = "sumsq")

  aov.sum
}



aov_stat_core <- function(aov.sum, type) {
  # get mean squared of residuals
  meansq.resid <- aov.sum[["meansq"]][nrow(aov.sum)]
  # get total sum of squares
  ss.total <- sum(aov.sum[["sumsq"]])
  # get sum of squares of residuals
  ss.resid <- aov.sum[["sumsq"]][nrow(aov.sum)]

  # number of terms in model
  n_terms <- nrow(aov.sum) - 1

  if (type == "omega") {
    # compute omega squared for each model term
    aovstat <- purrr::map_dbl(1:n_terms, function(x) {
      ss.term <- aov.sum[["sumsq"]][x]
      df.term <- aov.sum[["df"]][x]
      (ss.term - df.term * meansq.resid) / (ss.total + meansq.resid)
    })
  } else if (type == "eta") {
    # compute eta squared for each model term
    aovstat <-
      purrr::map_dbl(1:n_terms, ~ aov.sum[["sumsq"]][.x] / sum(aov.sum[["sumsq"]]))
  } else if (type %in% c("cohens.f", "peta")) {
    # compute partial eta squared for each model term
    aovstat <-
      purrr::map_dbl(1:n_terms, ~ aov.sum[["sumsq"]][.x] / (aov.sum[["sumsq"]][.x] + ss.resid))
  }

  # compute Cohen's F
  if (type == "cohens.f") aovstat <- sqrt(aovstat / (1 - aovstat))

  # give values names of terms
  names(aovstat) <- aov.sum[["term"]][1:n_terms]

  aovstat
}
