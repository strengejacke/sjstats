#' @title Effect size statistics for anova
#' @name eta_sq
#' @description Returns the (partial) eta-squared or omega-squared statistic
#'              for all terms in an anovas. \code{anova_stats()} returns
#'              a tidy summary, including all three statistics.
#'
#' @param model A fitted anova-model of class \code{aov} or \code{anova}. Other
#'              models are coerced to \code{\link[stats]{anova}}.
#' @param partial Logical, if \code{TRUE}, the partial eta-squared is returned.
#'
#' @return A numeric vector with the effect size statistics; for \code{anova_stats()},
#'         a tidy data frame with all three statistics.
#'
#' @note Interpret eta-squared like r-squared; a rule of thumb (Cohen):
#'         \itemize{
#'          \item .02 ~ small
#'          \item .13 ~ medium
#'          \item .26 ~ large
#'         }
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



#' @importFrom tibble tibble add_row
#' @importFrom sjmisc add_columns
#' @importFrom broom tidy
#' @importFrom stats anova
#' @rdname eta_sq
#' @export
anova_stats <- function(model) {
  # check that model inherits from correct class
  # else, try to coerce to anova table
  if (!inherits(model, c("aov", "anova"))) model <- stats::anova(model)

  # compute all model statstics
  etasq <- aov_stat(model, type = "eta")
  partial.etasq <- aov_stat(model, type = "peta")
  omegasq <- aov_stat(model, type = "omega")

  # bind as data frame
  tibble::tibble(etasq, partial.etasq, omegasq) %>%
    tibble::add_row(etasq = NA, partial.etasq = NA, omegasq = NA) %>%
    sjmisc::add_columns(broom::tidy(model))
}



#' @importFrom tibble has_name
#' @importFrom dplyr mutate
#' @importFrom rlang .data
aov_stat <- function(model, type) {
  # check that model inherits from correct class
  # else, try to coerce to anova table
  if (!inherits(model, c("aov", "anova"))) model <- stats::anova(model)

  # get summary table and number of terms in model
  aov.sum <- broom::tidy(model)
  n_terms <- nrow(aov.sum) - 1

  # for car::Anova, the meansq-column might be missing, so add it manually
  if (!tibble::has_name(aov.sum, "meansq"))
    aov.sum <- dplyr::mutate(aov.sum, meansq = .data$sumsq / .data$df)

  # get mean squared of residuals
  meansq.resid <- aov.sum[["meansq"]][nrow(aov.sum)]
  # get total sum of squares
  ss.total <- sum(aov.sum[["sumsq"]])
  # get sum of squares of residuals
  ss.resid <- aov.sum[["sumsq"]][nrow(aov.sum)]


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
  } else {
    # compute partial eta squared for each model term
    aovstat <-
      purrr::map_dbl(1:n_terms, ~ aov.sum[["sumsq"]][.x] / (aov.sum[["sumsq"]][.x] + ss.resid))
  }


  # give values names of terms
  names(aovstat) <- aov.sum[["term"]][1:n_terms]

  aovstat
}
