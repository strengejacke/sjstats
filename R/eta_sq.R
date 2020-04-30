#' @title Effect size statistics for anova
#' @name eta_sq
#' @description Returns the (partial) eta-squared, (partial) omega-squared,
#'   epsilon-squared statistic or Cohen's F for all terms in an anovas.
#'   \code{anova_stats()} returns a tidy summary, including all these statistics
#'   and power for each term.
#'
#' @param model A fitted anova-model of class \code{aov} or \code{anova}. Other
#'   models are coerced to \code{\link[stats]{anova}}.
#' @param partial Logical, if \code{TRUE}, the partial eta-squared is returned.
#' @param ci.lvl Scalar between 0 and 1. If not \code{NULL}, returns a data
#'   frame with effect sizes including lower and upper confidence intervals.
#' @param digits Amount of digits for returned values.
#'
#' @return A data frame with the term name(s) and effect size statistics; if
#'   \code{ci.lvl} is not \code{NULL}, a data frame including lower and
#'   upper confidence intervals is returned. For \code{anova_stats()}, a tidy
#'   data frame with all statistics is returned (excluding confidence intervals).
#'
#' @details See details in \code{\link[effectsize]{eta_squared}}.
#'
#' @references Levine TR, Hullett CR (2002): Eta Squared, Partial Eta Squared, and Misreporting of Effect Size in Communication Research (\href{https://www.msu.edu/~levinet/eta\%20squared\%20hcr.pdf}{pdf})
#'   \cr \cr
#'   Tippey K, Longnecker MT (2016): An Ad Hoc Method for Computing Pseudo-Effect Size for Mixed Model. (\href{http://www.scsug.org/wp-content/uploads/2016/11/Ad-Hoc-Method-for-Computing-Effect-Size-for-Mixed-Models_PROCEEDINGS-UPDATE-1.pdf}{pdf})
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
#' eta_sq(fit, partial = TRUE, ci.lvl = .8)
#'
#' anova_stats(car::Anova(fit, type = 2))
#' @importFrom effectsize eta_squared
#' @export
eta_sq <- function(model, partial = FALSE, ci.lvl = NULL) {
  out <- effectsize::eta_squared(model, partial = partial, ci = ci.lvl)

  if (isTRUE(partial)) {
    cname <- "partial.etasq"
  } else {
    cname <- "etasq"
  }

  out <- .fix_column_names(out, cname)
  class(out) <- c("sj_anova_stat", class(out))
  out
}




# helper -------------------


.fix_column_names <- function(out, cname) {
  colnames(out)[1] <- "term"
  colnames(out)[2] <- cname
  if ("CI_low" %in% colnames(out)) colnames(out)[which(colnames(out) == "CI_low")] <- "conf.low"
  if ("CI_high" %in% colnames(out)) colnames(out)[which(colnames(out) == "CI_high")] <- "conf.high"
  if ("CI" %in% colnames(out)) out$CI <- NULL
  out
}

