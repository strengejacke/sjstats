#' @title Eta-squared of fitted anova
#' @name eta_sq
#' @description Returns the eta-squared value for one-way-anovas.
#'
#' @param ... Fitted one-way-anova model or a dependent and grouping variable (see 'Examples').
#' @return The eta-squared value.
#'
#' @note Interpret eta-squared like r-squared or R-squared; a rule of thumb (Cohen):
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
#' fit <- aov(c12hour ~ as.factor(e42dep), data = efc)
#'
#' # print eta sqaured
#' eta_sq(fit)
#'
#' # grouping variable will be converted to factor autoamtically
#' eta_sq(efc$c12hour, efc$e42dep)
#'
#' @importFrom stats aov summary.lm
#' @export
eta_sq <- function(...) {
  # retrieve list of parameters
  input_list <- list(...)

  # check if fitted anova
  if (length(input_list) == 1 && inherits(input_list[[1]], "aov")) {
    # retrieve model
    fit <- input_list[[1]]
  } else if (length(input_list) == 2) {
    # retrieve variables
    depVar <- input_list[[1]]
    grpVar <- input_list[[2]]
    # convert to factor
    if (!is.factor(grpVar)) grpVar <- as.factor(grpVar)
    # fit anova
    fit <- stats::aov(depVar ~ grpVar)
  }
  # return eta squared
  stats::summary.lm(fit)$r.squared

  # 1 - var(fit$residuals, na.rm = T) / var(fit$model[,1], na.rm = T)
}
