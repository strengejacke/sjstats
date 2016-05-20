#' @title Root Mean Squared Error (RMSE)
#' @name rmse
#' @description Compute root mean squared error  of fitted linear (mixed effects) models.
#'
#' @param fit Fitted linear model of class \code{\link{lm}},
#'          \code{\link[lme4]{merMod}} (lme4) or \code{\link[nlme]{lme}} (nlme).
#' @param normalized Logical, use \code{TRUE} if normalized rmse should be returned.
#'
#' @return The root mean squared error of \code{fit}; or the normalized
#'           root mean squared error of \code{fit} if \code{normalized = TRUE}.
#'
#' @seealso \code{\link{cv}} for the coefficient of variation, and
#'            \code{\link{rse}} for the residual standard error.
#'
#' @references \itemize{
#'              \item \href{http://en.wikipedia.org/wiki/Root-mean-square_deviation}{Wikipedia: RMSD}
#'              \item \href{http://www.theanalysisfactor.com/assessing-the-fit-of-regression-models/}{Grace-Martin K: Assessing the Fit of Regression Models}
#'             }
#'
#' @note The RMSE is the square root of the variance of the residuals and indicates
#'         the absolute fit of the model to the data (difference between observed data
#'         to model's predicted values). "RMSE can be interpreted as the standard
#'         deviation of the unexplained variance, and has the useful property
#'         of being in the same units as the response variable. Lower values
#'         of RMSE indicate better fit. RMSE is a good measure of how accurately
#'         the model predicts the response, and is the most important criterion
#'         for fit if the main purpose of the model is prediction."
#'         (Grace-Martin K: Assessing the Fit of Regression Models).
#'         \cr \cr
#'         The normalized RMSE is the proportion of the RMSE related to the
#'         range of the response variable. Hence, lower values indicate
#'         less residual variance.
#'
#' @examples
#' data(efc)
#' fit <- lm(barthtot ~ c160age + c12hour, data = efc)
#' rmse(fit)
#'
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' rmse(fit)
#'
#' # normalized RMSE
#' library(nlme)
#' fit <- lme(distance ~ age, data = Orthodont)
#' rmse(fit, normalized = TRUE)
#'
#' @importFrom stats residuals
#' @export
rmse <- function(fit, normalized = FALSE) {
  # compute rmse
  rmse_val <- sqrt(mean(stats::residuals(fit) ^ 2, na.rm = TRUE))

  # if normalized, divide by range of response
  if (normalized) {
    if (any(class(fit) == "lmerMod") || any(class(fit) == "merModLmerTest")) {
      # check for package availability
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' needed for this function to work. Please install it.", call. = FALSE)
      }
      resp <- lme4::getME(fit, "y")
    } else if (any(class(fit) == "lme")) {
      # check for package availability
      if (!requireNamespace("nlme", quietly = TRUE)) {
        stop("Package 'nlme' needed for this function to work. Please install it.", call. = FALSE)
      }
      resp <- unname(nlme::getResponse(fit))
    } else {
      resp <- fit$model[[1]]
    }
    rmse_val <- rmse_val / (max(resp, na.rm = T) - min(resp, na.rm = T))
  }
  rmse_val
}


#' @title Residual Standard Error (RSE)
#' @name rse
#' @description Compute the residual standard error of fitted linear (mixed effects) models.
#'
#' @param fit Fitted linear model of class \code{\link{lm}} or
#'          \code{\link[lme4]{merMod}} (\pkg{lme4}).
#'
#' @return The residual standard error of \code{fit}.
#'
#' @seealso \code{\link{cv}} for the coefficient of variation, and
#'            \code{\link{rmse}} for the root mean squared error.
#'
#' @note The residual standard error is the square root of the residual
#'        sum of squares divided by the residual degrees of freedom.
#'
#' @examples
#' data(efc)
#' fit <- lm(barthtot ~ c160age + c12hour, data = efc)
#' rse(fit)
#'
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' rse(fit)
#'
#' @importFrom stats residuals df.residual
#' @export
rse <- function(fit) {
  # Residual standard error
  sqrt(sum(stats::residuals(fit) ^ 2, na.rm = T) / stats::df.residual(fit))
}
