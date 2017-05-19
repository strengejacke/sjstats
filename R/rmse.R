#' @title Compute model quality
#' @name rmse
#' @description Compute root mean squared error, residual standard error or
#'              mean square error of fitted linear (mixed effects) models.
#'
#' @param fit Fitted linear model of class \code{\link{lm}},
#'          \code{\link[lme4]{merMod}} (\pkg{lme4}) or \code{\link[nlme]{lme}} (\pkg{nlme}).
#' @param normalized Logical, use \code{TRUE} if normalized rmse should be returned.
#'
#' @seealso \code{\link{r2}} for R-squared or pseude-R-squared values, and
#'            \code{\link{cv}} for the coefficient of variation.
#'
#' @references \href{http://www.theanalysisfactor.com/assessing-the-fit-of-regression-models/}{Grace-Martin K: Assessing the Fit of Regression Models}
#'
#' @note \describe{
#'         \item{\strong{Root Mean Square Error}}{
#'         The RMSE is the square root of the variance of the residuals and indicates
#'         the absolute fit of the model to the data (difference between observed data
#'         to model's predicted values). \dQuote{RMSE can be interpreted as the standard
#'         deviation of the unexplained variance, and has the useful property
#'         of being in the same units as the response variable. Lower values
#'         of RMSE indicate better fit. RMSE is a good measure of how accurately
#'         the model predicts the response, and is the most important criterion
#'         for fit if the main purpose of the model is prediction.}
#'         \cite{(Grace-Martin K: Assessing the Fit of Regression Models)}
#'         \cr \cr
#'         The normalized RMSE is the proportion of the RMSE related to the
#'         range of the response variable. Hence, lower values indicate
#'         less residual variance.
#'         }
#'         \item{\strong{Residual Standard Error}}{
#'         The residual standard error is the square root of the residual
#'         sum of squares divided by the residual degrees of freedom.
#'         }
#'         \item{\strong{Mean Square Error}}{
#'         The mean square error is the mean of the sum of squared residuals,
#'         i.e. it measures the average of the squares of the errors. Lower
#'         values (closer to zero) indicate better fit.
#'         }
#'       }
#'
#' @examples
#' data(efc)
#' fit <- lm(barthtot ~ c160age + c12hour, data = efc)
#' rmse(fit)
#' rse(fit)
#'
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' rmse(fit)
#' mse(fit)
#'
#' # normalized RMSE
#' library(nlme)
#' fit <- lme(distance ~ age, data = Orthodont)
#' rmse(fit, normalized = TRUE)
#'
#' @importFrom stats residuals df.residual
#' @export
rmse <- function(fit, normalized = FALSE) {
  # compute rmse
  rmse_val <- sqrt(mse(fit))

  # if normalized, divide by range of response
  if (normalized) {
    # get response
    resp <- resp_val(fit)
    # cpmpute rmse, normalized
    rmse_val <- rmse_val / (max(resp, na.rm = T) - min(resp, na.rm = T))
  }

  rmse_val
}


#' @rdname rmse
#' @name rse
#' @export
rse <- function(fit) {
  # Residual standard error
  sqrt(sum(stats::residuals(fit) ^ 2, na.rm = T) / stats::df.residual(fit))
}


#' @rdname rmse
#' @name mse
#' @export
mse <- function(fit) {
  # Mean square error
  mean(stats::residuals(fit) ^ 2, na.rm = T)
}
