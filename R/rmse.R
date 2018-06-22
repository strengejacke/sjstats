#' @title Compute model quality
#' @name rmse
#' @description Compute root mean squared error, residual standard error or
#'              mean square error of fitted linear (mixed effects) models. Or
#'              error rate and binned residuals for logistic regression models.
#'
#' @param fit Fitted linear model of class \code{lm}, \code{merMod} (\pkg{lme4})
#'   or \code{lme} (\pkg{nlme}). For \code{error_rate()} and \code{binned_resid()},
#'   a \code{glm}-object with binomial-family.
#' @param normalized Logical, use \code{TRUE} if normalized rmse should be returned.
#' @param term Name of independent variable from \code{fit}. If not \code{NULL},
#'   average residuals for the categories of \code{term} are plotted; else,
#'   average residuals for the estimated probabilities of the response are
#'   plotted.
#' @param n.bins Numeric, the number of "bins".
#'
#' @seealso \code{\link{r2}} for R-squared or pseude-R-squared values, and
#'            \code{\link{cv}} for the coefficient of variation.
#'
#' @references
#'   Gelman A, Hill J (2007) Data Analysis Using Regression and Multilevel/Hierarchical Models. Cambridge, New York: Cambridge University Press
#'   \cr \cr
#'   \href{http://www.theanalysisfactor.com/assessing-the-fit-of-regression-models/}{Grace-Martin K: Assessing the Fit of Regression Models}
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
#'         \item{\strong{Error Rate}}{
#'         The error rate is a crude measure for model fit for logistic regression
#'         models. It is defined as the proportion of cases for which the
#'         deterministic prediction is wrong, i.e. the proportion where the the
#'         predicted probability is above 0.5, although y = 0 (and vice versa).
#'         In general, the error rate should be below 0.5 (i.e. 50\%), the
#'         closer to zero, the better. Furthermore, the error rate of the full
#'         model should be considerably below the null model's error rate
#'         (cf. Gelman and Hill 2007, pp. 99).
#'         }
#'         \item{\strong{Binned Residuals}}{
#'         (cf. Gelman and Hill 2007, pp. 97ff).
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
#' # Error Rate
#' efc$neg_c_7d <- ifelse(efc$neg_c_7 < median(efc$neg_c_7, na.rm = TRUE), 0, 1)
#' m <- glm(
#'   neg_c_7d ~ c161sex + barthtot + c172code,
#'   data = efc,
#'   family = binomial(link = "logit")
#' )
#' error_rate(m)
#'
#' # Binned residuals
#' binned_resid(m)
#' binned_resid(m, "barthtot")#'
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


#' @rdname rmse
#' @name error_rate
#' @export
error_rate <- function(fit) {
  m0 <- glm(
    formula = as.formula(sprintf("%s ~ 1", resp_var(fit))),
    family = stats::binomial(link = "logit"),
    data = model_frame(fit)
  )

  y1 <- resp_val(fit)
  y0 <- resp_val(m0)

  p1 <- stats::predict.glm(fit, type = "response")
  error1 <- mean((p1 > .5 & y1 == 0) | (p1 <= .5 & y1 == 1))

  p0 <- stats::predict.glm(m0, type = "response")
  error0 <- mean((p0 > .5 & y0 == 0) | (p0 <= .5 & y0 == 1))

  er <- list(error.model = error1, error.null = error0)
  class(er) <- c("sj_error_rate", class(er))

  er
}


#' @importFrom sjmisc recode_to to_value
#' @importFrom stats fitted
#' @importFrom purrr map_df
#' @importFrom sjlabelled get_label
#' @export
binned_resid <- function(fit, term = NULL, n.bins = NULL) {
  fv <- stats::fitted(fit)
  mf <- model_frame(fit)

  if (is.null(term))
    x <- fv
  else
    x <- mf[[term]]

  y <- sjmisc::recode_to(sjmisc::to_value(resp_val(fit))) - fv

  if (is.null(n.bins)) n.bins <- round(sqrt(length(x)))

  breaks.index <- floor(length(x) * (1:(n.bins - 1)) / n.bins)
  breaks <- unique(c(-Inf, sort(x)[breaks.index], Inf))

  x.binned <- as.numeric(cut(x, breaks))

  d <- suppressWarnings(
    purrr::map_df(1:n.bins, function(.x) {
      items <- (1:length(x))[x.binned == .x]
      x.range <- range(x[items], na.rm = TRUE)
      xbar <- mean(x[items], na.rm = TRUE)
      ybar <- mean(y[items], na.rm = TRUE)
      n <- length(items)
      sdev <- sd(y[items], na.rm = TRUE)
      data.frame(xbar = xbar, ybar = ybar, n = n, x.lo = x.range[1], x.hi = x.range[2], se = 2 * sdev / sqrt(n))
    })
  )

  d <- d[complete.cases(d), ]

  gr <- abs(d$ybar) > abs(d$se)
  d$group <- "yes"
  d$group[gr] <- "no"

  class(d) <- c("sj_binres", class(d))
  attr(d, "resp_var") <- resp_var(fit)
  attr(d, "term") <- term
  if (!is.null(term))
    attr(d, "term.label") <- sjlabelled::get_label(mf[[term]], def.value = term)
  else
    attr(d, "term.label") <- term

  d
}
