#' @title Compute model quality
#' @name rmse
#'
#' @description Compute various measures or tests to assess the model quality,
#'   like root mean squared error, residual standard error or mean square error
#'   of fitted linear (mixed effects) models. For logistic regression models,
#'   or mixed models with binary outcome, the error rate, binned residuals,
#'   Chi-square goodness-of-fit-test or the Hosmer-Lemeshow Goodness-of-fit-test
#'   can be performed.
#'
#' @param fit Fitted linear model of class \code{lm}, \code{merMod} (\pkg{lme4})
#'   or \code{lme} (\pkg{nlme}). For \code{error_rate()} and \code{binned_resid()},
#'   a \code{glm}-object with binomial-family. For \code{chisq_gof()}, a
#'   numeric vector, or a \code{glm}-object.
#' @param normalized Logical, use \code{TRUE} if normalized rmse should be returned.
#' @param term Name of independent variable from \code{fit}. If not \code{NULL},
#'   average residuals for the categories of \code{term} are plotted; else,
#'   average residuals for the estimated probabilities of the response are
#'   plotted.
#' @param n.bins Numeric, the number of bins to divide the data. For
#'   \code{hoslem_gof()}, the default is 10. For \code{binned_resid()}, if
#'   \code{n.bins = NULL}, the square root of the number of observations is
#'   taken.
#' @param prob Vector of probabilities (indicating the population probabilities)
#'   of the same length as \code{fit}'s amount of categories / factor levels.
#'   Use \code{nrow(table(fit))} to determine the amount of necessary values
#'   for \code{prob}. Only used, when \code{fit} is a vector, and not a
#'   \code{glm}-object.
#' @param weights Vector with weights, used to weight \code{fit}.
#'
#' @seealso \code{\link{r2}} for R-squared or pseude-R-squared values, and
#'            \code{\link{cv}} for the coefficient of variation.
#'
#' @references
#'   Gelman A, Hill J (2007) Data Analysis Using Regression and Multilevel/Hierarchical Models. Cambridge, New York: Cambridge University Press
#'   \cr \cr
#'   Hosmer, D. W., & Lemeshow, S. (2000). Applied Logistic Regression. Hoboken, NJ, USA: John Wiley & Sons, Inc. \doi{10.1002/0471722146}
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
#'         \item{\strong{Chi-squared Goodness-of-Fit Test}}{
#'         For vectors, this function is a convenient function for the
#'         \code{chisq.test()}, performing goodness-of-fit test.
#'         \cr \cr
#'         For \code{glm}-objects, this function performs a goodness-of-fit test
#'         based on the \code{X2GOFtest()} function of the \CRANpkg{binomTools}
#'         package. A well-fitting model shows no significant difference between
#'         the model and the observed data, i.e. the reported p-values should be
#'         greater than 0.05.
#'         }
#'         \item{\strong{Hosmer-Lemeshow Goodness-of-Fit Test}}{
#'         A well-fitting model shows no significant difference between
#'         the model and the observed data, i.e. the reported p-value should be
#'         greater than 0.05.
#'         }
#'       }
#'
#' @return \describe{
#'   \item{\code{chisq_gof()}}{
#'     For vectors, returns the object of the computed \code{\link[stats]{chisq.test}}.
#'     \cr \cr
#'     For \code{glm}-objects, an object of class \code{chisq_gof} with
#'     following values:
#'     \itemize{
#'       \item \code{p.value} the p-value for the goodness-of-fit test
#'       \item \code{z.score} the standardized z-score for the goodness-of-fit test
#'       \item \code{RSS} the residual sums of squares term
#'       \item \code{X2} the pearson chi-squared statistic
#'     }
#'   }
#'   \item{\code{hoslem_gof()}}{
#'     An object of class \code{hoslem_test} with
#'     following values:
#'     \itemize{
#'       \item \code{chisq} the Hosmer-Lemeshow chi-squared statistic
#'       \item \code{df} degrees of freedom
#'       \item \code{p.value} the p-value for the goodness-of-fit test
#'     }
#'   }
#' }
#'
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
#' # goodness-of-fit test for logistic regression
#' chisq_gof(m)
#'
#' # goodness-of-fit test for logistic regression
#' hoslem_gof(m)
#'
#' # goodness-of-fit test for vectors against probabilities
#' # differing from population
#' chisq_gof(efc$e42dep, c(0.3,0.2,0.22,0.28))
#' # equal to population
#' chisq_gof(efc$e42dep, prop.table(table(efc$e42dep)))#'
#'
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


#' @rdname rmse
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
