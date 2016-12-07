#' @title Chi-square goodness-of-fit-test
#' @name chisq_gof
#'
#' @description This method performs a Chi-square goodness-of-fit-test (GOF)
#'                either on a numeric vector against probabilities, or
#'                a Goodness-of-fit test for \code{\link{glm}}-objects for binary data.
#'
#' @param x Numeric vector, or a \code{\link{glm}}-object.
#' @param prob Vector of probabilities (indicating the population probabilities) of the same length
#'          as \code{x}'s amount of categories / factor levels. Use \code{nrow(table(x))} to
#'          determine the amount of necessary values for \code{prob}. Only used,
#'          when \code{x} is a vector, and not a \code{glm}-object.
#' @param weights Vector with weights, used to weight \code{x}.
#' @return For vectors, returns the object of the computed \code{\link{chisq.test}}.
#'           \cr \cr
#'           For \code{glm}-objects, an object of class \code{chisq_gof} with
#'           following values:
#'           \itemize{
#'            \item \code{p.value}	the p-value for the goodness-of-fit test
#'            \item \code{z.score} the standardized z-score for the goodness-of-fit test
#'            \item \code{RSS} the residual sums of squares term
#'            \item \code{X2} the pearson chi-squared statistic
#'           }
#'
#' @note For vectors, this function is a convenient function for the \code{\link{chisq.test}},
#'         performing goodness-of-fit test.
#'         \cr \cr
#'         For \code{glm}-objects, this function performs a goodness-of-fit test
#'         based on the \code{X2GOFtest} function of the \CRANpkg{binomTools} package.
#'         A well-fitting model shows no significant difference between
#'         the model and the observed data, i.e. the reported p-values should be
#'         greater than 0.05.
#'
#' @examples
#' data(efc)
#' # differing from population
#' chisq_gof(efc$e42dep, c(0.3,0.2,0.22,0.28))
#' # equal to population
#' chisq_gof(efc$e42dep, prop.table(table(efc$e42dep)))
#'
#' # goodness-of-fit test for logistic regression
#' efc$services <- ifelse(efc$tot_sc_e > 0, 1, 0)
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep, data = efc,
#'            family = binomial(link = "logit"))
#' chisq_gof(fit)
#'
#' @importFrom stats na.omit fitted resid formula as.formula lm pnorm chisq.test
#' @export
chisq_gof <- function(x, prob = NULL, weights = NULL) {
  if (inherits(x, "glm")) {

    # This is an adapted version from the
    # "binomTools" package. The "X2GOFtest()"
    # function did not work when model data frame
    # had missing values.
    y_hat <- stats::fitted(x)
    wt <- x$prior.weight
    vJ <- wt * y_hat * (1 - y_hat)
    cJ <- (1 - 2 * y_hat) / vJ
    X2 <- sum(stats::resid(x, type = "pearson") ^ 2)
    form <- stats::as.formula(x$formula)
    form[[2]] <- as.name("cJ")

    # use model matrix instead of data values,
    # because data may contain more variables
    # than needed, and due to missing may have
    # different row length
    dat <- stats::na.omit(x$model)
    dat$cJ <- cJ
    dat$vJ <- vJ
    RSS <- sum(stats::resid(stats::lm(form, data = dat, weights = vJ)) ^ 2)
    A <- 2 * (length(y_hat) - sum(1 / wt))
    z <- (X2 - x$df.residual) / sqrt(A + RSS)
    p.value <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    chi2gof <- list(p.value = p.value,
                    z.score = z,
                    RSS = RSS,
                    X2 = X2)
    class(chi2gof) <- "chi2gof"
  } else {
    # check if we have probs
    if (is.null(prob)) {
      warning("`prob` needs to be specified.", call. = F)
      return(invisible(NULL))
    }
    # performs a Chi-square goodnes-of-fit-test
    if (!is.null(weights)) x <- weight(x, weights)
    dummy <- as.vector(table(x))
    # goodness of fit-test. x is one-dimensional and
    # y not given
    chi2gof <- stats::chisq.test(dummy, p = prob)
  }
  return(chi2gof)
}


#' @title Hosmer-Lemeshow Goodness-of-fit-test
#' @name hoslem_gof
#'
#' @description This method performs a Hosmer-Lemeshow goodness-of-fit-test
#'                for generalized linear (mixed) models for binary data.
#'
#' @param x Fitted \code{\link{glm}} or \code{\link[lme4]{glmer}} model.
#' @param g Number of bins to divide the data. Default is 10.
#'
#' @return An object of class \code{hoslem_test} with
#'           following values:
#'           \itemize{
#'            \item \code{chisq} the Hosmer-Lemeshow chi-squared statistic
#'            \item \code{df} degrees of freedom
#'            \item \code{p.value} the p-value for the goodness-of-fit test
#'           }
#'
#' @note A well-fitting model shows no significant difference between
#'         the model and the observed data, i.e. the reported p-value should be
#'         greater than 0.05.
#'
#' @seealso \code{\link{r2}}
#'
#' @examples
#' data(efc)
#' # goodness-of-fit test for logistic regression
#' efc$services <- ifelse(efc$tot_sc_e > 0, 1, 0)
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep, data = efc,
#'            family = binomial(link = "logit"))
#' hoslem_gof(fit)
#'
#' @importFrom stats fitted pchisq quantile xtabs
#' @export
hoslem_gof <- function(x, g = 10) {
  # check for valid object class
  if (!inherits(x, c("glmerMod", "glm"))) {
    stop("'x' must be an object of class 'glm' or 'glmerMod'.", call. = F)
  }

  # mixed models (lme4)
  if (inherits(x, "glmerMod")) {
    # check for package availability
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("Package 'lme4' needed for this function to work. Please install it.", call. = FALSE)
    }
    y <- lme4::getME(x, "y")
    yhat <- stats::fitted(x)
  } else {
    y <- x$y
    yhat <- stats::fitted(x)
  }
  cutyhat <- cut(yhat,
                 breaks = stats::quantile(yhat, probs = seq(0, 1, 1 / g)),
                 include.lowest = TRUE)
  obs <- stats::xtabs(cbind(1 - y, y) ~ cutyhat)
  expect <- stats::xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq <- sum((obs - expect)^2 / expect)
  p.value <- 1 - stats::pchisq(chisq, g - 2)
  hoslem <- list(chisq = chisq,
                 df = g - 2,
                 p.value = p.value)
  class(hoslem) <- "hoslem_test"
  return(hoslem)
}
