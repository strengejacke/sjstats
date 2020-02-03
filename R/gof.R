#' @title Compute model quality
#' @name chisq_gof
#'
#' @description For logistic regression models, performs a Chi-squared
#'   goodness-of-fit-test.
#'
#' @param x A numeric vector or a \code{glm}-object.
#' @param prob Vector of probabilities (indicating the population probabilities)
#'   of the same length as \code{x}'s amount of categories / factor levels.
#'   Use \code{nrow(table(x))} to determine the amount of necessary values
#'   for \code{prob}. Only used, when \code{x} is a vector, and not a
#'   \code{glm}-object.
#' @param weights Vector with weights, used to weight \code{x}.
#'
#' @references
#'   Hosmer, D. W., & Lemeshow, S. (2000). Applied Logistic Regression. Hoboken, NJ, USA: John Wiley & Sons, Inc. \doi{10.1002/0471722146}
#'
#' @details For vectors, this function is a convenient function for the
#'         \code{chisq.test()}, performing goodness-of-fit test. For
#'         \code{glm}-objects, this function performs a goodness-of-fit test.
#'         A well-fitting model shows \emph{no} significant difference between the
#'         model and the observed data, i.e. the reported p-values should be
#'         greater than 0.05.
#'
#' @return For vectors, returns the object of the computed \code{\link[stats]{chisq.test}}.
#'     For \code{glm}-objects, an object of class \code{chisq_gof} with
#'     following values: \code{p.value}, the p-value for the goodness-of-fit test;
#'     \code{z.score}, the standardized z-score for the goodness-of-fit test;
#'     \code{rss}, the residual sums of squares term and \code{chisq}, the pearson
#'     chi-squared statistic.
#'
#' @examples
#' data(efc)
#' efc$neg_c_7d <- ifelse(efc$neg_c_7 < median(efc$neg_c_7, na.rm = TRUE), 0, 1)
#' m <- glm(
#'   neg_c_7d ~ c161sex + barthtot + c172code,
#'   data = efc,
#'   family = binomial(link = "logit")
#' )
#'
#' # goodness-of-fit test for logistic regression
#' chisq_gof(m)
#'
#' # goodness-of-fit test for vectors against probabilities
#' # differing from population
#' chisq_gof(efc$e42dep, c(0.3,0.2,0.22,0.28))
#'
#' # equal to population
#' chisq_gof(efc$e42dep, prop.table(table(efc$e42dep)))
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
    X2 <- sum(stats::resid(x, type = "pearson")^2)
    form <- stats::as.formula(x$formula)
    form[[2]] <- as.name("cJ")

    # use model matrix instead of data values,
    # because data may contain more variables
    # than needed, and due to missing may have
    # different row length
    dat <- stats::na.omit(x$model)
    dat$cJ <- cJ
    dat$vJ <- vJ

    RSS <- sum(stats::resid(stats::lm(form, data = dat, weights = vJ))^2)
    A <- 2 * (length(y_hat) - sum(1 / wt))
    z <- (X2 - x$df.residual) / sqrt(A + RSS)

    p.value <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

    chi2gof <- list(
      p.value = p.value,
      z.score = z,
      rss = RSS,
      chisq = X2
    )
    class(chi2gof) <- c("sj_chi2gof", "list")
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

  chi2gof
}
