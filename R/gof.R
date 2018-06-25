#' @rdname rmse
#' @importFrom stats na.omit fitted resid formula as.formula lm pnorm chisq.test
#' @export
chisq_gof <- function(fit, prob = NULL, weights = NULL) {
  if (inherits(fit, "glm")) {

    # This is an adapted version from the
    # "binomTools" package. The "X2GOFtest()"
    # function did not work when model data frame
    # had missing values.
    y_hat <- stats::fitted(fit)
    wt <- fit$prior.weight
    vJ <- wt * y_hat * (1 - y_hat)
    cJ <- (1 - 2 * y_hat) / vJ
    X2 <- sum(stats::resid(fit, type = "pearson") ^ 2)
    form <- stats::as.formula(fit$formula)
    form[[2]] <- as.name("cJ")

    # use model matrix instead of data values,
    # because data may contain more variables
    # than needed, and due to missing may have
    # different row length
    dat <- stats::na.omit(fit$model)
    dat$cJ <- cJ
    dat$vJ <- vJ

    RSS <- sum(stats::resid(stats::lm(form, data = dat, weights = vJ)) ^ 2)
    A <- 2 * (length(y_hat) - sum(1 / wt))
    z <- (X2 - fit$df.residual) / sqrt(A + RSS)

    p.value <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

    chi2gof <- list(
      p.value = p.value,
      z.score = z,
      RSS = RSS,
      X2 = X2
    )
    class(chi2gof) <- "chi2gof"
  } else {
    # check if we have probs
    if (is.null(prob)) {
      warning("`prob` needs to be specified.", call. = F)
      return(invisible(NULL))
    }

    # performs a Chi-square goodnes-of-fit-test
    if (!is.null(weights)) fit <- weight(fit, weights)
    dummy <- as.vector(table(fit))

    # goodness of fit-test. x is one-dimensional and
    # y not given
    chi2gof <- stats::chisq.test(dummy, p = prob)
  }

  chi2gof
}


#' @rdname rmse
#' @importFrom stats fitted pchisq quantile xtabs
#' @export
hoslem_gof <- function(fit, n.bins = 10) {
  # check for valid object class
  if (!inherits(fit, c("glmerMod", "glm"))) {
    stop("'x' must be an object of class 'glm' or 'glmerMod'.", call. = F)
  }

  # mixed models (lme4)
  if (inherits(fit, "glmerMod")) {
    y <- lme4::getME(fit, "y")
    yhat <- stats::fitted(fit)
  } else {
    y <- fit$y
    yhat <- stats::fitted(fit)
  }

  cutyhat <- cut(
    yhat,
    breaks = stats::quantile(yhat, probs = seq(0, 1, 1 / n.bins)),
    include.lowest = TRUE
  )

  obs <- stats::xtabs(cbind(1 - y, y) ~ cutyhat)
  expect <- stats::xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq <- sum((obs - expect)^2 / expect)
  p.value <- 1 - stats::pchisq(chisq, n.bins - 2)

  hoslem <- list(
    chisq = chisq,
    df = n.bins - 2,
    p.value = p.value
  )

  class(hoslem) <- c("hoslem_test", "list")
  hoslem
}
