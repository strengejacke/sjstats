#' @title Standardized beta coefficients and CI of linear and mixed models
#' @name std_beta
#' @description Returns the standardized beta coefficients, std. error and confidence intervals
#'                of a fitted linear (mixed) models.
#'
#' @param fit Fitted linear (mixed) model of class \code{lm}, \code{merMod}
#'   (\CRANpkg{lme4} package), \code{gls} or \code{stanreg}.
#' @param type If \code{fit} is of class \code{lm}, normal standardized coefficients
#'          are computed by default. Use \code{type = "std2"} to follow
#'          \href{http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf}{Gelman's (2008)}
#'          suggestion, rescaling the estimates by deviding them by two standard
#'          deviations, so resulting coefficients are directly comparable for
#'          untransformed binary predictors.
#' @param ci.lvl Numeric, the level of the confidence intervals.
#' @param ... Currently not used.
#'
#' @return A \code{tibble} with term names, standardized beta coefficients,
#'           standard error and confidence intervals of \code{fit}.
#'
#' @details \dQuote{Standardized coefficients refer to how many standard deviations a dependent variable will change,
#'         per standard deviation increase in the predictor variable. Standardization of the coefficient is
#'         usually done to answer the question of which of the independent variables have a greater effect
#'         on the dependent variable in a multiple regression analysis, when the variables are measured
#'         in different units of measurement (for example, income measured in dollars and family size
#'         measured in number of individuals)} \cite{(Source: Wikipedia)}
#'
#' @note For \code{\link[nlme]{gls}}-objects, standardized beta coefficients may be wrong
#'         for categorical variables (\code{factors}), because the \code{model.matrix} for
#'         \code{gls} objects returns the original data of the categorical vector,
#'         and not the 'dummy' coded vectors as for other classes. See, as example: \cr \cr
#'         \code{head(model.matrix(lm(neg_c_7 ~ as.factor(e42dep), data = efc, na.action = na.omit)))}
#'         \cr \cr and \cr \cr
#'         \code{head(model.matrix(nlme::gls(neg_c_7 ~ as.factor(e42dep), data = efc, na.action = na.omit)))}.
#'         \cr \cr
#'         In such cases, use \code{\link[sjmisc]{to_dummy}} to create dummies from
#'         factors.
#'
#' @references \href{http://en.wikipedia.org/wiki/Standardized_coefficient}{Wikipedia: Standardized coefficient}
#'             \cr \cr
#'             Gelman A. 2008. Scaling regression inputs by dividing by two standard deviations. Statistics in Medicine 27: 2865-2873 \url{http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf}
#'
#' @examples
#' # fit linear model
#' fit <- lm(Ozone ~ Wind + Temp + Solar.R, data = airquality)
#' # print std. beta coefficients
#' std_beta(fit)
#'
#' # print std. beta coefficients and ci, using
#' # 2 sd and center binary predictors
#' std_beta(fit, type = "std2")
#'
#' # std. beta for mixed models
#' library(lme4)
#' fit1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' std_beta(fit)
#'
#' @export
std_beta <- function(fit, ...) {
  UseMethod("std_beta")
}


#' @importFrom stats sd coef qnorm
#' @importFrom lme4 fixef getME
#' @importFrom dplyr slice
#' @rdname std_beta
#' @export
std_beta.merMod <- function(fit, ci.lvl = .95, ...) {
  # compute ci, two-ways
  ci <- 1 - ((1 - ci.lvl) / 2)

  # code from Ben Bolker, see
  # http://stackoverflow.com/a/26206119/2094622
  sdy <- stats::sd(lme4::getME(fit, "y"))
  sdx <- apply(lme4::getME(fit, "X"), 2, sd)
  sc <- lme4::fixef(fit) * sdx / sdy
  se.fixef <- stats::coef(summary(fit))[, "Std. Error"]
  se <- se.fixef * sdx / sdy

  data_frame(
    term = names(lme4::fixef(fit)),
    std.estimate = sc,
    std.error = se,
    conf.low = sc - stats::qnorm(ci) * se,
    conf.high = sc + stats::qnorm(ci) * se
  ) %>%
    dplyr::slice(-1)

}


#' @importFrom stats model.matrix coef terms qnorm sd
#' @importFrom sjlabelled as_numeric
#' @importFrom purrr map_if map_dbl
#' @rdname std_beta
#' @export
std_beta.lm <- function(fit, type = "std", ci.lvl = .95, ...) {
  std_beta_helper(fit = fit, type = type, ci.lvl = ci.lvl, se = summary(fit)$coef[, 2], ...)
}


#' @rdname std_beta
#' @export
std_beta.gls <- function(fit, type = "std", ci.lvl = .95, ...) {
  std_beta_helper(fit = fit, type = type, ci.lvl = ci.lvl, se = summary(fit)$tTable[, 2], ...)
}


std_beta_helper <- function(fit, type, ci.lvl, se, ...) {
  # compute ci, two-ways
  ci <- 1 - ((1 - ci.lvl) / 2)

  # has model intercept?
  tmp_i <- attr(stats::terms(fit), "intercept")
  has_intercept <- !is.null(tmp_i) & tmp_i == 1

  if (type == "std2") {
    fit <- std2(fit)
  }

  b <- stats::coef(fit)

  # remove intercept?
  if (has_intercept) b <- b[-1]

  if (type == "std2") {
    # stand. estimates need to be in variabel "beta"
    beta <- b
    # get standardized se
    beta.se <- summary(fit)$coefficients[, 2]
    # remove intercept?
    if (has_intercept) beta.se <- beta.se[-1]
  } else {
    # get data as data frame
    fit.data <- as.data.frame(stats::model.matrix(fit))
    # remove intercept?
    if (has_intercept) fit.data <- fit.data[, -1, drop = FALSE]

    # convert factor to numeric, else sd throws a warning
    fit.data <- fit.data %>%
      purrr::map_if(is.factor, ~sjlabelled::as_numeric(.x, keep.labels = F)) %>%
      as.data.frame()

    # get standard deviations for predictors
    sx <- purrr::map_dbl(fit.data, ~ stats::sd(.x, na.rm = TRUE))
    sy <- stats::sd(resp_val(fit), na.rm = TRUE)

    beta <- b * sx / sy

    # remove intercept?
    if (has_intercept) se <- se[-1]

    # compute standard error
    beta.se <- se * sx / sy
  }

  # return result
  data_frame(
    term = names(b),
    std.estimate = beta,
    std.error = beta.se,
    conf.low = beta - stats::qnorm(ci) * beta.se,
    conf.high = beta + stats::qnorm(ci) * beta.se
  )
}


#' @importFrom stats model.matrix sd
#' @export
std_beta.stanreg <- function(fit, ...) {
  if (!requireNamespace("rstanarm", quietly = TRUE))
    stop("Please install and load package `rstanarm` first.", call. = F)

  X <- stats::model.matrix(fit)
  sd_X <- apply(X, MARGIN = 2, FUN = stats::sd)[-1]
  sd_Y <- apply(rstanarm::posterior_predict(fit), MARGIN = 1, FUN = stats::sd)
  beta <- as.matrix(fit)[ , 2:ncol(X), drop = FALSE]
  sweep(
    sweep(beta, MARGIN = 2, STATS = sd_X, FUN = `*`),
    MARGIN = 1, STATS = sd_Y, FUN = `/`
  )
}


#' @importFrom dplyr n_distinct
#' @importFrom sjlabelled as_numeric
#' @importFrom sjmisc std
#' @importFrom stats weights
std2 <- function(x) {
  form <- stats::formula(x)
  data <- model_frame(x)
  terms <- pred_vars(x)
  resp <- resp_var(x)

  newdata <- purrr::map(colnames(data), function(.x) {
    v <- data[[.x]]

    if (.x %in% terms) {
      if (dplyr::n_distinct(v, na.rm = TRUE) == 2) {
        v <- sjlabelled::as_numeric(v)
        v <- v - mean(v, na.rm = T)
      } else if (is.numeric(v) && .x != "(weights)") {
        v <- sjmisc::std(v, robust = "2sd")
      }
    }

    v
  })

  newdata <- as.data.frame(newdata)
  colnames(newdata) <- colnames(data)

  w <- stats::weights(x)
  newdata <- newdata[, c(resp, terms)]
  if (!is.null(w)) {
    newdata <- cbind(newdata, w)
    lm(form, data = newdata, weights = w)
  } else {
    lm(form, data = newdata)
  }
}
