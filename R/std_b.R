#' @title Standardized beta coefficients and CI of linear and mixed models
#' @name std_beta
#' @description Returns the standardized beta coefficients, std. error and confidence intervals
#'                of a fitted linear (mixed) models.
#'
#' @param fit Fitted linear (mixed) model of class \code{lm} or
#'          \code{\link[lme4]{merMod}} (\CRANpkg{lme4} package).
#' @param type If \code{fit} is of class \code{lm}, normal standardized coefficients
#'          are computed by default. Use \code{type = "std2"} to follow
#'          \href{http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf}{Gelman's (2008)}
#'          suggestion, rescaling the estimates by deviding them by two standard
#'          deviations, so resulting coefficients are directly comparable for
#'          untransformed binary predictors.
#' @param ci.lvl Numeric, the level of the confidence intervals.
#'
#' @return A \code{tibble} with term names, standardized beta coefficients,
#'           standard error and confidence intervals of \code{fit}.
#'
#' @details \dQuote{Standardized coefficients refer to how many standard deviations a dependent variable will change,
#'         per standard deviation increase in the predictor variable. Standardization of the coefficient is
#'         usually done to answer the question of which of the independent variables have a greater effect
#'         on the dependent variable in a multiple regression analysis, when the variables are measured
#'         in different units of measurement (for example, income measured in dollars and family size
#'         measured in number of individuals).} \cite{(Source: Wikipedia)}
#'
#' @note For \code{\link[nlme]{gls}}-objects, standardized beta coefficients may be wrong
#'         for categorical variables (\code{factors}), because the \code{model.matrix} for
#'         \code{gls} objects returns the original data of the categorical vector,
#'         and not the 'dummy' coded vectors as for other classes. See, as example: \cr \cr
#'         \code{head(model.matrix(lm(neg_c_7 ~ as.factor(e42dep), data = efc, na.action = na.omit)))}
#'         \cr \cr and \cr \cr
#'         \code{head(model.matrix(nlme::gls(neg_c_7 ~ as.factor(e42dep), data = efc, na.action = na.omit)))}.
#'         \cr \cr
#'         In such cases, use \code{\link{to_dummy}} to create dummies from
#'         factors.
#'
#' @references \href{http://en.wikipedia.org/wiki/Standardized_coefficient}{Wikipedia: Standardized coefficient}
#'             \cr \cr
#'             Gelman A. 2008. Scaling regression inputs by dividing by two standard deviations. \emph{Statistics in Medicine 27: 2865–2873.} \url{http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf}
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
#' @importFrom stats model.matrix coef terms qnorm sd
#' @importFrom nlme getResponse
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map_if
#' @export
std_beta <- function(fit, type = "std", ci.lvl = .95) {
  # compute ci, two-ways
  ci <- 1 - ((1 - ci.lvl) / 2)

  # if we have merMod object (lme4), we need
  # other function to compute std. beta
  if (inherits(fit, c("lmerMod", "merModLmerTest")))
    return(sjs.stdmm(fit, ci.lvl))

  # has model intercept?
  tmp_i <- attr(stats::terms(fit), "intercept")
  has_intercept <- !is.null(tmp_i) & tmp_i == 1

  if (type == "std2" ) {
    # is package available?
    if (!requireNamespace("arm", quietly = TRUE)) {
      stop("Package `arm` needed for computing this type of standardized estimates. Please install it.", call. = FALSE)
    }
    # get standardized estimates
    b <- stats::coef(arm::standardize(fit))
  } else {
    # get coefficients
    b <- stats::coef(fit)
  }

  # remove intercept?
  if (has_intercept) b <- b[-1]

  if (type == "std2") {
    # stand. estimates need to be in variabel "beta"
    beta <- b
    # get standardized se
    beta.se <- summary(arm::standardize(fit))$coefficients[, 2]
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
      tibble::as_tibble()

    # get standard deviations for predictors
    sx <- purrr::map_dbl(fit.data, ~ stats::sd(.x, na.rm = TRUE))
    sy <- stats::sd(resp_val(fit), na.rm = TRUE)

    beta <- b * sx / sy

    if (inherits(fit, "gls"))
      se <- summary(fit)$tTable[, 2]
    else
      se <- summary(fit)$coef[, 2]

    # remove intercept?
    if (has_intercept) se <- se[-1]

    # compute standard error
    beta.se <- se * sx / sy
  }

  # return result
  tibble::tibble(
    term = names(b),
    std.estimate = beta,
    std.error = beta.se,
    conf.low = beta - stats::qnorm(ci) * beta.se,
    conf.high = beta + stats::qnorm(ci) * beta.se
  )
}


#' @importFrom stats sd coef
#' @importFrom lme4 fixef getME
#' @importFrom tibble tibble
sjs.stdmm <- function(fit, ci.lvl) {
  # compute ci, two-ways
  ci <- 1 - ((1 - ci.lvl) / 2)

  # code from Ben Bolker, see
  # http://stackoverflow.com/a/26206119/2094622
  sdy <- stats::sd(lme4::getME(fit, "y"))
  sdx <- apply(lme4::getME(fit, "X"), 2, sd)
  sc <- lme4::fixef(fit) * sdx / sdy
  se.fixef <- stats::coef(summary(fit))[, "Std. Error"]
  se <- se.fixef * sdx / sdy

  tibble::tibble(
    term = names(lme4::fixef(fit)),
    std.estimate = sc,
    std.error = se,
    conf.low = sc - stats::qnorm(ci) * se,
    conf.high = sc + stats::qnorm(ci) * se
  )
}
