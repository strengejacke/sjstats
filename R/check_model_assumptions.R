#' @title Check model assumptions
#' @name check_assumptions
#' @description \itemize{
#'                \item \code{outliers()} detects outliers in (generalized) linear models.
#'                \item \code{heteroskedastic()} checks a linear model for (non-)constant error variance.
#'                \item \code{autocorrelation()} checks for independence of errors.
#'                \item \code{normality()} checks linear models for (non-)normality of residuals.
#'                \item \code{multicollin()} checks predictors of linear models for multicollinearity.
#'                \item \code{check_assumptions()} checks all of the above assumptions.
#'              }
#'
#' @param x Fitted \code{lm} (for \code{outliers()}, may also be a \code{glm} model),
#'        or a (nested) data frame with a list-variable that contains fitted model
#'        objects.
#' @param model.column Name or index of the list-variable that contains the fitted
#'        model objects. Only applies, if \code{x} is a nested data frame (e.g
#'        with models fitted to \code{\link{bootstrap}} replicates).
#' @param iterations Numeric, indicates the number of iterations to remove
#'        outliers.
#' @param as.logical Logical, if \code{TRUE}, the values returned by
#'        \code{check_assumptions()} are \code{TRUE} or \code{FALSE}, indicating
#'        whether each violation of model assumotion holds true or not. If
#'        \code{FALSE} (the default), the p-value of the respective test-statistics
#'        is returned.
#' @param ... Other arguments, passed down to \code{\link[car]{durbinWatsonTest}}.
#'
#' @return A tibble with the respective statistics.
#'
#' @details These functions are wrappers that compute various test statistics,
#'          however, each of them returns a tibble instead of a list of values.
#'          Furthermore, all functions can also be applied to multiples models
#'          in stored in \emph{list-variables} (see 'Examples').
#'          \cr \cr
#'          \code{outliers()} wraps \code{\link[car]{outlierTest}} and iteratively
#'          removes outliers for \code{iterations} times, or if the r-squared value
#'          (for glm: the AIC) did not improve after removing outliers. The function
#'          returns a tibble with r-squared and AIC statistics for the original
#'          and updated model, as well as the update model itself (\code{$updated.model}),
#'          the number (\code{$removed.count}) and indices of the removed observations
#'          (\code{$removed.obs}).
#'          \cr \cr
#'          \code{heteroskedastic()} wraps \code{\link[car]{ncvTest}} and returns
#'          the p-value of the test statistics as tibble. A p-value < 0.05 indicates
#'          a non-constant variance (heteroskedasticity).
#'          \cr \cr
#'          \code{autocorrelation()} wraps \code{\link[car]{durbinWatsonTest}}
#'          and returns the p-value of the test statistics as tibble. A p-value
#'          < 0.05 indicates autocorrelated residuals. In such cases, robust
#'          standard errors (see \code{\link{robust}} return more accurate results
#'          for the estimates, or maybe a mixed model with error term for the
#'          cluster groups should be used.
#'          \cr \cr
#'          \code{normality()} calls \code{\link[stats]{shapiro.test}}
#'          and checks the standardized residuals for normal distribution.
#'          The p-value of the test statistics is returned as tibble. A p-value
#'          < 0.05 indicates a significant deviation from normal distribution.
#'          Note that this formal test almost always yields significant results
#'          for the distribution of residuals and visual inspection (e.g. qqplots)
#'          are preferable (see \code{\link[sjPlot]{plot_model}} with
#'          \code{type = "diag"}).
#'          \cr \cr
#'          \code{multicollin()} wraps \code{\link[car]{vif}} and returns
#'          the logical result as tibble. \code{TRUE}, if multicollinearity
#'          exists, else not. In case of multicollinearity, the names of independent
#'          variables that vioalte contribute to multicollinearity are printed
#'          to the console.
#'          \cr \cr
#'          \code{check_assumptions()} runs all of the above tests and returns
#'          a tibble with all test statistics included. In case the p-values
#'          are too confusing, use the \code{as.logical} argument, where all
#'          p-values are replaced with either \code{TRUE} (in case of violation)
#'          or \code{FALSE} (in case of model conforms to assumption of linar
#'          regression).
#'
#' @note These formal tests are very strict and in most cases violation of model
#'       assumptions are alerted, though the model is actually ok. It is
#'       preferable to check model assumptions based on visual inspection
#'       (see \code{\link[sjPlot]{plot_model}} with \code{type = "diag"}).
#'
#' @examples
#' data(efc)
#'
#' fit <- lm(barthtot ~ c160age + c12hour + c161sex + c172code, data = efc)
#' outliers(fit)
#' heteroskedastic(fit)
#' autocorrelation(fit)
#' normality(fit)
#' check_assumptions(fit)
#'
#' fit <- lm(barthtot ~ c160age + c12hour + c161sex + c172code + neg_c_7,
#'           data = efc)
#' outliers(fit)
#' check_assumptions(fit, as.logical = TRUE)
#'
#' # apply function to multiple models in list-variable
#' library(purrr)
#' library(dplyr)
#' tmp <- efc %>%
#'   bootstrap(50) %>%
#'   mutate(
#'     models = map(strap, ~lm(neg_c_7 ~ e42dep + c12hour + c161sex, data = .x))
#'   )
#'
#' # for list-variables, argument 'model.column' is the
#' # quoted name of the list-variable with fitted models
#' tmp %>% normality("models")
#' tmp %>% heteroskedastic("models")
#'
#' # Durbin-Watson-Test from package 'car' takes a little bit longer due
#' # to simulation of p-values...
#' \dontrun{
#' tmp %>% check_assumptions("models", as.logical = TRUE, reps = 100)}
#'
#' @importFrom tibble tibble
#' @export
check_assumptions <- function(x, model.column = NULL, as.logical = FALSE, ...) {
  # check assumptions for original
  hn <- suppressMessages(heteroskedastic(x, model.column)$heteroskedastic)
  mn <- suppressMessages(multicollin(x, model.column)$multicollin)
  nn <- suppressMessages(normality(x, model.column)$non.normality)
  an <- suppressMessages(autocorrelation(x, model.column, ...)$autocorrelation)

  # check whether user wants p-values or logical values
  if (as.logical) {
    hn <- hn < 0.05
    nn <- nn < 0.05
    an <- an < 0.05
  }

  tibble::tibble(
    heteroskedasticity = hn,
    multicollinearity = mn,
    non.normal.resid = nn,
    autocorrelation = an
  )
}


#' @rdname check_assumptions
#' @importFrom sjmisc is_empty
#' @importFrom stats update
#' @importFrom tibble tibble
#' @export
outliers <- function(x, iterations = 5) {
  # check package availability
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package `car` needed for this function to work. Please install it.", call. = F)
  }

  # copy current model
  model <- x
  # maximum loops
  maxcnt <- iterations
  outlier <- c()
  loop <- TRUE

  # check outlier for glm-models
  if (inherits(x, "glm")) {
    # get AIC-Value, for reference when to stop removing outliers.
    aic <- model$aic
    # start loop
    while (isTRUE(loop)) {
      # get outliers of model
      vars <- as.numeric(names(which(car::outlierTest(model, cutoff = Inf, n.max = Inf)$bonf.p < 1)))

      # no outliers found? then stop...
      if (sjmisc::is_empty(vars)) {
        loop <- FALSE
      } else {
        # remove outliers
        dummymodel <- stats::update(model, subset = -vars)
        # retrieve new AIC-value
        dummyaic <- dummymodel$aic
        # decrease maximum loops
        maxcnt <- maxcnt - 1

        # check whether AIC-value of updated model is larger
        # than previous AIC-value or if we have already all iterations done,
        if (dummyaic >= aic || maxcnt < 1) {
          loop <- FALSE
        } else {
          # else copy new model, which is the better one (according to AIC-value)
          model <- dummymodel
          # and get new AIC-value
          aic <- dummyaic
          # add outliers to final return value
          outlier <- c(outlier, vars)
        }
      }
    }

    # create return tibble
    rv <- tibble::tibble(
      models = c("original", "updated"),
      aic = c(stats::AIC(x), stats::AIC(model))
    )

  } else if (inherits(model, "lm")) {
    # get r2
    rs <- summary(model)$r.squared
    # start loop
    while (isTRUE(loop)) {
      # get outliers of model
      vars <- as.numeric(names(which(car::outlierTest(model, cutoff = Inf, n.max = Inf)$bonf.p < 1)))

      # do we have any outliers?
      if (sjmisc::is_empty(vars)) {
        loop <- FALSE
      } else {
        # retrieve variable numbers of outliers
        # update model by removing outliers
        dummymodel <- stats::update(model, subset = -vars)
        # retrieve new r2
        dummyrs <- summary(dummymodel)$r.squared
        # decrease maximum loops
        maxcnt <- maxcnt - 1

        # check whether r2 of updated model is lower
        # than previous r2 or if we have already all loop-steps done,
        # stop loop
        if (dummyrs < rs || maxcnt < 1) {
          loop <- FALSE
        } else {
          # else copy new model, which is the better one (according to r2)
          model <- dummymodel
          # and get new r2
          rs <- dummyrs
          # add outliers to final return value
          outlier <- c(outlier, vars)
        }
      }
    }

    # create return tibble
    rv <- tibble::tibble(
      models = c("original", "updated"),
      adjusted.r2 = c(summary(x)$adj.r.squared,
                      summary(model)$adj.r.squared),
      aic = c(stats::AIC(x), stats::AIC(model))
    )

  } else {
    stop("Model is not of class `lm` or `glm`.", call. = F)
  }

  # count removed cases
  removedcases <- length(outlier)

  # do we have any outliers?
  if (removedcases < 1) {
    message("No outliers detected.")
  } else {
    message(sprintf("%i outliers removed in updated model.", removedcases))

    structure(
      class = "sjstats_outliers",
      list(
        result = rv,
        updated.model = model,
        removed.count = removedcases,
        removed.obs = outlier,
        iterations = iterations - (maxcnt + 1)
      )
    )
  }
}


#' @rdname check_assumptions
#' @importFrom purrr map_dbl map
#' @export
heteroskedastic <- function(x, model.column = NULL) {
  # check package availability
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package `car` needed for this function to work. Please install it.", call. = F)
  }

  # check if we have list-variable, e.g. from nested data frames
  if (!is.null(model.column) && inherits(x[[model.column]], "list")) {

    p.val <- x[[model.column]] %>%
      # iterate all model columns in nested data frame
      purrr::map(~ .x) %>%
      # call ncvTest for each model, and just get p-value
      purrr::map_dbl(~ nonconstvar(.x))
  } else {
    # compute non-constant error variance test
    p.val <- nonconstvar(x)

    # print message, but not for nested models. only if 'x' is a single model
    if (p.val < 0.05) {
      message(sprintf("Heteroscedasticity (non-constant error variance) detected: p = %.3f", p.val))
    } else {
      message("Homoscedasticity (constant error variance) detected.")
    }
  }

  tibble::tibble(heteroskedastic = p.val)
}


#' @rdname check_assumptions
#' @export
autocorrelation <- function(x, model.column = NULL, ...) {
  # check package availability
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package `car` needed for this function to work. Please install it.", call. = F)
  }

  # check if we have list-variable, e.g. from nested data frames
  if (!is.null(model.column) && inherits(x[[model.column]], "list")) {

    p.val <-
      # iterate all model columns in nested data frame
      purrr::map(x[[model.column]], ~ .x) %>%
      # call ncvTest for each model, and just get p-value
      purrr::map_dbl(~ car::durbinWatsonTest(.x)$p)
  } else {
    # check for autocorrelation
    ts <- car::durbinWatsonTest(x, ...)
    p.val <- ts$p

    if (p.val < 0.05) {
      message(sprintf("Autocorrelated residuals detected: p = %.3f", p.val))
    } else {
      message("No autocorrelated residuals detected.")
    }
  }

  tibble::tibble(autocorrelation = p.val)
}


#' @rdname check_assumptions
#' @importFrom stats shapiro.test rstandard
#' @export
normality <- function(x, model.column = NULL) {
  # check if we have list-variable, e.g. from nested data frames
  if (!is.null(model.column) && inherits(x[[model.column]], "list")) {

    p.val <-
      # iterate all model columns in nested data frame
      purrr::map(x[[model.column]], ~ .x) %>%
      # call ncvTest for each model, and just get p-value
      purrr::map_dbl(~ stats::shapiro.test(stats::rstandard(.x))$p.value)
  } else {
    # check for normality of residuals
    ts <- stats::shapiro.test(stats::rstandard(x))
    p.val <- ts$p.value

    if (p.val < 0.05) {
      message(sprintf("Non-normality of residuals detected: p = %.3f", p.val))
    } else {
      message("Residuals are normally distributed.")
    }
  }

  tibble::tibble(non.normality = p.val)
}


#' @rdname check_assumptions
#' @importFrom purrr map_lgl
#' @export
multicollin <- function(x, model.column = NULL) {
  # check package availability
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package `car` needed for this function to work. Please install it.", call. = F)
  }

  # check if we have list-variable, e.g. from nested data frames
  if (!is.null(model.column) && inherits(x[[model.column]], "list")) {
    ts <-
      # iterate all model columns in nested data frame
      purrr::map(x[[model.column]], ~ .x) %>%
      # call vif for each model, and just get p-value
      purrr::map_lgl(~ any(sqrt(car::vif(.x)) > 2))
  } else {
    # check for autocorrelation
    ts <- sqrt(car::vif(x)) > 2

    if (any(ts)) {
      mp <- paste(sprintf("%s", names(ts)), collapse = ", ")
      message(paste0("Multicollinearity detected for following predictors: ", mp))
    } else {
      message("No multicollinearity detected.")
    }

    ts <- any(ts)
  }

  tibble::tibble(multicollin = ts)
}


#' @importFrom stats residuals df.residual pchisq anova
nonconstvar <- function(model) {
  sumry <- summary(model)

  residuals <- stats::residuals(model, type = "pearson")
  S.sq <- stats::df.residual(model) * (sumry$sigma) ^ 2 / sum(!is.na(residuals))

  .U <- (residuals ^ 2) / S.sq
  mod <- lm(.U ~ fitted.values(model))

  SS <- stats::anova(mod)$"Sum Sq"
  RegSS <- sum(SS) - SS[length(SS)]
  Chisq <- RegSS / 2

  stats::pchisq(Chisq, df = 1, lower.tail = FALSE)
}
