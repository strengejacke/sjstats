#' @title Check model assumptions
#' @name outliers
#' @description \itemize{
#'                \item \code{outliers()} detects outliers in (generalized) linear models.
#'                \item \code{heteroskedastic()} checks a linear model for (non-)constant error variance.
#'                \item \code{autocorrelation()} checks for independence of errors.
#'                \item \code{normality()} checks linear models for (non-)normality of residuals.
#'                \item \code{multicollin()} checks predictors of linear models for multicollinearity.
#'                \item \code{check_assumptions()} checks all of the above assumptions.
#'              }
#'
#' @param x Fitted \code{lm}; for \code{outliers()}, may also be a \code{glm} model.
#' @param iterations Numeric, indicates the number of iterations to remove
#'        outliers.
#' @param as.logical Logical, if \code{TRUE}, the values returned by
#'        \code{check_assumptions()} are \code{TRUE} or \code{FALSE}, indicating
#'        whether each violation of model assumotion holds true or not. If
#'        \code{FALSE} (the default), the p-value of the respective test-statistics
#'        is returned.
#'
#' @return A tibble with the respective statistics.
#'
#' @details These functions are wrappers for respective functions that compute
#'          the various test statistics, however, each of them returns a tibble
#'          instead of a list of values.
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
#'          are preferable (see \code{\link[sjPlot]{sjp.lm}} with \code{type = "ma"}).
#'          \cr \cr
#'          \code{multicollin()} wraps \code{\link[car]{vif}} and returns
#'          the logical result as tibble. If \code{TRUE}, multicollinearity
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
#'       assumptions are alerted, thought the model is actually ok. It is
#'       preferable to check model assumptions based on visual inspection
#'       (see \code{\link[sjPlot]{sjp.lm}} with \code{type = "ma"}).
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
  # remember how many cases have been removed
  removedcases <- 0
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
          # count removed cases
          removedcases <- removedcases + length(vars)
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
          # count removed cases
          removedcases <- removedcases + length(vars)
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

  # do we have any outliers?
  if (removedcases < 1) {
    message("No outliers detected.")
  } else {
    message(sprintf("%i outliers removed in updated model.", removedcases))

    structure(class = "sjstats_outliers",
              list(result = rv,
                   updated.model = model,
                   removed.count = removedcases,
                   removed.obs = vars,
                   iterations = iterations - (maxcnt + 1))
              )
  }
}


#' @rdname outliers
#' @export
heteroskedastic <- function(x) {
  # check package availability
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package `car` needed for this function to work. Please install it.", call. = F)
  }

  # compute non-constant error variance test
  ts <- car::ncvTest(x)
  p.val <- ts$p

  if (p.val < 0.05) {
    message(sprintf("Heteroscedasticity (non-constant error variance) detected: p = %.3f", p.val))
  } else {
    message("Homoscedasticity (constant error variance) detected.")
  }

  tibble::tibble(heteroskedastic = p.val)
}


#' @rdname outliers
#' @export
autocorrelation <- function(x) {
  # check package availability
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package `car` needed for this function to work. Please install it.", call. = F)
  }

  # check for autocorrelation
  ts <- car::durbinWatsonTest(x)
  p.val <- ts$p

  if (p.val < 0.05) {
    message(sprintf("Autocorrelated residuals detected: p = %.3f", p.val))
  } else {
    message("No autocorrelated residuals detected.")
  }

  tibble::tibble(autocorrelation = p.val)
}


#' @rdname outliers
#' @importFrom stats shapiro.test rstandard
#' @export
normality <- function(x) {
  # check for normality of residuals
  ts <- stats::shapiro.test(stats::rstandard(x))
  p.val <- ts$p.value

  if (p.val < 0.05) {
    message(sprintf("Non-normality of residuals detected: p = %.3f", p.val))
  } else {
    message("Residuals are normally distributed.")
  }

  tibble::tibble(non.normality = p.val)
}


#' @rdname outliers
#' @importFrom car vif
#' @export
multicollin <- function(x) {
  # check package availability
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package `car` needed for this function to work. Please install it.", call. = F)
  }

  # check for autocorrelation
  ts <- sqrt(car::vif(x)) > 2

  if (any(ts)) {
    mp <- paste(sprintf("%s", names(ts)), collapse = ", ")
    message(paste0("Multicollinearity detected for following predictors: ", mp))
  } else {
    message("No multicollinearity detected.")
  }

  tibble::tibble(multicollin = any(ts))
}


#' @export
#' @rdname outliers
check_assumptions <- function(x, as.logical = FALSE) {
  # first, check outliers
  models <- suppressMessages(outliers(x))

  if (is.null(models))
    result <- tibble::tibble()
  else
    result <- models$result

  # check assumptions for original
  hn <- suppressMessages(heteroskedastic(x)$heteroskedastic)
  mn <- suppressMessages(multicollin(x)$multicollin)
  nn <- suppressMessages(normality(x)$non.normality)
  an <- suppressMessages(autocorrelation(x)$autocorrelation)

  # check whether user wants p-values or logical values
  if (as.logical) {
    hn <- hn < 0.05
    nn <- nn < 0.05
    an <- an < 0.05
  }

  # if we have an updated model w/o outliers, check assumptions for
  # this model as well
  if (!is.null(models)) {
    hu <- suppressMessages(heteroskedastic(models$updated.model)$heteroskedastic)
    mu <- suppressMessages(multicollin(models$updated.model)$multicollin)
    nu <- suppressMessages(normality(models$updated.model)$non.normality)
    au <- suppressMessages(autocorrelation(models$updated.model)$autocorrelation)

    # check whether user wants p-values or logical values
    if (as.logical) {
      hu <- hu < 0.05
      nu <- nu < 0.05
      au <- au < 0.05
    }

    result <- tibble::add_column(
      result,
      heteroskedasticity = c(hn, hu),
      multicollinearity = c(mn, mu),
      non.normal.resid = c(nn, nu),
      autocorrelation = c(an, au)
    )
  } else {
    result <- tibble::tibble(
      heteroskedasticity = hn,
      multicollinearity = mn,
      non.normal.resid = nn,
      autocorrelation = an
    )
  }

  result
}