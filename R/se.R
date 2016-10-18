utils::globalVariables(c("strap", "models"))

#' @title Standard Error for variables or coefficients
#' @name se
#' @description Compute standard error for a variable, for all variables
#'                of a data frame, for joint random and fixed effects
#'                coefficients of mixed models, or for intraclass correlation
#'                coefficients (ICC).
#'
#' @param x (Numeric) vector, a data frame, a \code{merMod}-object
#'          as returned by the \code{\link[lme4]{lmer}}-method,
#'          or a list with estimate and p-value. For the latter case, the list
#'          must contain elements named \code{estimate} and \code{p.value}
#'          (see 'Examples' and 'Details').
#' @param nsim Numeric, the number of simulations for calculating the
#'          standard error for intraclass correlation coefficients, as
#'          obtained by the \code{\link{icc}}-function.
#'
#' @return The standard error of \code{x}, or for each variable
#'           if \code{x} is a data frame, or for the coefficients
#'           of a mixed model (see \code{\link[lme4]{coef.merMod}}).
#'
#' @note Computation of standard errors for coefficients of mixed models
#'         is based \href{http://stackoverflow.com/questions/26198958/extracting-coefficients-and-their-standard-error-from-lme}{on this code}.
#'
#' @details Unlike \code{\link[arm]{se.coef}}, which returns the standard error
#'            for fixed and random effects separately, this function computes
#'            the standard errors for joint (sums of) random and fixed
#'            effects coefficients. Hence, \code{se} returns the appropriate
#'            standard errors for \code{\link[lme4]{coef.merMod}}.
#'            \cr \cr
#'            The standard error for the \code{\link{icc}} is based on bootstrapping,
#'            thus, the \code{nsim}-argument is required. See 'Examples'.
#'            \cr \cr
#'            \code{se} also returns the standard error of an estimate (regression
#'            coefficient) and p-value, assuming a normal distribution to compute
#'            the z-score from the p-value (formula in short: \code{b / qnorm(p / 2)}).
#'
#'
#' @examples
#' se(rnorm(n = 100, mean = 3))
#'
#' data(efc)
#' se(efc[, 1:3])
#'
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' se(fit)
#'
#' # compute standard error from regression coefficient and p-value
#' se(list(estimate = .3, p.value = .002))
#'
#' \dontrun{
#' # compute standard error of ICC for the linear mixed model
#' icc(fit)
#' se(icc(fit))
#'
#' # the standard error for the ICC can be computed manually in this way,
#' # taking the fitted model example from above
#' library(dplyr)
#' dummy <- sleepstudy %>%
#'   # generate 100 bootstrap replicates of dataset
#'   bootstrap(100) %>%
#'   # run mixed effects regression on each bootstrap replicate
#'   mutate(models = lapply(.$strap, function(x) {
#'     lmer(Reaction ~ Days + (Days | Subject), data = x)
#'   })) %>%
#'   # compute ICC for each "bootstrapped" regression
#'   mutate(icc = unlist(lapply(.$models, icc)))
#' # now compute SE and p-values for the bootstrapped ICC, values
#' # may differ from above example due to random seed
#' boot_se(dummy, icc)
#' boot_p(dummy, icc)}
#'
#'
#' @importFrom stats qnorm
#' @export
se <- function(x, nsim = 100) {
  if (is_merMod(x)) {
    # return standard error for mixed models
    return(std_merMod(x))
  } else if (any(class(x) == "icc.lme4")) {
    # we have a ICC object, so do bootstrapping and compute SE for ICC
    return(std_e_icc(x, nsim))
  } else if (is.matrix(x) || is.data.frame(x)) {
    # init return variables
    stde <- c()
    stde_names <- c()
    # iterate all columns
    for (i in seq_len(ncol(x))) {
      # get and save standard error for each variable
      # of the data frame
      stde <- c(stde, std_e_helper(x[[i]]))
      # save column name as variable name
      stde_names <- c(stde_names, colnames(x)[i])
    }
    # set names to return vector
    names(stde) <- stde_names
    # return results
    return(stde)
  } else if (is.list(x)) {
    # compute standard error from regression coefficient and p-value
    return(x$estimate / abs(stats::qnorm(x$p.value / 2)))
  } else {
    return(std_e_helper(x))
  }
}

std_e_helper <- function(x) sqrt(var(x, na.rm = TRUE) / length(stats::na.omit(x)))

#' @importFrom stats coef setNames
#' @importFrom lme4 ranef
std_merMod <- function(fit) {
  se.merMod <- list()
  # get coefficients
  cc <- stats::coef(fit)
  # get names of intercepts
  inames <- names(cc)
  # variances of fixed effects
  fixed.vars <- diag(as.matrix(lme4::vcov.merMod(fit)))
  # extract variances of conditional modes
  r1 <- lme4::ranef(fit, condVar = TRUE)
  # we may have multiple random intercepts, iterate all
  for (i in 1:length(cc)) {
    cmode.vars <- t(apply(attr(r1[[i]], "postVar"), 3, diag))
    seVals <- sqrt(sweep(cmode.vars, 2, fixed.vars, "+"))
    # add results to return list
    se.merMod[[length(se.merMod) + 1]] <- stats::setNames(as.vector(seVals[1, ]),
                                                          c("intercept_se", "slope_se"))
  }
  # set names of list
  names(se.merMod) <- inames
  return(se.merMod)
}


#' @importFrom dplyr "%>%"
#' @importFrom stats model.frame
std_e_icc <- function(x, nsim) {
  # check whether model is still in environment?
  obj.name <- attr(x, ".obj.name", exact = T)
  if (!exists(obj.name, envir = globalenv()))
    stop(sprintf("Can't find merMod-object `%s` (that was used to compute the ICC) in the environment.", obj.name), call. = F)

  # get object, see whether formulas match
  fitted.model <- globalenv()[[obj.name]]
  model.formula <- attr(x, "formula", exact = T)
  if (!identical(model.formula, formula(fitted.model)))
    stop(sprintf("merMod-object `%s` was fitted with a different formula than ICC-model", obj.name), call. = F)

  # get model family, we may have glmer
  model.family <- attr(x, "family", exact = T)

  # check for all required arguments
  if (missing(nsim) || is.null(nsim)) nsim <- 100
  # get ICC, and compute bootstrapped SE, than return both
  bstr <- bootstr_icc_se(stats::model.frame(fitted.model), nsim, model.formula, model.family)

  # now compute SE and p-values for the bootstrapped ICC
  res <- data.frame(model = obj.name,
                    icc = as.vector(x),
                    std.err = boot_se(bstr, icc)[["std.err"]],
                    p.value = boot_p(bstr, icc)[["p.value"]])
  structure(class = "se.icc.lme4", list(result = res, bootstrap_data = bstr))
}

#' @importFrom dplyr mutate
#' @importFrom lme4 lmer glmer
#' @importFrom utils txtProgressBar
bootstr_icc_se <- function(.data, nsim, formula, model.family) {
  # create progress bar
  pb <- utils::txtProgressBar(min = 1, max = nsim, style = 3)

  # generate bootstraps
  dummy <- .data %>%
    bootstrap(nsim) %>%
    dplyr::mutate(models = lapply(strap, function(x) {
      # update progress bar
      utils::setTxtProgressBar(pb, x$resample.id)
      # check model family, then compute mixed model
      if (model.family == "gaussian")
        lme4::lmer(formula, data = x)
      else
        lme4::glmer(formula, data = x, family = model.family)
    })) %>%
    # compute ICC for each "bootstrapped" regression
    dplyr::mutate(icc = unlist(lapply(models, icc)))

  # close progresss bar
  close(pb)
  return(dummy)
}