#' @title Standard Error for variables
#' @name se
#' @description Compute standard error for a variable, for all variables
#'                of a data frame or for joint random and fixed effects
#'                coefficients of mixed models.
#'
#' @param x (Numeric) vector, a data frame or a \code{merMod}-object
#'          as returned by the \code{\link[lme4]{lmer}}-method.
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
#' @export
se <- function(x) {
  if (is_merMod(x)) {
    return(std_merMod(x))
  } else if (is.matrix(x) || is.data.frame(x)) {
    # init return variables
    stde <- c()
    stde_names <- c()
    # iterate all columns
    for (i in 1:ncol(x)) {
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
  } else {
    return(std_e_helper(x))
  }
}

std_e_helper <- function(x) sqrt(var(x, na.rm = TRUE) / length(stats::na.omit(x)))

#' @importFrom stats coef setNames
#' @importFrom lme4 ranef
std_merMod <- function(fit) {
  # check for package availability
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package `lme4` needed for this function to work. Please install it.", call. = FALSE)
  }
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
