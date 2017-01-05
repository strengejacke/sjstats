#' @title Expected and relative table values
#' @name table_values
#' @description This function calculates a table's cell, row and column percentages as
#'                well as expected values and returns all results as lists of tables.
#'
#' @param tab Simple \code{\link{table}} or \code{\link{ftable}} of which cell, row and column percentages
#'          as well as expected values are calculated. Tables of class \code{\link{xtabs}} and other will
#'          be coerced to \code{\link{ftable}} objects.
#' @param digits Amount of digits for the table percentage values.
#'
#' @return (Invisibly) returns a list with four tables:
#'         \enumerate{
#'          \item \code{cell} a table with cell percentages of \code{tab}
#'          \item \code{row} a table with row percentages of \code{tab}
#'          \item \code{col} a table with column percentages of \code{tab}
#'          \item \code{expected} a table with expected values of \code{tab}
#'         }
#'
#' @examples
#' tab <- table(sample(1:2, 30, TRUE), sample(1:3, 30, TRUE))
#' # show expected values
#' table_values(tab)$expected
#' # show cell percentages
#' table_values(tab)$cell
#'
#' @export
table_values <- function(tab, digits = 2) {
  # convert to ftable object
  if (!inherits(tab, "ftable")) tab <- ftable(tab)
  tab.cell <- round(100 * prop.table(tab), digits)
  tab.row <- round(100 * prop.table(tab, 1), digits)
  tab.col <- round(100 * prop.table(tab, 2), digits)
  tab.expected <- as.table(round(as.array(margin.table(tab, 1)) %*% t(as.array(margin.table(tab, 2))) / margin.table(tab)))

  # return results
  invisible(structure(class = "sjutablevalues",
                      list(cell = tab.cell,
                           row = tab.row,
                           col = tab.col,
                           expected = tab.expected)))
}


#' @title Levene-Test for One-Way-Anova
#' @name levene_test
#'
#' @description Plot results of Levene's Test for Equality of Variances for One-Way-Anova.
#'
#' @param depVar Dependent variable.
#' @param grpVar Grouping (independent) variable, as unordered factor.
#'
#' @examples
#' data(efc)
#' levene_test(efc$c12hour, efc$e42dep)
#'
#' @export
levene_test <- function(depVar, grpVar) {
  # check if grpVar is factor
  if (!is.factor(grpVar)) grpVar <- factor(grpVar)
  # remove missings
  df <- stats::na.omit(data.frame(depVar, grpVar))
  # calculate means
  means <- tapply(df$depVar, df$grpVar, mean)
  depVarNew <- abs(df$depVar - means[df$grpVar])
  message("\nLevene's Test for Homogeneity of Variances\n------------------------------------------")
  fit <- aov(depVarNew ~ df$grpVar)
  print(summary(fit))
  pval <- summary(fit)[[1]]['Pr(>F)'][1,1]
  # print "summary" of test
  message("\nConclusion:")
  if (pval > 0.05) {
    message("Groups are homogeneous. Everything's fine.\n")
  } else {
    message("Groups are not homogeneous!\n")
  }
}


# retrieve variance of random intercepts
# and residuals
lmer_var <- function(fit) {
  reva <- summary(fit)$varcor
  # retrieve only intercepts
  vars <- unlist(lapply(reva, function(x) x[[1]]))
  names(vars) <- names(reva)
  # residual variances
  if (inherits(fit, "glmerMod")) {
    # for logistic models, we use pi / 3
    resid_var <- (pi^2) / 3
  } else {
    # for linear models, we have a clear
    # residual variance
    resid_var <- attr(reva, "sc")^2
  }
  return(list('Between group variance' = vars,
              'Within group variance' = resid_var))
}


#' @importFrom stats pf
lm_pval_fstat <- function(x) {
  if (!inherits(x, "lm")) stop("Not an object of class `lm`.", call. = F)
  f <- summary(x)$fstatistic
  p <- stats::pf(f[1], f[2], f[3], lower.tail = F)
  return(as.vector(p))
}
