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
