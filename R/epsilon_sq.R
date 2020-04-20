#' @rdname eta_sq
#' @importFrom effectsize epsilon_squared
#' @export
epsilon_sq <- function(model, partial = FALSE, ci.lvl = NULL) {
  out <- effectsize::epsilon_squared(model, partial = partial, ci = ci.lvl)

  if (isTRUE(partial)) {
    cname <- "partial.epsilonsq"
  } else {
    cname <- "epsilonsq"
  }

  out <- .fix_column_names(out, cname)
  class(out) <- c("sj_anova_stat", class(out))
  out
}
