#' @rdname eta_sq
#' @importFrom effectsize omega_squared
#' @export
omega_sq <- function(model, partial = FALSE, ci.lvl = NULL) {
  out <- effectsize::omega_squared(model, partial = partial, ci = ci.lvl)

  if (isTRUE(partial)) {
    cname <- "partial.omegasq"
  } else {
    cname <- "omegasq"
  }

  out <- .fix_column_names(out, cname)
  class(out) <- c("sj_anova_stat", class(out))
  out
}
