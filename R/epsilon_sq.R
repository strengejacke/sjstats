#' @rdname eta_sq
#' @importFrom dplyr bind_cols mutate
#' @export
epsilon_sq <- function(model, ci.lvl = NULL, n = 1000, method = c("dist", "quantile")) {
  method <- match.arg(method)
  es <- aov_stat(model, type = "epsilon")

  x <- data_frame(
    term = names(es),
    es = es
  )


  if (!is.null(ci.lvl) && !is.na(ci.lvl)) {
    x <-
      es_boot_fun(
        model = model,
        type = "epsilon",
        ci.lvl = ci.lvl,
        n = n,
        boot.method = method
      )
  }

  colnames(x)[2] <- "epsilonsq"
  if (!is.null(attr(es, "stratum"))) x$stratum <- attr(es, "stratum")[1:nrow(x)]

  class(x) <- c("sj_anova_stat", class(x))

  x
}
