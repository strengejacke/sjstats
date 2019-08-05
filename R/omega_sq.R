#' @rdname eta_sq
#' @importFrom dplyr bind_cols mutate
#' @export
omega_sq <- function(model, partial = FALSE, ci.lvl = NULL, n = 1000, method = c("dist", "quantile")) {
  method <- match.arg(method)

  if (partial)
    type <- "pomega"
  else
    type <- "omega"


  es <- aov_stat(model, type = type)

  x <- data_frame(
    term = names(es),
    es = es
  )


  if (partial) {
    if (!is.null(ci.lvl) && !is.na(ci.lvl)) {
      x <-
        es_boot_fun(
          model = model,
          type = "pomega",
          ci.lvl = ci.lvl,
          n = n,
          boot.method = method
        )
    }
  } else {
    if (!is.null(ci.lvl) && !is.na(ci.lvl)) {
      x <- dplyr::bind_cols(x, omega_sq_ci(aov.sum = aov_stat_summary(model), ci.lvl = ci.lvl))
    }
  }

  colnames(x)[2] <- dplyr::case_when(
    type == "omega" ~ "omegasq",
    type == "pomega" ~ "partial.omegasq",
    TRUE ~ "effect.size"
  )

  if (!is.null(attr(es, "stratum"))) x$stratum <- attr(es, "stratum")[1:nrow(x)]

  class(x) <- c("sj_anova_stat", class(x))

  x
}



#' @importFrom purrr map_df
omega_sq_ci <- function(aov.sum, ci.lvl = .95) {
  rows <- nrow(aov.sum) - 1
  df.den <- aov.sum[["df"]][rows + 1]
  N <- sum(aov.sum[["df"]]) + 1

  purrr::map_df(
    1:rows,
    function(.x) {
      df.num = aov.sum[.x, "df"]
      test.stat <- aov.sum[.x, "statistic"]

      if (!is.na(test.stat)) {
        ci <- confint_ncg(
          F.value = test.stat,
          conf.level = ci.lvl,
          df.1 = df.num,
          df.2 = df.den
        )

        ci.low <- ci$Lower.Limit / (ci$Lower.Limit + N)
        ci.high <- ci$Upper.Limit / (ci$Upper.Limit + N)
      } else {
        ci.low <- ci.high <- NA
      }

      data.frame(
        conf.low = ci.low,
        conf.high = ci.high
      )
    }
  )
}
