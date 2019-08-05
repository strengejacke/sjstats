#' @title Effect size statistics for anova
#' @name eta_sq
#' @description Returns the (partial) eta-squared, (partial) omega-squared,
#'   epsilon-squared statistic or Cohen's F for all terms in an anovas.
#'   \code{anova_stats()} returns a tidy summary, including all these statistics
#'   and power for each term.
#'
#' @param model A fitted anova-model of class \code{aov} or \code{anova}. Other
#'   models are coerced to \code{\link[stats]{anova}}.
#' @param partial Logical, if \code{TRUE}, the partial eta-squared is returned.
#' @param digits Number of decimal points in the returned data frame.
#' @param ci.lvl Scalar between 0 and 1. If not \code{NULL}, returns a data
#'   frame with effect sizes including lower and upper confidence intervals.
#'
#' @inheritParams bootstrap
#' @inheritParams boot_ci
#'
#' @return A data frame with the term name(s) and effect size statistics; if
#'   \code{ci.lvl} is not \code{NULL}, a data frame including lower and
#'   upper confidence intervals is returned. For \code{anova_stats()}, a tidy
#'   data frame with all statistics is returned (excluding confidence intervals).
#'
#' @details For \code{eta_sq()} (with \code{partial = FALSE}), due to
#'   non-symmetry, confidence intervals are based on bootstrap-methods. In this
#'   case, \code{n} indicates the number of bootstrap samples to be drawn to
#'   compute the confidence intervals. Confidence intervals for partial
#'   omega-squared and epsilon-squared is also based on bootstrapping.
#'   \cr \cr
#'   Since bootstrapped confidence intervals are based on the bootstrap standard error
#'   (i.e. \code{mean(x) +/- qt(.975, df = length(x) - 1) * sd(x))}, bounds of
#'   the confidence interval may be negative. Use \code{method = "quantile"} to
#'   make sure that the confidence intervals are always positive.
#'
#' @references Levine TR, Hullett CR (2002): Eta Squared, Partial Eta Squared, and Misreporting of Effect Size in Communication Research (\href{https://www.msu.edu/~levinet/eta\%20squared\%20hcr.pdf}{pdf})
#'   \cr \cr
#'   Tippey K, Longnecker MT (2016): An Ad Hoc Method for Computing Pseudo-Effect Size for Mixed Model. (\href{http://www.scsug.org/wp-content/uploads/2016/11/Ad-Hoc-Method-for-Computing-Effect-Size-for-Mixed-Models_PROCEEDINGS-UPDATE-1.pdf}{pdf})
#'
#' @examples
#' # load sample data
#' data(efc)
#'
#' # fit linear model
#' fit <- aov(
#'   c12hour ~ as.factor(e42dep) + as.factor(c172code) + c160age,
#'   data = efc
#' )
#'
#' eta_sq(fit)
#' omega_sq(fit)
#' eta_sq(fit, partial = TRUE)
#' eta_sq(fit, partial = TRUE, ci.lvl = .8)
#'
#' anova_stats(car::Anova(fit, type = 2))
#'
#' @export
eta_sq <- function(model, partial = FALSE, ci.lvl = NULL, n = 1000, method = c("dist", "quantile")) {
  method <- match.arg(method)

  if (partial)
    type <- "peta"
  else
    type <- "eta"


  es <- aov_stat(model, type = type)

  x <- data_frame(
    term = names(es),
    es = es
  )

  if (partial) {
    if (!is.null(ci.lvl) && !is.na(ci.lvl)) {
      x <- dplyr::bind_cols(x, peta_sq_ci(aov.sum = aov_stat_summary(model), ci.lvl = ci.lvl))
    }
  } else {
    if (!is.null(ci.lvl) && !is.na(ci.lvl)) {
      x <-
        es_boot_fun(
          model = model,
          type = "eta",
          ci.lvl = ci.lvl,
          n = n,
          boot.method = method
        )
    }
  }

  colnames(x)[2] <- dplyr::case_when(
    type == "eta" ~ "etasq",
    type == "peta" ~ "partial.etasq",
    TRUE ~ "effect.size"
  )

  if (!is.null(attr(es, "stratum"))) x$stratum <- attr(es, "stratum")[1:nrow(x)]

  class(x) <- c("sj_anova_stat", class(x))

  x
}



#' @importFrom purrr map_df
peta_sq_ci <- function(aov.sum, ci.lvl = .95) {
  rows <- nrow(aov.sum) - 1
  df.den <- aov.sum[["df"]][rows + 1]

  purrr::map_df(
    1:rows,
    function(.x) {
      df.num = aov.sum[.x, "df"]
      test.stat <- aov.sum[.x, "statistic"]

      if (!is.na(test.stat)) {
        ci <- partial_eta_sq_ci(
          F.value = test.stat,
          df1 = df.num,
          df2 = df.den,
          conf.level = ci.lvl
        )

        data.frame(
          conf.low = ci$LL,
          conf.high = ci$UL
        )
      } else {
        data.frame(
          conf.low = NA,
          conf.high = NA
        )
      }
    }
  )
}


#' @importFrom purrr map map_df
#' @importFrom dplyr bind_cols mutate case_when pull
#' @importFrom stats anova formula aov
#' @importFrom sjmisc rotate_df
#' @importFrom insight get_data
es_boot_fun <- function(model, type, ci.lvl, n, boot.method = "dist") {

  if (inherits(model, "anova") || is.data.frame(model)) {
    if (type == "pomega")
      stop("Objects of class `anova` or `data.frame` not supported for partial Omega squared statistics.", call. = FALSE)
    else if (type == "eta")
      stop("Objects of class `anova` or `data.frame` not supported for Eta squared statistics.", call. = FALSE)
    else
      stop("Objects of class `anova` or `data.frame` not supported.", call. = FALSE)
  }


  es <- aov_stat(model = model, type = type)

  x <- data_frame(
    term = names(es),
    es = es
  )


  # need special handling for repeated measure anova here

  if (inherits(model, "aovlist")) {

    mdata <- insight::get_data(model)
    mformula <- stats::formula(attr(model, "terms"))

    # this is a bit sloppy, but I need to catch all exceptions here
    # if we have a 1-way-anova, map() could return a column with
    # one value per row (a vector). However, if the model has more
    # covariates/factors, map() returns a list-colum with 3 values
    # per row, which need to be spread into a 3 columns data frame.

    es <- mdata %>%
      bootstrap(n = n) %>%
      dplyr::mutate(es = purrr::map(
        .data$strap,
        function(i) {
          m <- stats::aov(mformula, data = i)
          dat <- aov_stat(m, type = type)
          sjmisc::rotate_df(as.data.frame(dat))
        }
      )) %>%
      dplyr::pull(2) %>%
      purrr::map_df(~ .x) %>%
      boot_ci(ci.lvl = ci.lvl, method = boot.method)

  } else {

    mdata <- insight::get_data(model)
    mformula <- stats::formula(model)

    # this is a bit sloppy, but I need to catch all exceptions here
    # if we have a 1-way-anova, map() could return a column with
    # one value per row (a vector). However, if the model has more
    # covariates/factors, map() returns a list-colum with 3 values
    # per row, which need to be spread into a 3 columns data frame.

    es <- mdata %>%
      bootstrap(n = n) %>%
      dplyr::mutate(es = purrr::map(
        .data$strap,
        function(i) {
          m <- lm(mformula, data = i)
          dat <- aov_stat(m, type = type)
          sjmisc::rotate_df(as.data.frame(dat))
        }
      )) %>%
      dplyr::pull(2) %>%
      purrr::map_df(~ .x) %>%
      boot_ci(ci.lvl = ci.lvl, method = boot.method)
  }


  x <- dplyr::bind_cols(x, es[1:nrow(x), -1, drop = FALSE])

  colnames(x)[2] <- dplyr::case_when(
    type == "eta" ~ "etasq",
    type == "epsilon" ~ "epsilonsq",
    type == "peta" ~ "partial.etasq",
    type == "omega" ~ "omegasq",
    type == "pomega" ~ "partial.omegasq",
    TRUE ~ "effect.size"
  )

  x
}
