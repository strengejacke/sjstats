#' @title Effect size statistics for anova
#' @name eta_sq
#' @description Returns the (partial) eta-squared, (partial) omega-squared statistic
#'   or Cohen's F for all terms in an anovas. \code{anova_stats()} returns
#'   a tidy summary, including all these statistics and power for each term.
#'
#' @param model A fitted anova-model of class \code{aov} or \code{anova}. Other
#'   models are coerced to \code{\link[stats]{anova}}.
#' @param partial Logical, if \code{TRUE}, the partial eta-squared is returned.
#' @param digits Number of decimal points in the returned data frame.
#' @param ci.lvl Scalar between 0 and 1. If not \code{NULL}, returns a data
#'   frame with effect sizes including lower and upper confidence intervals.
#'
#' @inheritParams bootstrap
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
#'   omega-squared is also based on bootstrapping.
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
eta_sq <- function(model, partial = FALSE, ci.lvl = NULL, n = 1000) {

  if (partial)
    type <- "peta"
  else
    type <- "eta"


  es <- aov_stat(model, type = type)

  x <- data_frame(
    term = var_names(names(es)),
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
          n = n
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



#' @rdname eta_sq
#' @importFrom dplyr bind_cols mutate
#' @export
omega_sq <- function(model, partial = FALSE, ci.lvl = NULL, n = 1000) {

  if (partial)
    type <- "pomega"
  else
    type <- "omega"


  es <- aov_stat(model, type = type)

  x <- data_frame(
    term = var_names(names(es)),
    es = es
  )


  if (partial) {
    if (!is.null(ci.lvl) && !is.na(ci.lvl)) {
      x <-
      es_boot_fun(
        model = model,
        type = "pomega",
        ci.lvl = ci.lvl,
        n = n
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


#' @rdname eta_sq
#' @export
cohens_f <- function(model) {
  es <- aov_stat(model, type = "cohens.f")

  data_frame(
    term = var_names(names(es)),
    cohens.f = es
  )
}



#' @importFrom sjmisc add_columns round_num
#' @importFrom broom tidy
#' @importFrom stats anova
#' @importFrom pwr pwr.f2.test
#' @rdname eta_sq
#' @export
anova_stats <- function(model, digits = 3) {
  # get tidy summary table
  aov.sum <- aov_stat_summary(model)

  # compute all model statstics
  etasq <- aov_stat_core(aov.sum, type = "eta")
  partial.etasq <- aov_stat_core(aov.sum, type = "peta")
  omegasq <- aov_stat_core(aov.sum, type = "omega")
  partial.omegasq <- aov_stat_core(aov.sum, type = "pomega")

  # compute power for each estimate
  cohens.f <- sqrt(partial.etasq / (1 - partial.etasq))

  # bind as data frame
  as <- dplyr::bind_rows(
    data.frame(etasq, partial.etasq, omegasq, partial.omegasq, cohens.f),
    data.frame(etasq = NA, partial.etasq = NA, omegasq = NA, partial.omegasq = NA, cohens.f = NA)
  ) %>%
    sjmisc::add_columns(aov.sum)

  # get nr of terms
  nt <- nrow(as) - 1

  # finally, compute power
  power <- c(
    pwr::pwr.f2.test(u = as$df[1:nt], v = as$df[nrow(as)], f2 = as$cohens.f[1:nt]^2)[["power"]],
    NA
  )

  add_cols(as, power = power) %>%
    sjmisc::round_num(digits = digits) %>%
    as.data.frame()
}



#' @importFrom dplyr mutate
#' @importFrom rlang .data
aov_stat <- function(model, type) {
  aov.sum <- aov_stat_summary(model)
  aov.res <- aov_stat_core(aov.sum, type)

  if (obj_has_name(aov.sum, "stratum"))
    attr(aov.res, "stratum") <- aov.sum[["stratum"]]

  aov.res
}


#' @importFrom stats anova residuals
#' @importFrom broom tidy
aov_stat_summary <- function(model) {
  # check if we have a mixed model
  mm <- is_merMod(model)
  ori.model <- model

  # check that model inherits from correct class
  # else, try to coerce to anova table
  if (!inherits(model, c("aov", "anova", "anova.rms", "aovlist"))) model <- stats::anova(model)

  # get summary table
  aov.sum <- as.data.frame(broom::tidy(model))

  # for mixed models, add information on residuals
  if (mm) {
    res <- stats::residuals(ori.model)
    aov.sum <- dplyr::bind_rows(
      aov.sum,
      data_frame(
        term = "Residuals",
        df = length(res) - sum(aov.sum[["df"]]),
        sumsq = sum(res^2, na.rm = TRUE),
        meansq = mse(ori.model),
        statistic = NA
      )
    )
  }


  # check if object has sums of square
  if (!obj_has_name(aov.sum, "sumsq")) {
    stop("Model object has no sums of squares. Cannot compute effect size statistic.", call. = FALSE)
  }


  # need special handling for rms-anova
  if (inherits(model, "anova.rms"))
    colnames(aov.sum) <- c("term", "df", "sumsq", "meansq", "statistic", "p.value")

  # for car::Anova, the meansq-column might be missing, so add it manually
  if (!obj_has_name(aov.sum, "meansq"))
    aov.sum <- add_cols(aov.sum, meansq = aov.sum$sumsq / aov.sum$df, .after = "sumsq")

  aov.sum$term <- var_names(aov.sum$term)

  aov.sum
}



aov_stat_core <- function(aov.sum, type) {
  # get mean squared of residuals
  meansq.resid <- aov.sum[["meansq"]][nrow(aov.sum)]
  # get total sum of squares
  ss.total <- sum(aov.sum[["sumsq"]])
  # get sum of squares of residuals
  ss.resid <- aov.sum[["sumsq"]][nrow(aov.sum)]

  # number of terms in model
  n_terms <- nrow(aov.sum) - 1

  # number of observations
  N <- sum(aov.sum[["df"]]) + 1


  if (type == "omega") {
    # compute omega squared for each model term
    aovstat <- purrr::map_dbl(1:n_terms, function(x) {
      ss.term <- aov.sum[["sumsq"]][x]
      df.term <- aov.sum[["df"]][x]
      (ss.term - df.term * meansq.resid) / (ss.total + meansq.resid)
    })
  } else if (type == "pomega") {
    # compute partial omega squared for each model term
    aovstat <- purrr::map_dbl(1:n_terms, function(x) {
      df.term <- aov.sum[["df"]][x]
      meansq.term <- aov.sum[["meansq"]][x]
      (df.term * (meansq.term - meansq.resid)) / (df.term * meansq.term + (N - df.term) * meansq.resid)
    })
  } else if (type == "eta") {
    # compute eta squared for each model term
    aovstat <-
      purrr::map_dbl(1:n_terms, ~ aov.sum[["sumsq"]][.x] / sum(aov.sum[["sumsq"]]))
  } else if (type %in% c("cohens.f", "peta")) {
    # compute partial eta squared for each model term
    aovstat <-
      purrr::map_dbl(1:n_terms, ~ aov.sum[["sumsq"]][.x] / (aov.sum[["sumsq"]][.x] + ss.resid))
  }

  # compute Cohen's F
  if (type == "cohens.f") aovstat <- sqrt(aovstat / (1 - aovstat))

  # give values names of terms
  names(aovstat) <- aov.sum[["term"]][1:n_terms]

  aovstat
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


#' @importFrom broom tidy
#' @importFrom purrr map map_df
#' @importFrom dplyr bind_cols mutate case_when pull
#' @importFrom stats anova formula aov
#' @importFrom sjmisc rotate_df
es_boot_fun <- function(model, type, ci.lvl, n) {

  es <- aov_stat(model = model, type = type)

  x <- data_frame(
    term = var_names(names(es)),
    es = es
  )


  # need special handling for repeated measure anova here

  if (inherits(model, "aovlist")) {

    mdata <- sjstats::model_frame(model)
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
      boot_ci()

  } else {

    mdata <- sjstats::model_frame(model)
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
      boot_ci()
  }


  x <- dplyr::bind_cols(x, es[1:nrow(x), -1])

  colnames(x)[2] <- dplyr::case_when(
    type == "eta" ~ "etasq",
    type == "peta" ~ "partial.etasq",
    type == "omega" ~ "omegasq",
    type == "pomega" ~ "partial.omegasq",
    TRUE ~ "effect.size"
  )

  x
}
