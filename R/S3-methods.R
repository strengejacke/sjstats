#' @importFrom nlme getData
#' @importFrom stats formula
#' @export
model.matrix.gls <- function(object, ...) {
  cbind(
    `(Intercept)` = 1,
    nlme::getData(object)[, all.vars(stats::formula(object))]
  )
}


#' @importFrom tibble tibble
#' @importFrom stats coef vcov pnorm
#' @importFrom dplyr case_when
#' @export
print.svyglm.nb <- function(x, se = c("robust", "model"), digits = 4, ...) {
  se <- match.arg(se)
  sm <- tidy_svyglm.nb(x, digits, v_se = se)[-1, -2]

  pan <- dplyr::case_when(
    sm$p.value < 0.001 ~ "<0.001 ***",
    sm$p.value < 0.01 ~ sprintf("%.*f ** ", digits, sm$p.value),
    sm$p.value < 0.05 ~ sprintf("%.*f *  ", digits, sm$p.value),
    sm$p.value < 0.1 ~ sprintf("%.*f .  ", digits, sm$p.value),
    TRUE ~  sprintf("%.*f    ", digits, sm$p.value)
  )

  sm$p.value <- pan
  print(sm, ...)

  # add dispersion parameter
  cat(sprintf("\nDispersion parameter Theta: %.*f", digits, attr(x, "nb.theta", exact = TRUE)))
  cat(sprintf("\n   Standard Error of Theta: %.*f", digits, attr(x, "nb.theta.se", exact = TRUE)))

  message(sprintf("\nShowing %s standard errors on link-scale (untransformed).", se))
}



#' @importFrom stats qnorm coef pnorm vcov
tidy_svyglm.nb <- function(x, digits = 4, v_se = c("robust", "model")) {
  v_se <- match.arg(v_se)

  if (!isNamespaceLoaded("survey"))
    requireNamespace("survey", quietly = TRUE)

  # keep original value, not rounded
  est <- stats::coef(x)
  se <- sqrt(diag(stats::vcov(x, stderr = v_se)))

  tibble::tibble(
    term = substring(names(stats::coef(x)), 5),
    estimate = round(est, digits),
    irr = round(exp(est), digits),
    std.error = round(se, digits),
    conf.low = round(exp(est - stats::qnorm(.975) * se), digits),
    conf.high = round(exp(est + stats::qnorm(.975) * se), digits),
    p.value = round(2 * stats::pnorm(abs(est / se), lower.tail = FALSE), digits)
  )
}



#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom tidyselect one_of
#' @export
model.frame.svyglm.nb <- function(formula, ...) {
  pred <- attr(formula, "nb.terms", exact = T)
  tibble::as_tibble(dplyr::select(formula$design$variables, tidyselect::one_of(pred)))
}



#' @export
family.svyglm.nb <- function(object, ...) {
  attr(object, "family", exact = TRUE)
}



#' @export
formula.svyglm.nb <- function(x, ...) {
  attr(x, "nb.formula", exact = TRUE)
}



#' @importFrom MASS glm.nb
#' @importFrom stats coef setNames predict.glm
#' @export
predict.svyglm.nb <- function(object, newdata = NULL,
                              type = c("link", "response", "terms"),
                              se.fit = FALSE, dispersion = NULL, terms = NULL,
                              na.action = na.pass, ...) {

  if (!isNamespaceLoaded("survey"))
    requireNamespace("survey", quietly = TRUE)

  fnb <- MASS::glm.nb(
    attr(object, "nb.formula", exact = TRUE),
    data = object$design$variables,
    weights = scaled.weights
  )

  cf <- stats::coef(fnb)
  names.cf <- names(cf)
  cf <- stats::coef(object)[-1]
  cf <- stats::setNames(cf, names.cf)
  fnb$coefficients <- cf

  stats::predict.glm(
    object = fnb,
    newdata = newdata,
    type = type,
    se.fit = se.fit,
    dispersion = dispersion,
    terms = terms,
    na.action = na.action,
    ...
  )
}


#' @importFrom MASS glm.nb
#' @importFrom stats coef setNames predict.glm
#' @export
residuals.svyglm.nb <- function(object, ...) {

  if (!isNamespaceLoaded("survey"))
    requireNamespace("survey", quietly = TRUE)

  fnb <- MASS::glm.nb(
    attr(object, "nb.formula", exact = TRUE),
    data = object$design$variables,
    weights = scaled.weights
  )

  y <- resp_val(fnb)
  mu <- stats::predict.glm(fnb, type = "response")
  wts <- fnb$prior.weights

  (y - mu) * sqrt(wts) / sqrt(fnb$family$variance(mu))
}


#' @importFrom stats terms formula
#' @export
terms.svyglm.nb <- function(x, ...) {

  if (!isNamespaceLoaded("survey"))
    requireNamespace("survey", quietly = TRUE)

  stats::terms(stats::formula(x), ...)
}


#' @importFrom purrr map flatten_df
#' @export
AIC.svyglm.nb <- function(object, ...) {
  ## FIXME this one just returns the AIC of the underlying glm.nb() model
  list(object, ...) %>%
    purrr::map(~ getaic(.x)) %>%
    purrr::flatten_df() %>%
    as.data.frame()
}


getaic <- function(x) {
  c(df = x$df, AIC = x$aic)
}


#' @export
deviance.svyglm.nb <- function(object, ...) {
  ## FIXME this one just returns the deviance of the underlying glm.nb() model
  object$deviance
}


#' @export
print.sjstats_r2 <- function(x, ...) {
  s3 <- NULL
  s4 <- NULL
  if (length(x) > 1) {
    if (identical(names(x[[2]]), "Nagelkerke")) {
      s1 <- "Cox & Snell's R-squared"
      s2 <- " Nagelkerke's R-squared"
    } else if (identical(names(x[[2]]), "Sums-of-Squares-r-squared")) {
      s1 <- "       R-squared (deviance)"
      s2 <- "R-squared (sums-of-squares)"
    } else if (identical(names(x[[2]]), "adj.R2")) {
      s1 <- "         R-squared"
      s2 <- "adjusted R-squared"
    } else if (identical(names(x[[2]]), "O2")) {
      s1 <- "    R-squared"
      s2 <- "Omega-squared"
    } else if (identical(names(x[[2]]), "R2(tau-11)")) {
      s1 <- "R-squared (tau-00)"
      s2 <- "R-squared (tau-11)"
      s3 <- "     Omega-squared"
      s4 <- "         R-squared"
    } else {
      return(NULL)
    }
    cat(sprintf("%s: %.4f\n%s: %.4f\n", s1, x[[1]], s2, x[[2]]))
    if (!is.null(s3)) {
      cat(sprintf("%s: %.4f\n%s: %.4f\n", s3, x[[3]], s4, x[[4]]))
    }
  } else {
    if (identical(names(x[[1]]), "D")) {
      s1 <- "Tjur's D"
    } else if (identical(names(x[[1]]), "Bayes R2")) {
      s1 <- "Bayes R2"
    } else {
      return(NULL)
    }
    cat(sprintf("%s: %.4f\n", s1, x[[1]]))
  }
}



#' @export
print.icc.lme4 <- function(x, comp, ...) {
  # print model information
  cat(sprintf("\n%s\n Family: %s (%s)\nFormula: %s\n\n",
              attr(x, "model", exact = T),
              attr(x, "family", exact = T),
              attr(x, "link", exact = T),
              paste(as.character(attr(x, "formula"))[c(2, 1, 3)], collapse = " ")))


  if (!missing(comp) && !is.null(comp) && comp == "var") {
    # get parameters
    tau.00 <- attr(x, "tau.00", exact = TRUE)
    tau.01 <- attr(x, "tau.01", exact = TRUE)
    tau.11 <- attr(x, "tau.11", exact = TRUE)
    rho.01 <- attr(x, "rho.01", exact = TRUE)

    # print within-group-variance sigma^2
    tmp <- sprintf("%.3f", attr(x, "sigma_2", exact = TRUE))
    cat(sprintf("      Within-group-variance: %8s\n", tmp))

    # print between-group-variance tau00
    for (i in seq_len(length(tau.00))) {
      tmp <- sprintf("%.3f", tau.00[i])
      cat(sprintf("     Between-group-variance: %8s (%s)\n",
                  tmp, names(tau.00)[i]))
    }

    # print random-slope-variance tau11
    for (i in seq_len(length(tau.11))) {
      tau.rs <- tau.11[i]
      # any random slope?
      if (!sjmisc::is_empty(tau.rs)) {
        tmp <- sprintf("%.3f", tau.rs)
        cat(sprintf("      Random-slope-variance: %8s (%s)\n",
                    tmp, names(tau.rs)))
      }
    }

    # print random-slope-covariance tau01
    for (i in seq_len(length(tau.01))) {
      tau.rs <- tau.01[i]
      # any random slope?
      if (!sjmisc::is_empty(tau.rs)) {
        tmp <- sprintf("%.3f", tau.rs)
        cat(sprintf(" Slope-Intercept-covariance: %8s (%s)\n",
                    tmp, names(tau.rs)))
      }
    }

    # print random-slope-correlation rho01
    for (i in seq_len(length(rho.01))) {
      rho.rs <- rho.01[i]
      # any random slope?
      if (!sjmisc::is_empty(rho.rs)) {
        tmp <- sprintf("%.3f", rho.rs)
        cat(sprintf("Slope-Intercept-correlation: %8s (%s)\n",
                    tmp, names(rho.rs)))
      }
    }
  } else {
    # get longest rand. effect name
    len <- max(nchar(names(x)))

    # print icc
    for (i in seq_len(length(x))) {
      # create info string
      infs <- sprintf("ICC (%s)", names(x[i]))
      # print info line, formatting all ICC values so they're
      # aligned properly
      cat(sprintf("%*s: %f\n",
                  len + 8,
                  infs,
                  as.vector(x[i])))
    }
  }
}


#' @importFrom tidyselect starts_with
#' @importFrom sjmisc remove_empty_cols
#' @importFrom cli cat_line
#' @importFrom crayon cyan blue red magenta green silver italic
#' @importFrom dplyr case_when
#' @export
print.icc.posterior <- function(x, ..., prob = .89, digits = 3) {
  # print model information
  cli::cat_line(crayon::italic("\n# Random Effect Variances and ICC\n"))
  cat(sprintf("Family: %s (%s)\nFormula: %s\n\n",
              attr(x, "family", exact = T),
              attr(x, "link", exact = T),
              as.character(attr(x, "formula"))[1]))

  x <- sjmisc::remove_empty_cols(x)

  get_re_col <- function(i, st) {
    dplyr::case_when(
      i == 1 ~ crayon::blue(st),
      i == 2 ~ crayon::cyan(st),
      i == 3 ~ crayon::magenta(st),
      i == 4 ~ crayon::green(st),
      i == 5 ~ crayon::red(st),
      TRUE ~ crayon::silver(st)
    )
  }

  cn <- colnames(x)
  cn.icc <- cn[tidyselect::starts_with("icc_", vars = cn)]
  cn.tau00 <- cn[tidyselect::starts_with("tau.00_", vars = cn)]

  # print icc

  for (i in seq_len(length(cn.icc))) {
    re.name <- substr(cn[i], 5, nchar(cn.icc[i]))

    cli::cat_line(get_re_col(i, sprintf("## %s", re.name)))

    # ICC
    ci <- hdi(x[[cn.icc[i]]], prob = prob)
    cli::cat_line(sprintf(
      "          ICC: %.*f (HDI %i%%: %.*f-%.*f)",
      digits,
      median(x[[cn.icc[i]]]),
      as.integer(round(prob * 100)),
      digits,
      ci[1],
      digits,
      ci[2]
    ))

    # ICC
    ci <- hdi(x[[cn.tau00[i]]], prob = prob)
    cli::cat_line(sprintf(
      "Between-group: %.*f (HDI %i%%: %.*f-%.*f)\n",
      digits,
      median(x[[cn.tau00[i]]]),
      as.integer(round(prob * 100)),
      digits,
      ci[1],
      digits,
      ci[2]
    ))
  }

  # print sigma squared

  ci <- hdi(x[["resid_var"]], prob = prob)
  infs <- crayon::red("## Residuals")
  cli::cat_line(sprintf(
    "%s\nWithin-group: %.*f (HDI %i%%: %.*f-%.*f)\n",
    infs,
    digits,
    median(x[["resid_var"]]),
    as.integer(round(prob * 100)),
    digits,
    ci[1],
    digits,
    ci[2]
  ))


  cn <- colnames(x)
  cn <- cn[tidyselect::starts_with("tau.11_", vars = cn)]

  if (!sjmisc::is_empty(cn)) cat(crayon::red("## Random-slope-variance\n"))

  # print Random-slope-variance

  for (i in seq_len(length(cn))) {
    tau.name <- substr(cn[i], 8, nchar(cn[i]))
    infs <- sprintf("%s", tau.name)
    ci <- hdi(x[[cn[i]]], prob = prob)
    cli::cat_line(sprintf(
      "%s: %.*f (HDI %i%%: %.*f-%.*f)",
      infs,
      digits,
      median(x[[cn[i]]]),
      as.integer(round(prob * 100)),
      digits,
      ci[1],
      digits,
      ci[2]
    ))
  }
}


#' @export
as.integer.sj_resample <- function(x, ...) {
  x$id
}



#' @export
as.data.frame.sj_resample <- function(x, ...) {
  x$data[x$id, , drop = FALSE]
}



#' @export
print.sj_resample <- function(x, ...) {
  n <- length(x$id)
  if (n > 12)
    id10 <- c(x$id[1:12], "...")
  else
    id10 <- x$id

  cat("<", paste0("id's of resample [", prettyNum(nrow(x$data), big.mark = ","), " x ",
                  prettyNum(ncol(x$data), big.mark = ","), "]"), "> ",
      paste(id10, collapse = ", "), "\n", sep = "")
}



#' @export
print.se.icc.lme4 <- function(x, ...) {
  cat("Standard Error of ICC\n")
  cat(sprintf("      Model: %s\n", x$result$model[[1]]))

  for (i in 1:nrow(x$result)) {
    cat(sprintf("        ICC: %.4f\n", x$result$icc[i]))
    cat(sprintf("  std. err.: %.4f\n", x$result$std.err[i]))
    cat(sprintf("    p-value: %.4f\n", x$result$p.value[i]))

    if (i < nrow(x$result)) cat("\n")
  }
}



#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @export
plot.sj_inequ_trend <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` required for plotting inequalities trends.", call. = F)
  }

  # add time indicator
  x$data$zeit <- seq_len(nrow(x$data))

  # get gather column names
  gather.cols1 <- colnames(x$data)[!colnames(x$data) %in% c("zeit", "lo", "hi")]
  gather.cols2 <- colnames(x$data)[!colnames(x$data) %in% c("zeit", "rr", "rd")]

  key_col <- "grp"
  value_col <- "y"

  # gather data to plot rr and rd
  dat1 <- tidyr::gather(x$data, !! key_col, !! value_col, !! gather.cols1)

  # gather data for raw prevalences
  dat2 <- tidyr::gather(x$data, !! key_col, !! value_col, !! gather.cols2)

  # Proper value names, for facet labels
  dat1$grp[dat1$grp == "rr"] <- "Rate Ratios"
  dat1$grp[dat1$grp == "rd"] <- "Rate Differences"

  # plot prevalences
  gp1 <- ggplot2::ggplot(dat2, ggplot2::aes_string(x = "zeit", y = "y", colour = "grp")) +
    ggplot2::geom_smooth(method = "loess", se = F) +
    ggplot2::labs(title = "Prevalance Rates for Lower and Higher SES Groups",
                  y = "Prevalances", x = "Time", colour = "") +
    ggplot2::scale_color_manual(values = c("darkblue", "darkred"), labels = c("High SES", "Low SES"))


  # plot rr and rd
  gp2 <- ggplot2::ggplot(dat1, ggplot2::aes_string(x = "zeit", y = "y", colour = "grp")) +
    ggplot2::geom_smooth(method = "loess", se = F) +
    ggplot2::facet_wrap(~grp, ncol = 1, scales = "free") +
    ggplot2::labs(title = "Proportional Change in Rate Ratios and Rate Differences",
                  colour = NULL, y = NULL, x = "Time") +
    ggplot2::guides(colour = FALSE)

  suppressMessages(graphics::plot(gp1))
  suppressMessages(graphics::plot(gp2))
}


#' @importFrom cli cat_line
#' @importFrom crayon blue cyan italic
#' @importFrom stats kruskal.test na.omit
#' @export
print.sj_mwu <- function(x, ...) {
  cli::cat_line(crayon::cyan(crayon::italic("\nMann-Whitney-U-Test\n")))
  # get data
  .dat <- x$df
  # print to console
  for (i in seq_len(nrow(.dat))) {
    # get value labels
    l1 <- .dat[i, "grp1.label"]
    l2 <- .dat[i, "grp2.label"]
    # do we have value labels?
    if (!is.null(l1) && !is.na(l1) %% !is.null(l2) && !is.na(l2)) {
      cli::cat_line(crayon::blue(sprintf("Groups %i = %s (n = %i) | %i = %s (n = %i):",
                  .dat[i, "grp1"], l1, .dat[i, "grp1.n"],
                  .dat[i, "grp2"], l2, .dat[i, "grp2.n"])))
    } else {
      cli::cat_line(crayon::blue(sprintf("Groups (%i|%i), n = %i/%i:",
                  .dat[i, "grp1"], .dat[i, "grp2"],
                  .dat[i, "grp1.n"], .dat[i, "grp2.n"])))
    }

    pval <- .dat[i, "p"]
    if (pval < 0.001) {
      pval <- 0.001
      p.string <- "<"
    } else {
      p.string <- "="
    }
    cat(sprintf("  U = %.3f, W = %.3f, p %s %.3f, Z = %.3f\n  effect-size r = %.3f\n  rank-mean(%i) = %.2f\n  rank-mean(%i) = %.2f\n\n",
                .dat[i, "u"], .dat[i, "w"], p.string, pval, .dat[i, "z"], .dat[i, "r"], .dat[i, "grp1"], .dat[i, "rank.mean.grp1"], .dat[i, "grp2"], .dat[i, "rank.mean.grp2"]))
  }

  # if we have more than 2 groups, also perfom kruskal-wallis-test
  if (length(unique(stats::na.omit(x$data$grp))) > 2) {
    cli::cat_line(crayon::cyan(crayon::italic("Kruskal-Wallis-Test\n")))
    kw <- stats::kruskal.test(x$data$x, x$data$grp)
    cat(sprintf("chi-squared = %.3f\n", kw$statistic))
    cat(sprintf("df = %i\n", kw$parameter))
    if (kw$p.value < 0.001) {
      p  <- 0.001
      p.string <- "<"
    } else {
      p <- kw$p.value
      p.string <- "="
    }
    cat(sprintf("p %s %.3f\n", p.string, p))
  }
}



#' @export
print.sj_splithalf <- function(x, ...) {
  cat(sprintf("\n   Split-Half Reliability: %.3f\n", x$splithalf))
  cat(sprintf("Spearman-Brown Adjustment: %.3f\n", x$spearmanbrown))
}



#' @export
print.sjstats_zcf <- function(x, ...) {
  cat(sprintf("   Observed zero-counts: %i\n", x$observed.zeros))
  cat(sprintf("  Predicted zero-counts: %i\n", x$predicted.zeros))
  cat(sprintf("                  Ratio: %.2f\n\n", x$ratio))

  lower <- 1 - x$tolerance
  upper <- 1 + x$tolerance

  if (x$ratio < lower)
    message("Model is underfitting zero-counts (probable zero-inflation).")
  else if (x$ratio > upper)
    message("Model is overfitting zero-counts.")
  else
    message("Model seems ok, ratio of observed and predicted zeros is within the tolerance range.")
}



#' @export
print.sjstats_ovderdisp <- function(x, ...) {
  cat("Overdispersion test\n\n")
  cat(sprintf("       dispersion ratio = %.4f\n", x$ratio))
  cat(sprintf("  Pearson's Chi-Squared = %.4f\n", x$chisq))
  cat(sprintf("                p-value = %.4f\n\n", x$p))

  if (x$p > 0.05)
    message("No overdispersion detected.")
  else
    message("Overdispersion detected.")
}



#' @export
print.sjstats_outliers <- function(x, ...) {
  print(x$result, ...)
}


#' @importFrom cli cat_line
#' @importFrom crayon cyan italic
#' @export
print.sj_xtab_stat <- function(x, ...) {
  # get length of method name, to align output
  l <- nchar(x$method)

  # is method shorter than p-value?
  if (l < 7) l <- 7

  # headline
  cli::cat_line(crayon::cyan(
    crayon::italic("\nMeasure of Association for Contingency Tables")
  ))

  # used fisher?
  if (x$fisher)
    cli::cat_line(crayon::cyan("                  (using Fisher's Exact Test)"))

  cat("\n")

  # print test statistic
  cat(sprintf("  %*s: %.4f\n", l, x$stat.name, x$statistic))
  cat(sprintf("  %*s: %.4f\n", l, x$method, x$estimate))

  # check if p <.001
  if (x$p.value < 0.001)
    cat(sprintf("  %*s: <0.001\n", l, "p-value", x$p.value))
  else
    cat(sprintf("  %*s: %.4f\n", l, "p-value", x$p.value))
}



#' @importFrom cli cat_line
#' @importFrom crayon cyan italic
#' @export
print.sjstats_pred_accuracy <- function(x, ...) {
  # headline
  cli::cat_line(crayon::cyan(crayon::italic("\nAccuracy of Model Predictions\n")))

  # statistics
  cat(sprintf("Accuracy: %.2f%%\n", 100 * x$accuracy))
  cat(sprintf("      SE: %.2f%%-points\n", 100 * x$std.error))
  cat(sprintf("  Method: %s", x$stat))
}



#' @export
print.sj_grpmean <- function(x, ...) {
  cat("\n")
  print_grpmean(x, ...)
}


#' @importFrom cli cat_line
#' @importFrom crayon blue
print_grpmean <- function(x, ...) {
  # headline
  cli::cat_line(crayon::blue(sprintf(
    "Grouped Means for %s by %s\n",
    attr(x, "dv.label", exact = TRUE),
    attr(x, "grp.label", exact = TRUE)
  )))

  # means
  print(as.data.frame(x), ...)

  # statistics
  cat(sprintf(
    "\nAnova: R2=%.3f; adj.R2=%.3f; F=%.3f; p=%.3f\n",
    attr(x, "r2", exact = TRUE),
    attr(x, "adj.r2", exact = TRUE),
    attr(x, "fstat", exact = TRUE),
    attr(x, "p.value", exact = TRUE)
  ))
}


#' @importFrom cli cat_line
#' @importFrom crayon cyan italic
#' @importFrom purrr walk
#' @export
print.sj_grpmeans <- function(x, ...) {

  cat("\n")
  purrr::walk(x, function(dat) {
    # get grouping title label
    grp <- attr(dat, "group", exact = T)

    # print title for grouping
    cli::cat_line(crayon::cyan(crayon::italic(sprintf("Grouped by:\n%s\n", grp))))

    # print grpmean-table
    print_grpmean(dat, ...)

    cat("\n\n")
  })
}


#' @export
print.sj_revar <- function(x, ...) {
  # get parameters
  xn <- names(x)
  tau.00 <- x[str_ends_with(xn, "tau.00")]
  tau.01 <- x[str_ends_with(xn, "tau.01")]
  tau.11 <- x[str_ends_with(xn, "tau.11")]
  rho.01 <- x[str_ends_with(xn, "rho.01")]
  sigma_2 <- x[str_ends_with(xn, "sigma_2")]

  # print within-group-variance sigma^2
  tmp <- sprintf("%.3f", sigma_2)
  cat(sprintf("      Within-group-variance: %8s\n", tmp))

  # print between-group-variance tau00
  for (i in seq_len(length(tau.00))) {
    tmp <- sprintf("%.3f", tau.00[i])
    cat(sprintf(
      "     Between-group-variance: %8s (%s)\n",
      tmp,
      substr(names(tau.00)[i], start = 1, stop = nchar(names(tau.00)[i]) - 7)
    ))
  }

  # print random-slope-variance tau11
  for (i in seq_len(length(tau.11))) {
    tau.rs <- tau.11[i]
    # any random slope?
    if (!sjmisc::is_empty(tau.rs)) {
      tmp <- sprintf("%.3f", tau.rs)
      cat(sprintf(
        "      Random-slope-variance: %8s (%s)\n",
        tmp,
        substr(names(tau.rs), start = 1, stop = nchar(names(tau.rs)) - 7)
      ))
    }
  }

  # print random-slope-covariance tau01
  for (i in seq_len(length(tau.01))) {
    tau.rs <- tau.01[i]
    # any random slope?
    if (!sjmisc::is_empty(tau.rs)) {
      tmp <- sprintf("%.3f", tau.rs)
      cat(sprintf(
        " Slope-Intercept-covariance: %8s (%s)\n",
        tmp,
        substr(names(tau.rs), start = 1, stop = nchar(names(tau.rs)) - 7)
      ))
    }
  }

  # print random-slope-correlation rho01
  for (i in seq_len(length(rho.01))) {
    rho.rs <- rho.01[i]
    # any random slope?
    if (!sjmisc::is_empty(rho.rs)) {
      tmp <- sprintf("%.3f", rho.rs)
      cat(sprintf(
        "Slope-Intercept-correlation: %8s (%s)\n",
        tmp,
        substr(names(rho.rs), start = 1, stop = nchar(names(rho.rs)) - 7)
      ))
    }
  }
}


#' @importFrom rlang .data
#' @importFrom sjmisc rotate_df
#' @importFrom dplyr case_when
#' @importFrom purrr map_df
#' @importFrom tibble add_column
#' @export
print.sjstats.pca_rotate <- function(x, cutoff = .1, ...) {

  xs <- attr(x, "variance", exact = TRUE)

  rn <- rownames(x)

  x <- x %>%
    round(4) %>%
    purrr::map_df(~ dplyr::case_when(
      abs(.x) < cutoff ~ "",
      TRUE ~ as.character(.x)
    )) %>%
    as.data.frame() %>%
    tibble::add_column(variable = rn, .before = 1)

  xs <- xs %>%
    round(3) %>%
    as.data.frame() %>%
    sjmisc::rotate_df()

  colnames(xs) <- sprintf("PC%i", 1:ncol(xs))
  rownames(xs) <- c("Proportion variance", "Cumulative variance", "Proportion explained", "Cumulative explained")

  print(x, quote = FALSE, ...)
  cat("\n")
  print(xs, ...)
}


#' @export
print.sjstats.pca <- function(x, ...) {

  x <- as.data.frame(round(x, 4))
  rownames(x) <- c("Standard deviation", "Eigenvalue", "Proportion variance", "Cumulative variance")

  print(x, ...)
}
