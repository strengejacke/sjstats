#' @importFrom stats formula
#' @export
model.matrix.gls <- function(object, ...) {
  if (!requireNamespace("nlme"))
    stop("Package `nlme` is required, please install it first.", call. = FALSE)

  cbind(
    `(Intercept)` = 1,
    nlme::getData(object)[, all.vars(stats::formula(object))]
  )
}


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



#' @importFrom stats coef vcov pnorm
#' @importFrom dplyr case_when
#' @export
print.svyglm.zip <- function(x, se = c("robust", "model"), digits = 4, ...) {
  se <- match.arg(se)
  sm <- tidy_svyglm.zip(x, digits, v_se = se)[-1, ]

  pan <- dplyr::case_when(
    sm$p.value < 0.001 ~ "<0.001 ***",
    sm$p.value < 0.01 ~ sprintf("%.*f ** ", digits, sm$p.value),
    sm$p.value < 0.05 ~ sprintf("%.*f *  ", digits, sm$p.value),
    sm$p.value < 0.1 ~ sprintf("%.*f .  ", digits, sm$p.value),
    TRUE ~  sprintf("%.*f    ", digits, sm$p.value)
  )

  sm$p.value <- pan
  print(sm, ...)

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

  data_frame(
    term = substring(names(stats::coef(x)), 5),
    estimate = round(est, digits),
    irr = round(exp(est), digits),
    std.error = round(se, digits),
    conf.low = round(exp(est - stats::qnorm(.975) * se), digits),
    conf.high = round(exp(est + stats::qnorm(.975) * se), digits),
    p.value = round(2 * stats::pnorm(abs(est / se), lower.tail = FALSE), digits)
  )
}



#' @importFrom stats qnorm coef pnorm vcov
tidy_svyglm.zip <- function(x, digits = 4, v_se = c("robust", "model")) {
  v_se <- match.arg(v_se)

  if (!isNamespaceLoaded("survey"))
    requireNamespace("survey", quietly = TRUE)

  # keep original value, not rounded
  est <- stats::coef(x)
  se <- sqrt(diag(stats::vcov(x, stderr = v_se)))

  data_frame(
    term = substring(names(stats::coef(x)), 5),
    estimate = round(est, digits),
    std.error = round(se, digits),
    conf.low = round(exp(est - stats::qnorm(.975) * se), digits),
    conf.high = round(exp(est + stats::qnorm(.975) * se), digits),
    p.value = round(2 * stats::pnorm(abs(est / se), lower.tail = FALSE), digits)
  )
}



#' @importFrom dplyr select
#' @export
model.frame.svyglm.nb <- function(formula, ...) {
  pred <- attr(formula, "nb.terms", exact = T)
  dplyr::select(formula$design$variables, string_one_of(pattern = pred, x = colnames(formula$design$variables)))
}



#' @importFrom dplyr select
#' @export
model.frame.svyglm.zip <- function(formula, ...) {
  pred <- attr(formula, "zip.terms", exact = T)
  dplyr::select(formula$design$variables, string_one_of(pattern = pred, x = colnames(formula$design$variables)))
}



#' @export
family.svyglm.nb <- function(object, ...) {
  attr(object, "family", exact = TRUE)
}



#' @export
formula.svyglm.nb <- function(x, ...) {
  attr(x, "nb.formula", exact = TRUE)
}



#' @export
formula.svyglm.zip <- function(x, ...) {
  attr(x, "zip.formula", exact = TRUE)
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
#' @importFrom insight get_response
#' @export
residuals.svyglm.nb <- function(object, ...) {

  if (!isNamespaceLoaded("survey"))
    requireNamespace("survey", quietly = TRUE)

  fnb <- MASS::glm.nb(
    attr(object, "nb.formula", exact = TRUE),
    data = object$design$variables,
    weights = scaled.weights
  )

  y <- insight::get_response(fnb)
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


#' @importFrom insight print_color
#' @export
print.tidy_stan <- function(x, ...) {
  insight::print_color("\nSummary Statistics of Stan-Model\n\n", "blue")
  digits <- attr(x, "digits")

  for (i in x) {
    insight::print_color(paste0("# ", attr(i, "main_title")), "blue")
    cat(" ")
    insight::print_color(attr(i, "sub_title"), "red")
    cat("\n\n")
    rem <- which(colnames(i) %in% c("Parameter", "Group", "Response", "Function"))
    i <- i[, -rem]
    colnames(i)[1] <- "Parameter"
    i$ESS <- as.character(i$ESS)
    i$pd <- sprintf("%.1f%%", 100 * i$pd)
    i[] <- lapply(i, function(.j) {
      if (is.numeric(.j)) .j <- sprintf("%.*f", digits, .j)
      .j
    })
    print.data.frame(i, quote = FALSE, row.names = FALSE)
    cat("\n\n")
  }
}


#' @importFrom sjmisc trim
clean_term_name <- function(x) {
  x <- sjmisc::trim(x)
  format(x, width = max(nchar(x)))
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


#' @importFrom stats kruskal.test na.omit
#' @export
print.sj_mwu <- function(x, ...) {
  insight::print_color("\n# Mann-Whitney-U-Test\n\n", "blue")
  # get data
  .dat <- x$df
  # print to console
  for (i in seq_len(nrow(.dat))) {
    # get value labels
    l1 <- .dat[i, "grp1.label"]
    l2 <- .dat[i, "grp2.label"]
    # do we have value labels?
    if (!is.null(l1) && !is.na(l1) %% !is.null(l2) && !is.na(l2)) {
      insight::print_color(
        sprintf(
          "Groups %i = %s (n = %i) | %i = %s (n = %i):\n",
          .dat[i, "grp1"],
          l1,
          .dat[i, "grp1.n"],
          .dat[i, "grp2"],
          l2,
          .dat[i, "grp2.n"]
        ), "cyan"
      )
    } else {
      insight::print_color(
        sprintf("Groups (%i|%i), n = %i/%i:\n",
                .dat[i, "grp1"], .dat[i, "grp2"],
                .dat[i, "grp1.n"], .dat[i, "grp2.n"]),
        "cyan"
      )
    }

    pval <- .dat[i, "p"]
    if (pval < 0.001) {
      pval <- 0.001
      p.string <- "<"
    } else {
      p.string <- "="
    }

    cat(sprintf(
      "  U = %.3f, W = %.3f, p %s %.3f, Z = %.3f\n",
      .dat[i, "u"], .dat[i, "w"], p.string, pval, .dat[i, "z"]
    ))

    string_es <- "effect-size r"
    string_r <- sprintf("%.3f", .dat[i, "r"])
    string_group1 <- sprintf("rank-mean(%i)", .dat[i, "grp1"])
    string_group2 <- sprintf("rank-mean(%i)", .dat[i, "grp2"])
    string_rm1 <- sprintf("%.2f", .dat[i, "rank.mean.grp1"])
    string_rm2 <- sprintf("%.2f", .dat[i, "rank.mean.grp2"])

    space1 <- max(nchar(c(string_es, string_group1, string_group2)))
    space2 <- max(nchar(c(string_r, string_rm1, string_rm2)))

    cat(
      sprintf("  %*s = %*s\n", space1, string_es, space2 + 1, string_r),
      sprintf(" %*s = %*s\n", space1, string_group1, space2, string_rm1),
      sprintf(" %*s = %*s\n\n", space1, string_group2, space2, string_rm2)
    )
  }

  # if we have more than 2 groups, also perfom kruskal-wallis-test
  if (length(unique(stats::na.omit(x$data$grp))) > 2) {
    insight::print_color("# Kruskal-Wallis-Test\n\n", "blue")
    kw <- stats::kruskal.test(x$data$dv, x$data$grp)
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
print.sj_outliers <- function(x, ...) {
  print(x$result, ...)
}


#' @export
print.sj_xtab_stat <- function(x, ...) {
  # get length of method name, to align output
  l <- max(nchar(c(x$method, x$stat.name, "p-value")))

  # headline
  insight::print_color("\n# Measure of Association for Contingency Tables\n", "blue")

  # used fisher?
  if (x$fisher)
    insight::print_color("                  (using Fisher's Exact Test)\n", "blue")

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



#' @export
print.sj_grpmean <- function(x, ...) {
  cat("\n")
  print_grpmean(x, ...)
}


print_grpmean <- function(x, ...) {
  # headline
  insight::print_color(sprintf(
    "# Grouped Means for %s by %s\n\n",
    attr(x, "dv.label", exact = TRUE),
    attr(x, "grp.label", exact = TRUE)
  ), "blue")

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


#' @importFrom purrr walk
#' @export
print.sj_grpmeans <- function(x, ...) {

  cat("\n")
  purrr::walk(x, function(dat) {
    # get grouping title label
    grp <- attr(dat, "group", exact = T)

    # print title for grouping
    insight::print_color(sprintf("Grouped by:\n%s\n\n", grp), "cyan")

    # print grpmean-table
    print_grpmean(dat, ...)

    cat("\n\n")
  })
}


#' @export
print.sj_mediation <- function(x, digits = 2, ...) {
  insight::print_color("\n# Causal Mediation Analysis for Stan Model\n\n", "blue")
  insight::print_color(sprintf(
    "  Treatment: %s\n   Mediator: %s\n   Response: %s\n",
    attr(x, "treatment", exact = TRUE),
    attr(x, "mediator", exact = TRUE),
    attr(x, "response", exact = TRUE)
  ), "cyan")

  cat("\n")

  prop.med <- 100 * x[5, 2:4]
  x <- x[c(1, 2, 4), ]

  x$value <- format(round(x$value, digits = digits))
  x$hdi.low <- format(round(x$hdi.low, digits = digits))
  x$hdi.high <- format(round(x$hdi.high, digits = digits))
  prop.med <- format(round(prop.med, digits = digits))

  # ensure minimum width for column header
  if (max(nchar(x$value)) < 8) x$value <- format(x$value, width = 8, justify = "right")

  indent.width1 <- max(nchar(x$value)) + 17
  indent.width2 <- max(nchar(x$hdi.low)) + max(nchar(x$hdi.high)) + 4

  cat(sprintf(
    "%s%s\n",
    format("Estimate", width = indent.width1, justify = "right"),
    format(sprintf("HDI (%i%%)", as.integer(100 * attr(x, "prob", exact = TRUE))), width = indent.width2, justify = "right")
  ))

  cat(sprintf("  Direct effect: %s [%s %s]\n", x$value[1], x$hdi.low[1], x$hdi.high[1]))
  cat(sprintf("Indirect effect: %s [%s %s]\n", x$value[2], x$hdi.low[2], x$hdi.high[2]))
  cat(sprintf("   Total effect: %s [%s %s]\n", x$value[3], x$hdi.low[3], x$hdi.high[3]))

  insight::print_color(
    sprintf(
      "\nProportion mediated: %s%% [%s%% %s%%]\n",
      prop.med[1], prop.med[2], prop.med[3])
  , "red")

  if (prop.med[1] < 0)
    message("\nDirect and indirect effects have opposite directions. The proportion mediated is not meaningful.")
}


#' @export
print.sj_pval <- function(x, digits = 3, summary = FALSE, ...) {

  if (summary) {
    df.kr <- attr(x, "df.kr", exact = TRUE)
    t.kr <- attr(x, "t.kr", exact = TRUE)

    if (!is.null(df.kr)) x$df <- df.kr
    if (!is.null(t.kr)) x$statistic <- t.kr
  }

  x <- purrr::map_if(x, is.numeric, round, digits = digits)
  print.data.frame(as.data.frame(x), ..., row.names = TRUE)
}


#' @export
summary.sj_pval <- function(object, digits = 3, summary = FALSE, ...) {
  print(object, digits, summary = TRUE)
}


#' @export
print.sj_chi2gof <- function(x, ...) {
  insight::print_color("\n# Chi-squared Goodness-of-Fit Test\n\n", "blue")

  v1 <- sprintf("%.3f", x$chisq)
  v2 <- sprintf("%.3f", x$z.score)
  v3 <- sprintf("%.3f", x$p.value)

  space <- max(nchar(c(v1, v2, v3)))

  cat(sprintf("  Chi-squared: %*s\n", space, v1))
  cat(sprintf("      z-score: %*s\n", space, v2))
  cat(sprintf("      p-value: %*s\n\n", space, v3))

  if (x$p.value >= 0.05)
    message("Summary: model seems to fit well.")
  else
    message("Summary: model does not fit well.")
}


#' @export
print.sj_check_assump <- function(x, ...) {
  insight::print_color("\n# Checking Model-Assumptions\n\n", "blue")
  cat(sprintf("  Model: %s", attr(x, "formula", exact = TRUE)))

  insight::print_color("\n\n                          violated    statistic\n", "red")

  v1 <- ifelse(x$heteroskedasticity < 0.05, "yes", "no")
  v2 <- ifelse(x$multicollinearity > 4, "yes", "no")
  v3 <- ifelse(x$non.normal.resid < 0.05, "yes", "no")
  v4 <- ifelse(x$autocorrelation < 0.05, "yes", "no")

  s1 <- sprintf("p = %.3f", x$heteroskedasticity)
  s2 <- sprintf("vif = %.3f", x$multicollinearity)
  s3 <- sprintf("p = %.3f", x$non.normal.resid)
  s4 <- sprintf("p = %.3f", x$autocorrelation)

  cat(sprintf("  Heteroskedasticity      %8s  %11s\n", v1, s1))
  cat(sprintf("  Non-normal residuals    %8s  %11s\n", v3, s3))
  cat(sprintf("  Autocorrelated residuals%8s  %11s\n", v4, s4))
  cat(sprintf("  Multicollinearity       %8s  %11s\n", v2, s2))
}


#' @export
print.sj_ttest <- function(x, ...) {
  insight::print_color(sprintf("\n%s (%s)\n", x$method, x$alternative), "blue")

  group <- attr(x, "group.name", exact = TRUE)
  xn <- attr(x, "x.name", exact = TRUE)
  yn <- attr(x, "y.name", exact = TRUE)

  if (!is.null(group))
    verbs <- c("of", "by")
  else
    verbs <- c("between", "and")

  st <- sprintf("# t=%.2f  df=%i  p-value=%.3f\n\n", x$statistic, as.integer(x$df), x$p.value)

  if (!is.null(yn)) {
    insight::print_color(sprintf("\n# comparison %s %s %s %s\n", verbs[1], xn, verbs[2], yn), "cyan")
  }

  insight::print_color(st, "cyan")


  if (!is.null(yn)) {
      if (!is.null(group)) {
      l1 <- sprintf("mean in group %s", group[1])
      l2 <- sprintf("mean in group %s", group[2])
    } else {
      l1 <- sprintf("mean of %s", xn)
      l2 <- sprintf("mean of %s", yn)
    }

    l3 <- "difference of mean"

    slen <- max(nchar(c(l1, l2, l3)))

    cat(sprintf("  %s: %.3f\n", format(l1, width = slen), x$estimate[1]))
    cat(sprintf("  %s: %.3f\n", format(l2, width = slen), x$estimate[2]))
    cat(sprintf("  %s: %.3f [%.3f  %.3f]\n", format(l3, width = slen), x$estimate[1] - x$estimate[2], x$ci[1], x$ci[2]))
  } else {
    cat(sprintf("  mean of %s: %.3f [%.3f, %.3f]\n", xn, x$estimate[1], x$ci[1], x$ci[2]))
  }

  cat("\n")
}


#' @export
print.sj_wmwu <- function(x, ...) {
  insight::print_color(sprintf("\n%s (%s)\n", x$method, x$alternative), "blue")

  group <- attr(x, "group.name", exact = TRUE)
  xn <- attr(x, "x.name", exact = TRUE)

  insight::print_color(sprintf("\n# comparison of %s by %s\n", xn, group), "cyan")
  insight::print_color(sprintf("# Chisq=%.2f  df=%i  p-value=%.3f\n\n", x$statistic, as.integer(x$parameter), x$p.value), "cyan")
  cat(sprintf("  difference in mean rank score: %.3f\n\n", x$estimate))
}


#' @export
print.sj_wcor <- function(x, ...) {
  insight::print_color(sprintf("\nWeighted %s\n\n", x$method), "blue")

  if (!is.null(x$ci)) {
    cilvl <- sprintf("%.2i%%", as.integer(100 * x$ci.lvl))
    cat(sprintf("  estimate [%s CI]: %.3f [%.3f %.3f]\n", cilvl, x$estimate, x$ci[1], x$ci[2]))
    cat(sprintf("            p-value: %.3f\n\n", x$p.value))
  } else {
    cat(sprintf("  estimate: %.3f\n", x$estimate))
    cat(sprintf("   p-value: %.3f\n\n", x$p.value))
  }
}


#' @importFrom sjmisc round_num
#' @export
print.sj_anova_stat <- function(x, digits = 3, ...) {
  print.data.frame(sjmisc::round_num(x, digits), ..., row.names = TRUE)
}
