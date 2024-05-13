#' @export
print.svyglm.nb <- function(x, se = c("robust", "model"), digits = 4, ...) {
  se <- match.arg(se)
  sm <- tidy_svyglm.nb(x, digits, v_se = se)[-1, -2]

  pan <- ifelse(sm$p.value < 0.001, "<0.001 ***",
    ifelse(sm$p.value < 0.01, sprintf("%.*f ** ", digits, sm$p.value), # nolint
      ifelse(sm$p.value < 0.05, sprintf("%.*f *  ", digits, sm$p.value), # nolint
        ifelse(sm$p.value < 0.1, sprintf("%.*f .  ", digits, sm$p.value), # nolint
          sprintf("%.*f    ", digits, sm$p.value)
        )
      )
    )
  )

  sm$p.value <- pan
  print(sm, ...)

  # add dispersion parameter
  cat(sprintf("\nDispersion parameter Theta: %.*f", digits, attr(x, "nb.theta", exact = TRUE)))
  cat(sprintf("\n   Standard Error of Theta: %.*f", digits, attr(x, "nb.theta.se", exact = TRUE)))

  message(sprintf("\nShowing %s standard errors on link-scale (untransformed).", se))
}


#' @export
print.svyglm.zip <- function(x, se = c("robust", "model"), digits = 4, ...) {
  se <- match.arg(se)
  sm <- tidy_svyglm.zip(x, digits, v_se = se)[-1, ]

  pan <- ifelse(sm$p.value < 0.001, "<0.001 ***",
    ifelse(sm$p.value < 0.01, sprintf("%.*f ** ", digits, sm$p.value), # nolint
      ifelse(sm$p.value < 0.05, sprintf("%.*f *  ", digits, sm$p.value), # nolint
        ifelse(sm$p.value < 0.1, sprintf("%.*f .  ", digits, sm$p.value), # nolint
          sprintf("%.*f    ", digits, sm$p.value)
        )
      )
    )
  )

  sm$p.value <- pan
  print(sm, ...)

  message(sprintf("\nShowing %s standard errors on link-scale (untransformed).", se))
}


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
    conf.low = round(exp(est - stats::qnorm(0.975) * se), digits),
    conf.high = round(exp(est + stats::qnorm(0.975) * se), digits),
    p.value = round(2 * stats::pnorm(abs(est / se), lower.tail = FALSE), digits)
  )
}


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
    conf.low = round(exp(est - stats::qnorm(0.975) * se), digits),
    conf.high = round(exp(est + stats::qnorm(0.975) * se), digits),
    p.value = round(2 * stats::pnorm(abs(est / se), lower.tail = FALSE), digits)
  )
}


#' @export
model.frame.svyglm.nb <- function(formula, ...) {
  pred <- attr(formula, "nb.terms", exact = TRUE)
  formula$design$variables[intersect(pred, colnames(formula$design$variables))]
}


#' @export
model.frame.svyglm.zip <- function(formula, ...) {
  pred <- attr(formula, "zip.terms", exact = TRUE)
  formula$design$variables[intersect(pred, colnames(formula$design$variables))]
}


#' @importFrom stats family
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


#' @export
predict.svyglm.nb <- function(object, newdata = NULL,
                              type = c("link", "response", "terms"),
                              se.fit = FALSE, dispersion = NULL, terms = NULL,
                              na.action = stats::na.pass, ...) {
  insight::check_if_installed(c("survey", "MASS"))

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


#' @export
terms.svyglm.nb <- function(x, ...) {

  if (!isNamespaceLoaded("survey"))
    requireNamespace("survey", quietly = TRUE)

  stats::terms(stats::formula(x), ...)
}


#' @export
AIC.svyglm.nb <- function(object, ...) {
  ## FIXME this one just returns the AIC of the underlying glm.nb() model
  aics <- lapply(list(object, ...), getaic)
  as.data.frame(do.call(rbind, aics))
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

  cat("<", paste0(
    "id's of resample [", prettyNum(nrow(x$data), big.mark = ","), " x ",
    prettyNum(ncol(x$data), big.mark = ","), "]"
  ), "> ",
  toString(id10), "\n",
  sep = ""
  )
}


#' @export
plot.sj_inequ_trend <- function(x, ...) {
  .data <- NULL
  insight::check_if_installed("ggplot2")

  # add time indicator
  x$data$zeit <- seq_len(nrow(x$data))

  # get gather column names
  gather.cols1 <- colnames(x$data)[!colnames(x$data) %in% c("zeit", "lo", "hi")]
  gather.cols2 <- colnames(x$data)[!colnames(x$data) %in% c("zeit", "rr", "rd")]

  # gather data to plot rr and rd
  dat1 <- datawizard::data_to_long(x$data, select = gather.cols1, names_to = "grp", values_to = "y")

  # gather data for raw prevalences
  dat2 <- datawizard::data_to_long(x$data, select = gather.cols1, names_to = "grp", values_to = "y")

  # Proper value names, for facet labels
  dat1$grp[dat1$grp == "rr"] <- "Rate Ratios"
  dat1$grp[dat1$grp == "rd"] <- "Rate Differences"

  # plot prevalences
  gp1 <- ggplot2::ggplot(dat2, ggplot2::aes_string(x = "zeit", y = "y", colour = "grp")) +
    ggplot2::geom_smooth(method = "loess", se = FALSE) +
    ggplot2::labs(title = "Prevalance Rates for Lower and Higher SES Groups",
                  y = "Prevalances", x = "Time", colour = "") +
    ggplot2::scale_color_manual(values = c("darkblue", "darkred"), labels = c("High SES", "Low SES"))


  # plot rr and rd
  gp2 <- ggplot2::ggplot(dat1, ggplot2::aes_string(x = "zeit", y = "y", colour = "grp")) +
    ggplot2::geom_smooth(method = "loess", se = FALSE) +
    ggplot2::facet_wrap(~grp, ncol = 1, scales = "free") +
    ggplot2::labs(title = "Proportional Change in Rate Ratios and Rate Differences",
                  colour = NULL, y = NULL, x = "Time") +
    ggplot2::guides(colour = "none")

  suppressMessages(graphics::plot(gp1))
  suppressMessages(graphics::plot(gp2))
}


#' @export
print.sj_xtab_stat <- function(x, ...) {
  # get length of method name, to align output
  l <- max(nchar(c(x$method, x$stat.name, "p-value", "Observations")))

  # headline
  insight::print_color("\n# Measure of Association for Contingency Tables\n", "blue")

  # used fisher?
  if (x$fisher)
    insight::print_color("                  (using Fisher's Exact Test)\n", "blue")

  cat("\n")

  # print test statistic
  cat(sprintf("  %*s: %.4f\n", l, x$stat.name, x$statistic))
  cat(sprintf("  %*s: %.4f\n", l, x$method, x$estimate))
  cat(sprintf("  %*s: %g\n", l, "df", x$df))
  cat(sprintf("  %*s: %s\n", l, "p-value", insight::format_p(x$p.value, stars = TRUE, name = NULL)))
  cat(sprintf("  %*s: %g\n", l, "Observations", x$n_obs))
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


#' @export
print.sj_anova_stat <- function(x, digits = 3, ...) {
  x$p.value <- insight::format_p(x$p.value, name = NULL)
  cat(insight::export_table(x, digits = digits, protect_integers = TRUE))
}
