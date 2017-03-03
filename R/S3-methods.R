#' @importFrom nlme getData getCovariateFormula
#' @export
model.matrix.gls <- function(object, ...) {
  mm <- cbind(`(Intercept)` = 1,
              nlme::getData(object)[, all.vars(nlme::getCovariateFormula(object))])
  return(mm)
}

#' @importFrom nlme getResponse getData getCovariateFormula
#' @importFrom sjmisc get_label
#' @export
model.frame.gls <- function(formula, ...) {
  if (!inherits(formula, "gls")) {
    stop("`formula` needs to be an object of class `gls`.", call. = F)
    return(NULL)
  } else {
    y <- nlme::getResponse(formula)
    mf <- cbind(y, nlme::getData(formula)[, all.vars(nlme::getCovariateFormula(formula))])
    colnames(mf)[1] <- sjmisc::get_label(y, def.value = "Response")
    return(mf)
  }
}


#' @export
model.frame.gee <- function(formula, ...) {
  # get function call
  tmp <- as.list(formula$call)
  model.vars <- all.vars(tmp$formula)
  all.data <- get(deparse(tmp$data))
  all.data[, model.vars]
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
    } else {
      return(NULL)
    }
    cat(sprintf("%s: %.4f\n", s1, x[[1]]))
  }
}


#' @export
print.icc.lme4 <- function(x, comp, ...) {
  # print model information
  cat(sprintf("%s\n Family: %s (%s)\nFormula: %s\n\n",
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
  cat(sprintf("      Model: %s\n", x$result$model))
  cat(sprintf("        ICC: %.4f\n", x$result$icc))
  cat(sprintf("  std. err.: %.4f\n", x$result$std.err))
  cat(sprintf("    p-value: %.4f\n", x$result$p.value))
}


#' @importFrom tidyr gather_
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

  # gather data to plot rr and rd
  dat1 <- tidyr::gather_(x$data, key_col = "grp", value_col = "y", gather_cols = gather.cols1)

  # gather data for raw prevalences
  dat2 <- tidyr::gather_(x$data, key_col = "grp", value_col = "y", gather_cols = gather.cols2)

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


#' @export
print.sj_mwu <- function(x, ...) {
  cat("Mann-Whitney-U-Test\n")
  cat("-------------------\n")
  cat("showing statistics between groups (x|y)\n\n")
  # get data
  .dat <- x$df
  # print to console
  for (i in seq_len(nrow(.dat))) {
    # get value labels
    l1 <- .dat[i, "grp1.label"]
    l2 <- .dat[i, "grp2.label"]
    # do we have value labels?
    if (!is.null(l1) && !is.na(l1) %% !is.null(l2) && !is.na(l2)) {
      cat(sprintf("Groups %i = %s (n = %i) | %i = %s (n = %i):\n",
                  .dat[i, "grp1"], l1, .dat[i, "grp1.n"],
                  .dat[i, "grp2"], l2, .dat[i, "grp2.n"]))
    } else {
      cat(sprintf("Groups (%i|%i), n = %i/%i:\n",
                  .dat[i, "grp1"], .dat[i, "grp2"],
                  .dat[i, "grp1.n"], .dat[i, "grp2.n"]))
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
    cat("\nPerforming Kruskal-Wallis-Test\n")
    cat("------------------------------\n")
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
  cat(sprintf("   Split-Half Reliability: %.3f\n", x$splithalf))
  cat(sprintf("Spearman-Brown Adjustment: %.3f\n", x$spearmanbrown))
}

#' @export
print.sjstats_zcf <- function(x, ...) {
  cat(sprintf("   Observed zero-counts: %i\n", x$observed.zeros))
  cat(sprintf("  Predicted zero-counts: %i\n", x$predicted.zeros))
  cat(sprintf("                  Ratio: %.2f\n", x$ratio))
  cat(ifelse(x$ratio < 1,
             "Model is underfitting zero-counts.\n",
             "Model is overfitting zero-counts.\n"))
}

#' @export
print.sjstats_outliers <- function(x, ...) {
  print(x$result, ...)
}

#' @export
print.sj_xtab_stat <- function(x, ...) {
  # get length of method name, to align output
  l <- nchar(x$method)

  # is method shorter than p-value?
  if (l < 7) l <- 7

  # headline
  cat("Measure of Association for Contingency table\n")

  # used fisher?
  if (x$fisher)
    cat("                 (using Fisher's Exact Test)\n")
  else
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
