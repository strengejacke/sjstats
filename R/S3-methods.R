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



#' @importFrom dplyr select
#' @export
model.frame.svyglm.nb <- function(formula, ...) {
  pred <- attr(formula, "nb.terms", exact = T)
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
#' @importFrom rlang .data
#' @importFrom dplyr filter slice select
#' @importFrom sjmisc var_rename trim
#' @export
print.tidy_stan <- function(x, ...) {

  insight::print_color("\n# Summary Statistics of Stan-Model\n\n", "blue")

  # check if data has certain terms, so we know if we print
  # zero inflated or multivariate response models

  zi <- string_starts_with(pattern = "b_zi_", x = x$term)
  resp.cor <- string_starts_with(pattern = "rescor__", x = x$term)
  ran.eff <- obj_has_name(x, "random.effect")
  multi.resp <- obj_has_name(x, "response")
  cumulative <- obj_has_name(x, "response.level")

  if (cumulative) x <- sjmisc::var_rename(x, response.level = "response")

  x$term <- gsub("^b_", "", x$term)
  x$term <- gsub(pattern = "^bsp_mo", replacement = "", x = x$term)
  simplex <- string_starts_with(pattern = "simo_mo", x = x$term)
  x$term <- gsub(
    pattern = "^simo_mo(.*)(\\.)(.*)(\\.)",
    replacement = "\\1 \\[\\3\\]",
    x = x$term
  )
  x$simplex <- ""
  if (!sjmisc::is_empty(simplex))
    x$simplex[simplex] <- "simplex"

  if (is.null(attr(x, "trans"))) {
    x <- get_hdi_data(x, digits = as.numeric(attr(x, "digits")))
  } else {
    x <- get_hdi_data(x, digits = as.numeric(attr(x, "digits")), ci_pattern = "ci.", ci_name = "CI")
  }

  # print zero-inflated models ----

  if (!sjmisc::is_empty(zi)) {

    # if we have random effects, make sure that any random effects
    # from zero-inlfated model are correctly handled here

    if (ran.eff) {
      zi <- union(zi, string_contains(pattern = "__zi", x = x$term))
    }

    x.zi <- dplyr::slice(x, !! zi)
    x <- dplyr::slice(x, -!! zi)

    x.zi$term <- gsub("(^b_zi_)", "", x.zi$term)
    x.zi$term <- gsub("(^zi_)", "", x.zi$term)

    x$term <- clean_term_name(x$term)
    x.zi$term <- clean_term_name(x.zi$term)

    if (ran.eff) {
      print_stan_ranef(x, zeroinf = TRUE)
    } else {
      simplex <- which(x$simplex == "simplex")

      if (!sjmisc::is_empty(simplex)) {
        x.sp <- dplyr::slice(x, !! simplex)
        x <- dplyr::slice(x, -!! simplex)
      } else
        x.sp <- NULL

      insight::print_color("## Conditional Model:\n\n", "blue")
      x <- trim_hdi(x)
      colnames(x)[1] <- ""

      x %>%
        dplyr::select(-.data$simplex) %>%
        as.data.frame() %>%
        print(..., row.names = FALSE)
      cat("\n")

      if (!is.null(x.sp))
        print_stan_simplex(x.sp)
    }


    if (ran.eff) {
      print_stan_zeroinf_ranef(x.zi)
    } else {
      simplex <- which(x$simplex == "simplex")

      if (!sjmisc::is_empty(simplex)) {
        x.sp <- dplyr::slice(x.zi, !! simplex)
        x.zi <- dplyr::slice(x.zi, -!! simplex)
      } else
        x.sp <- NULL

     insight::print_color("## Zero-Inflated Model:\n\n", "blue")
      x.zi <- trim_hdi(x.zi)
      colnames(x.zi)[1] <- ""

      x.zi %>%
        dplyr::select(-.data$simplex) %>%
        as.data.frame() %>%
        print(..., row.names = FALSE)
      cat("\n")

      if (!is.null(x.sp))
        print_stan_simplex(x.sp)
    }


    # print multivariate response models ----

  } else if (!sjmisc::is_empty(resp.cor) || multi.resp || cumulative) {

    # get the residual correlation information from data
    x.cor <- dplyr::slice(x, !! resp.cor)

    # if we have any information, remove it from remaining summary
    if (!sjmisc::is_empty(resp.cor))
      x <- dplyr::slice(x, -!! resp.cor)

    # first, print summary for each response model
    responses <- unique(x$response)
    resp.string <- ""
    if (cumulative) resp.string <- "-level"

    for (resp in responses) {

      if (ran.eff) {
        print_stan_mv_re(x, resp)
      } else {
        insight::print_color(sprintf("## Response%s: %s\n\n", resp.string, insight::print_color(resp, "red")), "blue")

        xr <- x %>%
          dplyr::filter(.data$response == !! resp) %>%
          dplyr::select(-1, -.data$simplex) %>%
          dplyr::mutate(term = clean_term_name(.data$term)) %>%
          trim_hdi() %>%
          as.data.frame()

        colnames(xr)[1] <- ""
        print(xr, ..., row.names = FALSE)

        cat("\n")
      }
    }

    # finally, if we had information on residual correlation,
    # print this as well. if "set_rescor"(FALSE)", this part
    # would be missing

    if (nrow(x.cor) > 0) {
      x.cor$term <- clean_term_name(x.cor$term)
      x.cor$term <- gsub("rescor__", "", x = x.cor$term, fixed = TRUE)
      x.cor$term <- gsub("__", "-", x = x.cor$term, fixed = TRUE)

      insight::print_color(sprintf("## Residual Correlations\n\n", resp), "cyan")

      x.cor <- x.cor %>%
        dplyr::select(-1, -.data$simplex) %>%
        sjmisc::var_rename(term = "correlation") %>%
        trim_hdi() %>%
        as.data.frame()

      colnames(x.cor)[1] <- ""
      print(x.cor, ..., row.names = FALSE)
    }
  } else if (ran.eff) {

    # print random effects models ----

    print_stan_ranef(x)

  } else {

    simplex <- which(x$simplex == "simplex")

    if (!sjmisc::is_empty(simplex)) {
      x.sp <- dplyr::slice(x, !! simplex)
      x <- dplyr::slice(x, -!! simplex)
    } else
      x.sp <- NULL

    colnames(x)[1] <- ""
    x %>%
      dplyr::select(-.data$simplex) %>%
      as.data.frame() %>%
      print(..., row.names = FALSE)

    if (!is.null(x.sp)) {
      cat("\n")
      print_stan_simplex(x.sp, ...)
    }
  }
}


print_stan_simplex <- function(x, ...) {
  insight::print_color("## Simplex Parameters:\n\n", "blue")
  x <- trim_hdi(x)

  if (colnames(x)[1] != "term")
    x <- dplyr::select(x, -1)

  colnames(x)[1] <- ""

  x %>%
    dplyr::select(-.data$simplex) %>%
    as.data.frame() %>%
    print(..., row.names = FALSE)
  cat("\n")
}


#' @importFrom sjmisc is_empty
#' @importFrom dplyr slice select filter mutate
print_stan_ranef <- function(x, zeroinf = FALSE) {
  # find fixed effects - if type = "all"
  fe <- which(x$random.effect == "")

  if (!sjmisc::is_empty(fe)) {

    if (!zeroinf)
      insight::print_color("## Fixed effects:\n\n", "blue")
    else
      insight::print_color("## Conditional Model: Fixed effects\n\n", "blue")

    # if we have fixed and random effects, get information for
    # fixed effects, and then remove these summary lines from
    # the data frame, so only random effects remain for later

    x.fe <- dplyr::slice(x, !! fe)
    x <- dplyr::slice(x, -!! fe)

    simplex <- which(x.fe$simplex == "simplex")

    if (!sjmisc::is_empty(simplex)) {
      x.sp <- dplyr::slice(x.fe, !! simplex)
      x.fe <- dplyr::slice(x.fe, -!! simplex)
    } else
      x.sp <- NULL

    x.fe$term <- clean_term_name(x.fe$term)
    x.fe <- trim_hdi(x.fe)

    colnames(x.fe)[2] <- ""

    x.fe %>%
      dplyr::select(-1, -.data$simplex) %>%
      as.data.frame() %>%
      print(row.names = FALSE)

    cat("\n")

    if (!is.null(x.sp))
      print_stan_simplex(x.sp)
  }

  # iterate all random effects
  re <- unique(x$random.effect)

  # remove random effects from zero inflated model
  re.zi <- string_contains(pattern = "__zi", x = re)
  if (!sjmisc::is_empty(re.zi)) re <- re[-re.zi]

  for (r in re) {

    if (!zeroinf)
      insight::print_color(sprintf("## Random effect %s\n\n", insight::print_color(r, "red")), "blue")
    else
      insight::print_color(sprintf("## Conditional Model: Random effect %s\n\n", insight::print_color(r, "red")), "blue")

    xr <- x %>%
      dplyr::filter(.data$random.effect == !! r) %>%
      dplyr::select(-1, -.data$simplex) %>%
      dplyr::mutate(term = clean_term_name(.data$term)) %>%
      trim_hdi() %>%
      as.data.frame()

    colnames(xr)[1] <- ""
    print(xr, row.names = FALSE)

    cat("\n")
  }
}


#' @importFrom sjmisc is_empty
trim_hdi <- function(x) {
  cn <- colnames(x)
  hdi.cols <- string_starts_with(pattern = "HDI", x = cn)

  if (!sjmisc::is_empty(hdi.cols)) {
    x <- x %>%
      purrr::map_at(hdi.cols, function(i) {
        spaces <- grep(pattern = "[ ", i, fixed = TRUE)
        if (length(spaces) == length(i))
          i <- gsub("[ ", "[", i, fixed = TRUE)

        spaces <- grep(pattern = " ]", i, fixed = TRUE)
        if (length(spaces) == length(i))
          i <- gsub(" ]", "]", i, fixed = TRUE)

        i
      }) %>%
      as.data.frame()

    colnames(x) <- cn
  }
  x
}


#' @importFrom sjmisc is_empty
#' @importFrom dplyr slice select filter mutate
print_stan_mv_re <- function(x, resp) {
  # find fixed effects - is type = "all"
  fe <- which(x$random.effect == "")

  if (!sjmisc::is_empty(fe)) {

    insight::print_color(sprintf("## Fixed effects for response: %s\n\n", insight::print_color(resp, "red")), "blue")

    x.fe <- dplyr::slice(x, !! fe)
    x <- dplyr::slice(x, -!! fe)

    xr <- x.fe %>%
      dplyr::filter(.data$response == !! resp) %>%
      dplyr::select(-1:-2, -.data$simplex) %>%
      dplyr::mutate(term = clean_term_name(.data$term)) %>%
      trim_hdi() %>%
      as.data.frame()

    colnames(xr)[1] <- ""
    print(xr, row.names = FALSE)

    cat("\n")
  }

  # iterate all random effects
  re <- unique(x$random.effect)

  for (r in re) {

    find.re <- which(x$random.effect == r & x$response == resp)

    if (!sjmisc::is_empty(find.re)) {
      insight::print_color(sprintf("## Random effect %s", insight::print_color(resp, "red")), "blue")
      insight::print_color(sprintf(" for response %s\n\n", insight::print_color(resp, "red")), "blue")

      xr <- x %>%
        dplyr::filter(.data$random.effect == !! r, .data$response == !! resp) %>%
        dplyr::select(-1:-2) %>%
        dplyr::mutate(term = clean_term_name(.data$term)) %>%
        trim_hdi() %>%
        as.data.frame()

      colnames(xr)[1] <- ""
      print(xr, row.names = FALSE)

      cat("\n")
    }
  }
}


#' @importFrom sjmisc is_empty
#' @importFrom dplyr slice select filter mutate
print_stan_zeroinf_ranef <- function(x) {
  # find fixed effects - is type = "all"
  fe <- which(x$random.effect == "")

  if (!sjmisc::is_empty(fe)) {

    insight::print_color("## Zero-Inflated Model: Fixed effects\n\n", "blue")

    # if we have fixed and random effects, get information for
    # fixed effects, and then remove these summary lines from
    # the data frame, so only random effects remain for later

    x.fe <- dplyr::slice(x, !! fe)
    x <- dplyr::slice(x, -!! fe)
    x.fe$term <- clean_term_name(x.fe$term)
    x.fe <- trim_hdi(x.fe)

    colnames(x.fe)[2] <- ""

    x.fe %>%
      dplyr::select(-1, -.data$simplex) %>%
      as.data.frame() %>%
      print(row.names = FALSE)

    cat("\n")
  }

  # iterate all random effects
  re <- unique(x$random.effect)

  # remove random effects from zero inflated model
  re.zi <- string_contains(pattern = "__zi", x = re)

  if (!sjmisc::is_empty(re.zi)) {

    re <- re[re.zi]

    x$random.effect <- gsub("__zi", "", x$random.effect, fixed = TRUE)
    x$term <- gsub("__zi", "", x$term, fixed = TRUE)
    re <- gsub("__zi", "", re, fixed = TRUE)

    for (r in re) {
      insight::print_color(sprintf("## Zero-Inflated Model: Random effect %s\n\n", insight::print_color(r, "red")), "blue")

      xr <- x %>%
        dplyr::filter(.data$random.effect == !! r) %>%
        dplyr::select(-1, -.data$simplex) %>%
        dplyr::mutate(term = clean_term_name(.data$term)) %>%
        trim_hdi() %>%
        as.data.frame()

      colnames(xr)[1] <- ""
      print(xr, row.names = FALSE)

      cat("\n")
    }
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



#' @export
print.sj_se_icc <- function(x, ...) {
  if (x$adjusted)
    cat("Standard Error of adjusted ICC\n")
  else
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
print.sj_pred_accuracy <- function(x, ...) {
  # headline
  insight::print_color("\n# Accuracy of Model Predictions\n\n", "blue")

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
  cat(.colour("blue", "\n# Causal Mediation Analysis for Stan Model\n\n"))
  cat(.colour("cyan", sprintf(
    "  Treatment: %s\n   Mediator: %s\n   Response: %s\n",
    attr(x, "treatment", exact = TRUE),
    attr(x, "mediator", exact = TRUE),
    attr(x, "response", exact = TRUE)
  )))

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

  cat(.colour("red",
    sprintf(
      "\nProportion mediated: %s%% [%s%% %s%%]\n",
      prop.med[1], prop.med[2], prop.med[3]))
  )

  if (prop.med[1] < 0)
    message("\nDirect and indirect effects have opposite directions. The proportion mediated is not meaningful.")
}


#' @importFrom purrr map_at map_df
#' @importFrom dplyr bind_cols select
get_hdi_data <- function(x, digits, ci_pattern = "hdi.", ci_name = "HDI") {
  cn <- colnames(x)
  prob <- attr(x, "prob", exact = TRUE)

  hdi.cols <- string_starts_with(pattern = ci_pattern, x = cn)

  # convert all to character, with fixed fractional part
  x <- x %>%
    purrr::map_at(hdi.cols, ~ sprintf("%.*f", digits, .x)) %>%
    as.data.frame(stringsAsFactors = FALSE)

  # left-justify term column
  x$term <- format(x$term, justify = "left")

  ci_cols <- hdi.cols[seq(1, length(hdi.cols), by = 2)]
  ci_pos <- as.vector(regexpr("_", cn))

  dummy <- purrr::map(ci_cols, function(i) {
    # get length of longest value, for proper formatting
    ml1 <- max(nchar(x[[i]]))
    ml2 <- max(nchar(x[[i + 1]]))

    tmp <- data_frame(hdi = sprintf("[%*s %*s]", ml1, x[[i]], ml2, x[[i + 1]]))

    if (ci_pos[i] < 0) {
      if (!is.null(prob))
        interv <- 100 * prob[1]
      else
        interv <- 90
    } else
      interv <- as.numeric(substr(cn[i], ci_pos[i] + 1, nchar(cn[i])))

    colnames(tmp) <- sprintf("%s(%.9g%%)", ci_name, interv)
    tmp
  }) %>%
    dplyr::bind_cols()

  # colnames(dummy) <- new_cn

  x <- dplyr::select(x, !! -hdi.cols)
  dat <- dplyr::bind_cols(x[, 1:(hdi.cols[1] - 1), drop = FALSE], dummy)

  if (ncol(x) >= hdi.cols[1])
    dat <- dplyr::bind_cols(dat, x[, hdi.cols[1]:ncol(x), drop = FALSE])

  dat
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
  cat(.colour("blue", "\n# Chi-squared Goodness-of-Fit Test\n\n"))

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
  cat(.colour("blue", "\n# Checking Model-Assumptions\n\n"))
  cat(sprintf("  Model: %s", attr(x, "formula", exact = TRUE)))

  cat(.colour("red", "\n\n                          violated    statistic\n"))

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
  cat(.colour("blue", sprintf("\n%s (%s)\n", x$method, x$alternative)))

  group <- attr(x, "group.name", exact = TRUE)
  xn <- attr(x, "x.name", exact = TRUE)
  yn <- attr(x, "y.name", exact = TRUE)

  if (!is.null(group))
    verbs <- c("of", "by")
  else
    verbs <- c("between", "and")

  st <- sprintf("# t=%.2f  df=%i  p-value=%.3f\n\n", x$statistic, as.integer(x$df), x$p.value)

  if (!is.null(yn)) {
    cat(.colour("cyan", sprintf("\n# comparison %s %s %s %s\n", verbs[1], xn, verbs[2], yn)))
  }

  cat(.colour("cyan", st))


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
    cat(sprintf("  mean of %s: %.3f [%.3f  %.3f]\n", xn, x$estimate[1], x$ci[1], x$ci[2]))
  }

  cat("\n")
}


#' @export
print.sj_wmwu <- function(x, ...) {
  cat(.colour("blue", sprintf("\n%s (%s)\n", x$method, x$alternative)))

  group <- attr(x, "group.name", exact = TRUE)
  xn <- attr(x, "x.name", exact = TRUE)

  cat(.colour("cyan", sprintf("\n# comparison of %s by %s\n", xn, group)))
  cat(.colour("cyan", sprintf("# Chisq=%.2f  df=%i  p-value=%.3f\n\n", x$statistic, as.integer(x$parameter), x$p.value)))
  cat(sprintf("  difference in mean rank score: %.3f\n\n", x$estimate))
}


#' @export
print.sj_wcor <- function(x, ...) {
  cat(.colour("blue", sprintf("\nWeighted %s\n\n", x$method)))

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


#' @export
getME.brmsfit <- function(object, name, ...) {
  rv <- NULL
  if (name == "X") {
    rv <- as.matrix(cbind(1, model_frame(object)[insight::find_predictors(object, effects = "fixed")]))
    colnames(rv)[1] = "Intercept"
  }
  rv
}


#' @export
getME.stanreg <- function(object, name, ...) {
  if (!requireNamespace("rstanarm", quietly = TRUE))
    stop("Package `rstanarm` needed for this function to work. Please install it.", call. = FALSE)

  rv <- NULL
  if (name == "X") {
    rv <- rstanarm::get_x(object)
  }
  rv
}
