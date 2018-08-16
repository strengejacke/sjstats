#' @importFrom nlme getData
#' @importFrom stats formula
#' @export
model.matrix.gls <- function(object, ...) {
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


#' @importFrom purrr map_chr map2_chr
#' @export
print.sj_r2 <- function(x, digits = 3, ...) {
  cat("\nR-Squared for Generalized Linear Mixed Model\n\n")
  print_icc_and_r2(x, digits, ...)
}


#' @importFrom purrr map_chr map2_chr
#' @export
print.sj_icc <- function(x, digits = 4, ...) {
  cat("\nIntra-Class Correlation Coefficient for Generalized Linear Mixed Model\n\n")
  print_icc_and_r2(x, digits, ...)
}


print_icc_and_r2 <- function(x, digits, ...) {
  # print model information
  cat(crayon::blue(
    sprintf("Family : %s (%s)\nFormula: %s\n\n",
            attr(x, "family", exact = T),
            attr(x, "link", exact = T),
            paste(as.character(attr(x, "formula"))[c(2, 1, 3)], collapse = " ")
    )))

  labels <- purrr::map_chr(x, ~ names(.x))
  width <- max(nchar(labels))

  out <- purrr::map2_chr(
    x, labels,
    ~ sprintf("%*s: %.*f", width, .y, digits, .x)
  )

  cat(paste0(out, collapse = "\n"))
  cat("\n\n")
}


#' @export
print.sj_icc_merMod <- function(x, comp, ...) {
  # print model information
  cat(sprintf("\n%s\n\n", attr(x, "model", exact = T)))

  cat(crayon::blue(
    sprintf("Family : %s (%s)\nFormula: %s\n\n",
            attr(x, "family", exact = T),
            attr(x, "link", exact = T),
            paste(as.character(attr(x, "formula"))[c(2, 1, 3)], collapse = " ")
    )))

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
      cat(sprintf("     Between-group-variance: %8s (%s)\n", tmp, names(tau.00)[i]))
    }

    # print random-slope-variance tau11
    for (i in seq_len(length(tau.11))) {
      tau.rs <- tau.11[i]
      # any random slope?
      if (!sjmisc::is_empty(tau.rs)) {
        tmp <- sprintf("%.3f", tau.rs)
        cat(sprintf("      Random-slope-variance: %8s (%s)\n", tmp, names(tau.rs)))
      }
    }

    # print random-slope-covariance tau01
    for (i in seq_len(length(tau.01))) {
      tau.rs <- tau.01[i]
      # any random slope?
      if (!sjmisc::is_empty(tau.rs)) {
        tmp <- sprintf("%.3f", tau.rs)
        cat(sprintf(" Slope-Intercept-covariance: %8s (%s)\n", tmp, names(tau.rs)))
      }
    }

    # print random-slope-correlation rho01
    for (i in seq_len(length(rho.01))) {
      rho.rs <- rho.01[i]
      # any random slope?
      if (!sjmisc::is_empty(rho.rs)) {
        tmp <- sprintf("%.3f", rho.rs)
        cat(sprintf("Slope-Intercept-correlation: %8s (%s)\n", tmp, names(rho.rs)))
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
      cat(sprintf("%*s: %.4f\n",
                  len + 8,
                  infs,
                  as.vector(x[i])))
    }
  }
}


#' @importFrom rlang .data
#' @importFrom dplyr filter slice select
#' @importFrom crayon blue cyan red
#' @importFrom sjmisc var_rename trim
#' @export
print.tidy_stan <- function(x, ...) {

  cat(crayon::blue("\n# Summary Statistics of Stan-Model\n\n"))

  # check if data has certain terms, so we know if we print
  # zero inflated or multivariate response models

  zi <- string_starts_with(pattern = "b_zi_", x = x$term)
  resp.cor <- string_starts_with(pattern = "rescor__", x = x$term)
  ran.eff <- obj_has_name(x, "random.effect")
  multi.resp <- obj_has_name(x, "response")
  cumulative <- obj_has_name(x, "response.level")

  if (cumulative) x <- sjmisc::var_rename(x, response.level = "response")

  x$term <- gsub("b_", "", x$term, fixed = TRUE)

  x <- get_hdi_data(x, digits = as.numeric(attr(x, "digits")))


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
      cat(crayon::blue("## Conditional Model:\n\n"))
      x <- trim_hdi(x)
      colnames(x)[1] <- ""

      x %>%
        as.data.frame() %>%
        print(..., row.names = FALSE)
      cat("\n")
    }


    if (ran.eff) {
      print_stan_zeroinf_ranef(x.zi)
    } else {
      cat(crayon::blue("## Zero-Inflated Model:\n\n"))
      x.zi <- trim_hdi(x.zi)
      colnames(x.zi)[1] <- ""

      x.zi %>%
        as.data.frame() %>%
        print(..., row.names = FALSE)
      cat("\n")
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
        cat(crayon::blue(sprintf("## Response%s: %s\n\n", resp.string, crayon::red(resp))))

        xr <- x %>%
          dplyr::filter(.data$response == !! resp) %>%
          dplyr::select(-1) %>%
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

      cat(crayon::cyan(sprintf("## Residual Correlations\n\n", resp)))

      x.cor <- x.cor %>%
        dplyr::select(-1) %>%
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
    colnames(x)[1] <- ""
    x %>%
      as.data.frame() %>%
      print(..., row.names = FALSE)
  }
}


#' @importFrom sjmisc is_empty
#' @importFrom crayon blue red
#' @importFrom dplyr slice select filter mutate
print_stan_ranef <- function(x, zeroinf = FALSE) {
  # find fixed effects - is type = "all"
  fe <- which(x$random.effect == "")

  if (!sjmisc::is_empty(fe)) {

    if (!zeroinf)
      cat(crayon::blue("## Fixed effects:\n\n"))
    else
      cat(crayon::blue("## Conditional Model: Fixed effects\n\n"))

    # if we have fixed and random effects, get information for
    # fixed effects, and then remove these summary lines from
    # the data frame, so only random effects remain for later

    x.fe <- dplyr::slice(x, !! fe)
    x <- dplyr::slice(x, -!! fe)
    x.fe$term <- clean_term_name(x.fe$term)
    x.fe <- trim_hdi(x.fe)

    colnames(x.fe)[2] <- ""

    x.fe %>%
      dplyr::select(-1) %>%
      as.data.frame() %>%
      print(row.names = FALSE)

    cat("\n")
  }

  # iterate all random effects
  re <- unique(x$random.effect)

  # remove random effects from zero inflated model
  re.zi <- string_contains(pattern = "__zi", x = re)
  if (!sjmisc::is_empty(re.zi)) re <- re[-re.zi]

  for (r in re) {

    if (!zeroinf)
      cat(crayon::blue(sprintf("## Random effect %s\n\n", crayon::red(r))))
    else
      cat(crayon::blue(sprintf("## Conditional Model: Random effect %s\n\n", crayon::red(r))))

    xr <- x %>%
      dplyr::filter(.data$random.effect == !! r) %>%
      dplyr::select(-1) %>%
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
#' @importFrom crayon blue red
#' @importFrom dplyr slice select filter mutate
print_stan_mv_re <- function(x, resp) {
  # find fixed effects - is type = "all"
  fe <- which(x$random.effect == "")

  if (!sjmisc::is_empty(fe)) {

    cat(crayon::blue(sprintf("## Fixed effects for response: %s\n\n", crayon::red(resp))))

    x.fe <- dplyr::slice(x, !! fe)
    x <- dplyr::slice(x, -!! fe)

    xr <- x.fe %>%
      dplyr::filter(.data$response == !! resp) %>%
      dplyr::select(-1:-2) %>%
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
      cat(crayon::blue(sprintf("## Random effect %s", crayon::red(r))))
      cat(crayon::blue(sprintf(" for response %s\n\n", crayon::red(resp))))

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
#' @importFrom crayon blue red
#' @importFrom dplyr slice select filter mutate
print_stan_zeroinf_ranef <- function(x) {
  # find fixed effects - is type = "all"
  fe <- which(x$random.effect == "")

  if (!sjmisc::is_empty(fe)) {

    cat(crayon::blue("## Zero-Inflated Model: Fixed effects\n\n"))

    # if we have fixed and random effects, get information for
    # fixed effects, and then remove these summary lines from
    # the data frame, so only random effects remain for later

    x.fe <- dplyr::slice(x, !! fe)
    x <- dplyr::slice(x, -!! fe)
    x.fe$term <- clean_term_name(x.fe$term)
    x.fe <- trim_hdi(x.fe)

    colnames(x.fe)[2] <- ""

    x.fe %>%
      dplyr::select(-1) %>%
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
      cat(crayon::blue(sprintf("## Zero-Inflated Model: Random effect %s\n\n", crayon::red(r))))

      xr <- x %>%
        dplyr::filter(.data$random.effect == !! r) %>%
        dplyr::select(-1) %>%
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


#' @importFrom sjmisc remove_empty_cols
#' @importFrom crayon cyan blue red magenta green silver
#' @importFrom dplyr case_when
#' @export
print.sj_icc_brms <- function(x, digits = 2, ...) {
  # print model information
  cat("\n# Random Effect Variances and ICC\n\n")
  cat(sprintf(
    "Family: %s (%s)\n\n",
    attr(x, "family", exact = T),
    attr(x, "link", exact = T)
  ))

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

  prob <- attr(x, "prob", exact = TRUE)
  cn <- names(x)

  # print icc

  for (i in seq_len(length(cn))) {
    re.name <- substr(cn[i], 5, nchar(cn[i]))

    cat(get_re_col(i, sprintf("## %s\n", re.name)))

    icc.val <- sprintf("%.*f", digits, x[i])
    tau.val <- sprintf("%.*f", digits, attr(x, "tau.00", exact = TRUE)[i])
    ml <- max(nchar(icc.val), nchar(tau.val))

    ci.icc <- attr(x, "hdi.icc", exact = TRUE)[[i]]
    ci.icc.lo <- sprintf("%.*f", digits, ci.icc[1])
    ci.icc.hi <- sprintf("%.*f", digits, ci.icc[2])

    ci.tau <- attr(x, "hdi.tau.00", exact = TRUE)[[i]]
    ci.tau.lo <- sprintf("%.*f", digits, ci.tau[1])
    ci.tau.hi <- sprintf("%.*f", digits, ci.tau[2])

    ml.ci <- max(nchar(ci.icc.lo), nchar(ci.tau.lo))
    mh.ci <- max(nchar(ci.icc.hi), nchar(ci.tau.hi))

    # ICC
    cat(sprintf(
      "          ICC: %*s  HDI %i%%: [%*s %*s]\n",
      ml,
      icc.val,
      as.integer(round(prob * 100)),
      ml.ci,
      ci.icc.lo,
      mh.ci,
      ci.icc.hi
    ))

    # Tau00
    cat(sprintf(
      "Between-group: %*s  HDI %i%%: [%*s %*s]\n\n",
      ml,
      tau.val,
      as.integer(round(prob * 100)),
      ml.ci,
      ci.tau.lo,
      mh.ci,
      ci.tau.hi
    ))
  }

  # print sigma squared

  ci <- attr(x, "hdi.sigma_2", exact = TRUE)
  infs <- crayon::red("## Residuals")
  cat(sprintf(
    "%s\nWithin-group: %.*f  HDI %i%%: [%.*f %.*f]\n",
    infs,
    digits,
    attr(x, "sigma_2", exact = TRUE),
    as.integer(round(prob * 100)),
    digits,
    ci[1],
    digits,
    ci[2]
  ))


  rsv <- attr(x, "tau.11", exact = TRUE)
  if (!is.null(rsv)) cat(crayon::red("\n## Random-slope-variance\n"))

  # print Random-slope-variance

  for (i in seq_len(length(rsv))) {
    infs <- sprintf("%s", substr(names(rsv[i]), 8, nchar(names(rsv[i]))))
    ci <- attr(x, "hdi.tau.11")[[i]]
    cat(sprintf(
      "%s: %.*f  HDI %i%%: [%.*f %.*f]\n",
      infs,
      digits,
      rsv[i],
      as.integer(round(prob * 100)),
      digits,
      ci[1],
      digits,
      ci[2]
    ))
  }
}


#' @importFrom sjmisc remove_empty_cols
#' @importFrom crayon cyan blue red magenta green silver
#' @importFrom dplyr case_when
#' @export
print.sj_icc_stanreg <- function(x, digits = 2, ...) {
  # print model information
  cat("\n# Random Effect Variances and ICC\n\n")
  cat(sprintf(
    "Family: %s (%s)\n\n",
    attr(x, "family", exact = T),
    attr(x, "link", exact = T)
  ))

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

  prob <- attr(x, "prob", exact = TRUE)
  cn <- names(x)

  # print icc

  for (i in seq_len(length(cn))) {

    cat(get_re_col(i, sprintf("## %s\n", cn[i])))

    icc.val <- sprintf("%.*f", digits, x[i])
    tau.val <- sprintf("%.*f", digits, attr(x, "tau.00", exact = TRUE)[i])
    ml <- max(nchar(icc.val), nchar(tau.val))

    ci.icc <- attr(x, "hdi.icc", exact = TRUE)[[i]]
    ci.icc.lo <- sprintf("%.*f", digits, ci.icc[1])
    ci.icc.hi <- sprintf("%.*f", digits, ci.icc[2])

    ci.tau <- attr(x, "hdi.tau.00", exact = TRUE)[[i]]
    ci.tau.lo <- sprintf("%.*f", digits, ci.tau[1])
    ci.tau.hi <- sprintf("%.*f", digits, ci.tau[2])

    ml.ci <- max(nchar(ci.icc.lo), nchar(ci.tau.lo))
    mh.ci <- max(nchar(ci.icc.hi), nchar(ci.tau.hi))

    # ICC
    cat(sprintf(
      "          ICC: %*s  HDI %i%%: [%*s %*s]\n",
      ml,
      icc.val,
      as.integer(round(prob * 100)),
      ml.ci,
      ci.icc.lo,
      mh.ci,
      ci.icc.hi
    ))

    # Tau00
    cat(sprintf(
      "Between-group: %*s  HDI %i%%: [%*s %s]\n\n",
      ml,
      tau.val,
      as.integer(round(prob * 100)),
      ml.ci,
      ci.tau.lo,
      ci.tau.hi
    ))
  }

  # print sigma squared

  ci <- attr(x, "hdi.sigma_2", exact = TRUE)
  infs <- crayon::red("## Residuals")
  cat(sprintf(
    "%s\nWithin-group: %.*f  HDI %i%%: [%.*f %.*f]\n",
    infs,
    digits,
    attr(x, "sigma_2", exact = TRUE),
    as.integer(round(prob * 100)),
    digits,
    ci[1],
    digits,
    ci[2]
  ))


  rsv <- attr(x, "tau.11", exact = TRUE)
  if (!is.null(rsv)) cat(crayon::cyan("\n## Random-slope-variance\n"))

  # print Random-slope-variance

  for (i in seq_len(length(rsv))) {
    infs <- sprintf("%s", names(rsv[i]))
    ci <- attr(x, "hdi.tau.11")[[i]]
    cat(sprintf(
      "%s: %.*f  HDI %i%%: [%.*f %.*f]\n",
      infs,
      digits,
      rsv[i],
      as.integer(round(prob * 100)),
      digits,
      ci[1],
      digits,
      ci[2]
    ))
  }


  rsicov <- attr(x, "tau.01", exact = TRUE)
  rsicor <- attr(x, "rho.01", exact = TRUE)
  if (!is.null(rsicov)) cat(crayon::cyan("\n## Random-slope-intercept covariances\n"))

  # print Random-slope-variance

  for (i in seq_len(length(rsicov))) {
    infs <- sprintf("%s", gsub("^Sigma\\[(.*),\\(Intercept\\)\\]" , "\\1", names(rsicov[i])))

    cov.val <- sprintf("%.*f", digits, rsicov[i])
    cor.val <- sprintf("%.*f", digits, rsicor[i])

    ml <- max(nchar(cov.val), nchar(cor.val))

    ci.cov <- attr(x, "hdi.tau.01")[[i]]
    ci.cor <- attr(x, "hdi.rho.01")[[i]]

    ci.cov.lo <- sprintf("%.*f", digits, ci.cov[1])
    ci.cov.hi <- sprintf("%.*f", digits, ci.cov[2])

    ci.cor.lo <- sprintf("%.*f", digits, ci.cor[1])
    ci.cor.hi <- sprintf("%.*f", digits, ci.cor[2])

    ml.ci <- max(nchar(ci.cov.lo), nchar(ci.cor.lo))
    mh.ci <- max(nchar(ci.cov.hi), nchar(ci.cov.hi))

    cat(sprintf(
      " Covariance %s: %*s  HDI %i%%: [%*s %*s]\n",
      infs,
      ml,
      cov.val,
      as.integer(round(prob * 100)),
      ml.ci,
      ci.cov.lo,
      mh.ci,
      ci.cov.hi
    ))

    infs <- sprintf("%s", gsub("^Sigma\\[(.*),\\(Intercept\\)\\]" , "\\1", names(rsicor[i])))

    cat(sprintf(
      "Correlation %s: %*s  HDI %i%%: [%*s %*s]\n",
      infs,
      ml,
      cor.val,
      as.integer(round(prob * 100)),
      ml.ci,
      ci.cor.lo,
      mh.ci,
      ci.cor.hi
    ))
  }
}


#' @export
print.icc_ppd <- function(x, digits = 2, ...) {
  # print model information
  cat("\n# Random Effect Variances and ICC\n\n")

  reform <- attr(x, "re.form", exact = TRUE)
  if (is.null(reform))
    reform <- "all random effects"
  else
    reform <- deparse(reform)

  cat(sprintf(
    "Family: %s (%s)\nConditioned on: %s\n\n",
    attr(x, "family", exact = T),
    attr(x, "link", exact = T),
    reform
  ))

  prob <- attr(x, "prob", exact = TRUE)

  cat(crayon::blue("## Variance Ratio (comparable to ICC)\n"))

  icc.val <- sprintf("%.*f", digits, x["icc"])

  ci.icc <- attr(x, "hdi.icc", exact = TRUE)
  ci.icc.lo <- sprintf("%.*f", digits, ci.icc[1])
  ci.icc.hi <- sprintf("%.*f", digits, ci.icc[2])

  # ICC
  cat(sprintf(
    "Ratio: %s  HDI %i%%: [%s %s]\n",
    icc.val,
    as.integer(round(prob * 100)),
    ci.icc.lo,
    ci.icc.hi
  ))

  cat(crayon::blue("\n## Variances of Posterior Predicted Distribution\n"))

  null.model <- sprintf("%.*f", digits, x["tau.00"])

  ci.null <- attr(x, "hdi.tau.00", exact = TRUE)
  ci.null.lo <- sprintf("%.*f", digits, ci.null[1])
  ci.null.hi <- sprintf("%.*f", digits, ci.null[2])

  full.model <- sprintf("%.*f", digits, x["total.var"])

  ci.full <- attr(x, "hdi.total", exact = TRUE)
  ci.full.lo <- sprintf("%.*f", digits, ci.full[1])
  ci.full.hi <- sprintf("%.*f", digits, ci.full[2])

  ml <- max(nchar(null.model), nchar(full.model))
  ml.ci <- max(nchar(ci.full.lo), nchar(ci.null.lo))
  mh.ci <- max(nchar(ci.full.hi), nchar(ci.null.hi))

  # Conditioned on fixed effects
  cat(sprintf(
    "Conditioned on fixed effects: %*s  HDI %i%%: [%*s %*s]\n",
    ml,
    null.model,
    as.integer(round(prob * 100)),
    ml.ci,
    ci.null.lo,
    mh.ci,
    ci.null.hi
  ))

  # Conditioned on random effects
  cat(sprintf(
    "Conditioned on rand. effects: %*s  HDI %i%%: [%*s %*s]\n",
    ml,
    full.model,
    as.integer(round(prob * 100)),
    ml.ci,
    ci.full.lo,
    mh.ci,
    ci.full.hi
  ))

  cat(crayon::red("\n## Difference in Variances\n"))

  res <- sprintf("%.*f", digits, x["resid.var"])

  ci.res <- attr(x, "hdi.resid", exact = TRUE)
  ci.res.lo <- sprintf("%.*f", digits, ci.res[1])
  ci.res.hi <- sprintf("%.*f", digits, ci.res[2])

  # ICC
  cat(sprintf(
    "Difference: %s  HDI %i%%: [%s %s]\n",
    res,
    as.integer(round(prob * 100)),
    ci.res.lo,
    ci.res.hi
  ))
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


#' @importFrom crayon blue
#' @importFrom stats kruskal.test na.omit
#' @export
print.sj_mwu <- function(x, ...) {
  cat(crayon::blue("\n# Mann-Whitney-U-Test\n\n"))
  # get data
  .dat <- x$df
  # print to console
  for (i in seq_len(nrow(.dat))) {
    # get value labels
    l1 <- .dat[i, "grp1.label"]
    l2 <- .dat[i, "grp2.label"]
    # do we have value labels?
    if (!is.null(l1) && !is.na(l1) %% !is.null(l2) && !is.na(l2)) {
      cat(crayon::cyan(
        sprintf(
          "Groups %i = %s (n = %i) | %i = %s (n = %i):\n",
          .dat[i, "grp1"],
          l1,
          .dat[i, "grp1.n"],
          .dat[i, "grp2"],
          l2,
          .dat[i, "grp2.n"]
        )
      ))
    } else {
      cat(crayon::cyan(
        sprintf("Groups (%i|%i), n = %i/%i:\n",
                .dat[i, "grp1"], .dat[i, "grp2"],
                .dat[i, "grp1.n"], .dat[i, "grp2.n"])
      ))
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
    cat(crayon::blue("# Kruskal-Wallis-Test\n\n"))
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



#' @importFrom crayon blue
#' @export
print.sj_splithalf <- function(x, ...) {
  cat(crayon::blue("\n# Internal Consistency\n\n"))
  cat(sprintf("   Split-Half Reliability: %.3f\n", x$splithalf))
  cat(sprintf("Spearman-Brown Adjustment: %.3f\n", x$spearmanbrown))
}



#' @importFrom crayon blue
#' @export
print.sj_zcf <- function(x, ...) {
  cat(crayon::blue("\n# Zero-Count overfitting\n\n"))
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



#' @importFrom crayon blue
#' @export
print.sj_ovderdisp <- function(x, ...) {
  cat(crayon::blue("\n# Overdispersion test\n\n"))
  cat(sprintf("       dispersion ratio = %.4f\n", x$ratio))
  cat(sprintf("  Pearson's Chi-Squared = %.4f\n", x$chisq))
  cat(sprintf("                p-value = %.4f\n\n", x$p))

  if (x$p > 0.05)
    message("No overdispersion detected.")
  else
    message("Overdispersion detected.")
}



#' @export
print.sj_outliers <- function(x, ...) {
  print(x$result, ...)
}


#' @importFrom crayon blue
#' @export
print.sj_xtab_stat <- function(x, ...) {
  # get length of method name, to align output
  l <- nchar(x$method)

  # is method shorter than p-value?
  if (l < 7) l <- 7

  # headline
  cat(crayon::blue("\n# Measure of Association for Contingency Tables\n"))

  # used fisher?
  if (x$fisher)
    cat(crayon::blue("                  (using Fisher's Exact Test)\n"))

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



#' @importFrom crayon blue
#' @export
print.sj_pred_accuracy <- function(x, ...) {
  # headline
  cat(crayon::blue("\n# Accuracy of Model Predictions\n\n"))

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


#' @importFrom crayon blue
print_grpmean <- function(x, ...) {
  # headline
  cat(crayon::blue(sprintf(
    "# Grouped Means for %s by %s\n\n",
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


#' @importFrom crayon cyan
#' @importFrom purrr walk
#' @export
print.sj_grpmeans <- function(x, ...) {

  cat("\n")
  purrr::walk(x, function(dat) {
    # get grouping title label
    grp <- attr(dat, "group", exact = T)

    # print title for grouping
    cat(crayon::cyan(sprintf("Grouped by:\n%s\n\n", grp)))

    # print grpmean-table
    print_grpmean(dat, ...)

    cat("\n\n")
  })
}


#' @export
print.sj_revar <- function(x, ...) {
  # get parameters
  xn <- names(x)
  tau.00 <- x[string_ends_with("tau.00", xn)]
  tau.01 <- x[string_ends_with("tau.01", xn)]
  tau.11 <- x[string_ends_with("tau.11", xn)]
  rho.01 <- x[string_ends_with("rho.01", xn)]
  sigma_2 <- x[string_ends_with("sigma_2", xn)]

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
#' @export
print.sj_pca_rotate <- function(x, cutoff = .1, ...) {

  xs <- attr(x, "variance", exact = TRUE)

  rn <- rownames(x)

  x <- x %>%
    round(4) %>%
    purrr::map_df(~ dplyr::case_when(
      abs(.x) < cutoff ~ "",
      TRUE ~ as.character(.x)
    )) %>%
    as.data.frame() %>%
    add_cols(variable = rn, .after = -1)

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
print.sj_pca <- function(x, ...) {
  x <- as.data.frame(round(x, 4))
  rownames(x) <- c("Standard deviation", "Eigenvalue", "Proportion variance", "Cumulative variance")

  print(x, ...)
}


#' @importFrom crayon blue
#' @importFrom purrr map_if
#' @export
print.sj_rope <- function(x, digits = 1, ...) {
  cat(crayon::blue("\n# Proportions of samples inside and outside the ROPE\n\n"))

  # nicer column names for output
  colnames(x) <- c("term", "inside", "outside")

  # left-justify term column
  x$term <- format(x$term, justify = "left")

  x <- x %>%
    purrr::map_if(is.numeric, ~ sprintf("%.*f%%", digits, .x)) %>%
    as.data.frame()

  colnames(x)[1] <- ""
  print(x, ..., row.names = FALSE)
}


#' @importFrom crayon blue
#' @export
print.sj_hdi <- function(x, digits = 2, ...) {
  cat(crayon::blue("\n# Highest Density Interval\n\n"))

  dat <- get_hdi_data(x, digits)
  colnames(dat)[1] <- ""

  print(as.data.frame(dat), ..., row.names = FALSE)
}


#' @importFrom crayon blue cyan
#' @export
print.sj_equi_test <- function(x, ...) {
  cat(crayon::blue("\n# Test for Practical Equivalence of Model Predictors\n\n"))
  cat(crayon::cyan(sprintf(
    "  Effect Size: %.2f\n         ROPE: [%.2f %.2f]\n",
    attr(x, "eff_size", exact = TRUE),
    attr(x, "rope", exact = TRUE)[1],
    attr(x, "rope", exact = TRUE)[2]
  )))

  if (!is.null(attr(x, "nsamples", exact = TRUE))) {
    cat(crayon::cyan(sprintf(
      "      Samples: %i\n",
      attr(x, "nsamples", exact = TRUE)
    )))
  }

  cat("\n")

  dat <- get_hdi_data(x, digits = 2)
  dat[["inside.rope"]] <- sprintf("%.2f", dat[["inside.rope"]])
  colnames(dat) <- c("", "H0", "%inROPE", "HDI(95%)")

  print(as.data.frame(dat), ..., row.names = FALSE)

  if (isTRUE(attr(x, "critical"))) {
    message("\n(*) the number of effective samples may be insufficient for some parameters")
  }
}


#' @importFrom crayon blue cyan
#' @export
print.sj_mediation <- function(x, digits = 2, ...) {
  cat(crayon::blue("\n# Causal Mediation Analysis for Stan Model\n\n"))
  cat(crayon::cyan(sprintf(
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

  cat(crayon::red(
    sprintf(
      "\nProportion mediated: %s%% [%s%% %s%%]\n",
      prop.med[1], prop.med[2], prop.med[3]))
  )

  if (prop.med[1] < 0)
    message("\nDirect and indirect effects have opposite directions. The proportion mediated is not meaningful.")
}


#' @importFrom purrr map_at map_df
#' @importFrom dplyr bind_cols select
get_hdi_data <- function(x, digits) {
  cn <- colnames(x)
  prob <- attr(x, "prob", exact = TRUE)

  hdi.cols <- string_starts_with(pattern = "hdi.", x = cn)

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
        interv <- round(100 * prob[1])
      else
        interv <- 90
    } else
      interv <- round(100 * as.numeric(substr(cn[i], ci_pos[i] + 1, nchar(cn[i]))))

    colnames(tmp) <- sprintf("HDI(%d%%)", interv)
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
print.sj_error_rate <- function(x, ...) {
  cat(crayon::blue("\n# Error Rate of Logistic Regression Model\n\n"))
  cat(sprintf("  Full model: %.2f%%\n", 100 * x$error.model))
  cat(sprintf("  Null model: %.2f%%\n", 100 * x$error.null))

  cat(crayon::blue("\n# Likelihood-Ratio-Test\n\n"))

  v1 <- sprintf("%.3f", x$lrt.chisq)
  v2 <- sprintf("%.3f", x$lrt.p)

  space <- max(nchar(c(v1, v2)))

  cat(sprintf("  Chi-squared: %*s\n", space, v1))
  cat(sprintf("      p-value: %*s\n\n", space, v2))
}


#' @importFrom rlang .data
#' @export
print.sj_binres <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` required.", call. = F)
  }

  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package `scales` required.", call. = F)
  }

  x$se.lo <- -x$se
  if (length(unique(x$group)) > 1)
    ltitle <- "Within error bounds"
  else
    ltitle <- NULL

  term <- attr(x, "term", exact = TRUE)
  term.label <- attr(x, "term.label", exact = TRUE)

  if (is.null(term))
    xtitle <- sprintf("Estimated Probability of %s", attr(x, "resp_var", exact = TRUE))
  else
    xtitle = term.label

  p <- ggplot2::ggplot(data = x, ggplot2::aes_string(x = "xbar")) +
    ggplot2::geom_abline(slope = 0, intercept = 0, colour = "grey80")

  if (!is.null(term)) {
    p <- p +
      ggplot2::stat_smooth(
        ggplot2::aes_string(y = "ybar"),
        method = "loess",
        se = FALSE,
        colour = "#00b159",
        size = .5
      )
  }

  p <- p +
    ggplot2::geom_ribbon(ggplot2::aes_string(ymin = -Inf, ymax = "se.lo"), alpha = .1 , fill = "grey70") +
    ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "se", ymax = Inf), alpha = .1 , fill = "grey70") +
    ggplot2::geom_line(ggplot2::aes_string(y = "se"), colour = "grey70") +
    ggplot2::geom_line(ggplot2::aes_string(y = "se.lo"), colour = "grey70") +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("#d11141", "#00aedb")) +
    ggplot2::labs(
      y = "Average residual",
      x = xtitle,
      colour = ltitle
    )

  if (is.null(term)) {
    p <- p + ggplot2::scale_x_continuous(labels = scales::percent)
  }

  if (is.null(ltitle)) {
    p <- p + ggplot2::geom_point(ggplot2::aes_string(y = "ybar"))
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes_string(y = "ybar", colour = "group"))
  }

  suppressWarnings(graphics::plot(p))
}


#' @export
print.sj_chi2gof <- function(x, ...) {
  cat(crayon::blue("\n# Chi-squared Goodness-of-Fit Test\n\n"))

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


#' @importFrom crayon blue
#' @export
print.sj_hoslem <- function(x, ...) {
  cat(crayon::blue("\n# Hosmer-Lemeshow Goodness-of-Fit Test\n\n"))

  v1 <- sprintf("%.3f", x$chisq)
  v2 <- sprintf("%i    ", x$df)
  v3 <- sprintf("%.3f", x$p.value)

  space <- max(nchar(c(v1, v2, v3)))

  cat(sprintf("  Chi-squared: %*s\n", space, v1))
  cat(sprintf("           df: %*s\n", space, v2))
  cat(sprintf("      p-value: %*s\n\n", space, v3))

  if (x$p.value >= 0.05)
    message("Summary: model seems to fit well.")
  else
    message("Summary: model does not fit well.")
}


#' @importFrom crayon blue
#' @export
print.sj_check_assump <- function(x, ...) {
  cat(crayon::blue("\n# Checking Model-Assumptions\n\n"))
  cat(sprintf("  Model: %s", attr(x, "formula", exact = TRUE)))

  cat(crayon::red("\n\n                          violated    statistic\n"))

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


#' @importFrom crayon blue
#' @export
print.sj_item_diff <- function(x, ...) {
  cat(crayon::blue("\n# Item Difficulty\n\n"))

  items <- attr(x, "items", exact = TRUE)
  ideal <- attr(x, "ideal.difficulty", exact = TRUE)
  spaces <- max(nchar(items))

  cat(crayon::red(sprintf("  %*s  ideal\n", spaces + 10, "difficulty")))

  for (i in 1:length(items))
    cat(sprintf("  %*s      %.2f   %.2f\n", spaces, items[i], x[i], ideal[i]))
}


#' @importFrom crayon blue cyan
#' @export
print.sj_ttest <- function(x, ...) {
  cat(crayon::blue(sprintf("\n%s (%s)\n", x$method, x$alternative)))

  group <- attr(x, "group.name", exact = TRUE)
  xn <- attr(x, "x.name", exact = TRUE)
  yn <- attr(x, "y.name", exact = TRUE)

  if (!is.null(group))
    verbs <- c("of", "by")
  else
    verbs <- c("between", "and")

  st <- sprintf("# t=%.2f  df=%i  p-value=%.3f\n\n", x$statistic, as.integer(x$df), x$p.value)

  if (!is.null(yn)) {
    cat(crayon::cyan(sprintf("\n# comparison %s %s %s %s\n", verbs[1], xn, verbs[2], yn)))
  }

  cat(crayon::cyan(st))


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


#' @importFrom crayon blue cyan
#' @export
print.sj_wmwu <- function(x, ...) {
  cat(crayon::blue(sprintf("\n%s (%s)\n", x$method, x$alternative)))

  group <- attr(x, "group.name", exact = TRUE)
  xn <- attr(x, "x.name", exact = TRUE)

  cat(crayon::cyan(sprintf("\n# comparison of %s by %s\n", xn, group)))
  cat(crayon::cyan(sprintf("# Chisq=%.2f  df=%i  p-value=%.3f\n\n", x$statistic, as.integer(x$parameter), x$p.value)))
  cat(sprintf("  difference in mean rank score: %.3f\n\n", x$estimate))
}


#' @importFrom sjmisc round_num
#' @export
print.sj_anova_stat <- function(x, digits = 3, ...) {
  print.data.frame(sjmisc::round_num(x, digits), ..., row.names = TRUE)
}
