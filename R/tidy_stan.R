#' @title Tidy summary output for stan models
#' @name tidy_stan
#'
#' @description Returns a tidy summary output for stan models.
#'
#' @param x A \code{stanreg}, \code{stanfit} or \code{brmsfit} object.
#' @param trans Name of a function or character vector naming a function, used
#'        to apply transformations on the estimate and HDI-values. The
#'        values for standard errors are \emph{not} transformed!
#' @param digits Amount of digits to round numerical values in the output.
#'
#' @inheritParams hdi
#'
#' @return A tidy data frame, summarizing \code{x}, with consistent column names.
#'         To distinguish multiple HDI values, column names for the HDI get a suffix
#'         when \code{prob} has more than one element.
#'
#' @details The returned data frame has an additonal class-attribute,
#'    \code{tidy_stan}, to pass the result to its own \code{print()}-method.
#'    The \code{print()}-method create a cleaner output, especially for multilevel,
#'    zero-inflated or multivariate response models, where - for instance -
#'    the conditional part of a model is printed separately from the zero-inflated
#'    part, or random and fixed effects are printed separately.
#'    \cr \cr
#'    The returned data frame gives information on:
#'    \itemize{
#'      \item{The Bayesian point estimate (column \emph{estimate}, which is by
#'            default the posterior median; other statistics are also possible,
#'            see argument \code{typical}).}
#'      \item{
#'        The standard error (which is actually the \emph{median absolute deviation}).
#'      }
#'      \item{
#'        The HDI (see \code{\link{hdi}}). Computation for HDI is based on the
#'        code from Kruschke 2015, pp. 727f.
#'      }
#'      \item{
#'        The ratio of effective numbers of samples, \emph{neff_ratio}, (i.e.
#'        effective number of samples divided by total number of samples).
#'        This ratio ranges from 0 to 1, and should be close to 1. The closer
#'        this ratio comes to zero means that the chains may be inefficient,
#'        but possibly still okay.
#'      }
#'      \item{
#'        The Rhat statistics. When Rhat is above 1, it usually indicates that
#'        the chain has not yet converged, indicating that the drawn samples
#'        might not be trustworthy. Drawing more iteration may solve this issue.
#'      }
#'      \item{
#'        The Monte Carlo standard error (see \code{\link{mcse}}). It is defined
#'        as standard deviation of the chains divided by their effective sample
#'        size and \dQuote{provides a quantitative suggestion of how big the
#'        estimation noise is} (\emph{Kruschke 2015, p.187}).
#'      }
#'    }
#'
#' @references Kruschke JK. \emph{Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan.} 2nd edition. Academic Press, 2015
#' \cr \cr
#' Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB. \emph{Bayesian data analysis.} 3rd ed. Boca Raton: Chapman & Hall/CRC, 2013
#' \cr \cr
#' Gelman A, Rubin DB. \emph{Inference from iterative simulation using multiple sequences.} Statistical Science 1992;7: 457–511
#' \cr \cr
#' McElreath R. \emph{Statistical Rethinking. A Bayesian Course with Examples in R and Stan.} Chapman and Hall, 2015
#'
#' @examples
#' \dontrun{
#' if (require("rstanarm")) {
#'   fit <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)
#'   tidy_stan(fit)
#'   tidy_stan(fit, prob = c(.89, .5))
#' }}
#'
#' @importFrom purrr map flatten_dbl map_dbl modify_if
#' @importFrom dplyr bind_cols select mutate slice inner_join n_distinct
#' @importFrom stats mad formula
#' @importFrom bayesplot rhat neff_ratio
#' @importFrom sjmisc is_empty trim
#' @export
tidy_stan <- function(x, prob = .89, typical = "median", trans = NULL, type = c("fixed", "random", "all"), digits = 2) {

  # only works for rstanarm- or brms-models
  if (!inherits(x, c("stanreg", "stanfit", "brmsfit")))
    stop("`x` needs to be a stanreg- or brmsfit-object.", call. = F)

  # check arguments
  type <- match.arg(type)

  # get data frame and family info
  mod.dat <- as.data.frame(x)
  faminfo <- model_family(x)

  # for brmsfit models, we need to remove some columns here to
  # match data rows later
  if (inherits(x, "brmsfit")) mod.dat <- brms_clean(mod.dat)

  # compute HDI
  out.hdi <- hdi(x, prob = prob, trans = trans, type = type)

  # get statistics
  nr <- bayesplot::neff_ratio(x)

  # we need names of elements, for correct removal

  if (inherits(x, "brmsfit")) {
    cnames <- make.names(names(nr))
    keep <- cnames %in% out.hdi$term
  } else {
    keep <- names(nr) %in% out.hdi$term
  }

  nr <- nr[keep]
  ratio <- data_frame(
    term = names(nr),
    ratio = nr
  )


  rh <- bayesplot::rhat(x)

  if (inherits(x, "brmsfit")) {
    cnames <- make.names(names(rh))
    keep <- cnames %in% out.hdi$term
  } else {
    keep <- names(rh) %in% out.hdi$term
  }

  rh <- rh[keep]
  rhat <- data_frame(
    term = names(rh),
    rhat = rh
  )

  if (inherits(x, "brmsfit")) {
    ratio$term <- make.names(ratio$term)
    rhat$term <- make.names(rhat$term)
  }

  se <- mcse(x, type = type)
  se <- se[se$term %in% out.hdi$term, ]

  est <- purrr::map_dbl(mod.dat, ~ sjstats::typical_value(.x, fun = typical))

  out <- data_frame(
    term = names(est),
    estimate = est,
    std.error = purrr::map_dbl(mod.dat, stats::mad)
  ) %>%
    dplyr::inner_join(
      out.hdi,
      by = "term"
    ) %>%
    dplyr::inner_join(
      ratio,
      by = "term"
    ) %>%
    dplyr::inner_join(
      rhat,
      by = "term"
    ) %>%
    dplyr::inner_join(
      se,
      by = "term"
    )


  # transform estimate, if requested

  if (!is.null(trans)) {
    trans <- match.fun(trans)
    out$estimate <- trans(out$estimate)
  }


  # check if we need to remove random or fixed effects
  out <- remove_effects_from_stan(out, type, is.brms = inherits(x, "brmsfit"))


  # create variable to index random effects

  if (type == "random" || type == "all") {

    out <- add_cols(out, random.effect = "", .after = -1)

    # find random intercepts

    ri <- grep("b\\[\\(Intercept\\) (.*)\\]", out$term)

    if (!sjmisc::is_empty(ri)) {
      out$random.effect[ri] <- "(Intercept)"
      out$term[ri] <- gsub("b\\[\\(Intercept\\) (.*)\\]", "\\1", out$term[ri])

      # check if we have multiple intercepts or nested groups
      # if yes, append group name to intercept label

      multi.grps <- gsub(pattern = ":([^:]+)", "\\2", out$term[ri])

      if (dplyr::n_distinct(multi.grps) > 1) {
        out$random.effect[ri] <- sprintf("(Intercept: %s)", multi.grps)
      }
    }


    ## TODO if random intercept variable has levels with "/" or ".", these are also dot-separated
    ## and will be confused with nested groups

    # find random intercepts

    ri1 <- grep("r_(.*)\\.(.*)\\.", out$term)
    ri2 <- which(gsub("r_(.*)\\.(.*)\\.", "\\2", out$term) == "Intercept")

    if (!sjmisc::is_empty(ri1)) {
      ri <- intersect(ri1, ri2)
      out$random.effect[ri] <- "(Intercept)"
      out$term[ri] <- gsub("r_(.*)\\.(.*)\\.", "\\1", out$term[ri])

      # check if we have multiple intercepts or nested groups
      # if yes, append group name to intercept label

      multi.grps <- gsub(pattern = "\\.([^\\.]+)$", "\\2", out$term[ri])

      if (dplyr::n_distinct(multi.grps) > 1) {
        out$random.effect[ri] <- sprintf("(Intercept: %s)", multi.grps)
      }
    }


    # find residual variance for random intercept

    rsig1 <- which(gsub("(Sigma)\\[(.*)\\,(.*)\\]", "\\1", out$term) == "Sigma")
    rsig2 <- which(gsub("(Sigma)\\[(.*)\\,(.*)\\]", "\\3", out$term) == "(Intercept)")

    if (!sjmisc::is_empty(rsig1) && !sjmisc::is_empty(rsig2)) {
      rs <- intersect(rsig1, rsig2)
      out$random.effect[rs] <- "(Intercept)"

      out$term[rs] <- gsub(
        pattern = ":(Intercept)",
        replacement = "",
        sprintf("sigma (%s)", gsub("(Sigma)\\[(.*)\\,(.*)\\]", "\\2", out$term)[rs]),
        fixed = TRUE
      )

      rs <- setdiff(rsig1, rsig2)

      if (!sjmisc::is_empty(rs)) {
        out$random.effect[rs] <- gsub("(Sigma)\\[(.*)\\,(.*)\\]", "\\3", out$term)[rs]
        out$term[rs] <- "sigma"
      }
    }


    ## TODO extract Sigma for stanmvreg random effects


    # find random slopes

    rs1 <- grep("b\\[(.*) (.*)\\]", out$term)
    rs2 <- which(gsub("b\\[(.*) (.*)\\]", "\\1", out$term) != "(Intercept)")

    if (!sjmisc::is_empty(rs1)) {
      rs <- intersect(rs1, rs2)
      rs.string <- gsub("b\\[(.*) (.*)\\]", "\\1", out$term[rs])
      out$random.effect[rs] <- rs.string
      out$term[rs] <- gsub("b\\[(.*) (.*)\\]", "\\2", out$term[rs])
    }


    # find random slopes

    rs1 <- grep("r_(.*)\\.(.*)\\.", out$term)
    rs2 <- which(gsub("r_(.*)\\.(.*)\\.", "\\2", out$term) != "Intercept")

    if (!sjmisc::is_empty(rs1)) {
      rs <- intersect(rs1, rs2)
      rs.string <- gsub("r_(.*)\\.(.*)\\.", "\\2", out$term[rs])
      out$random.effect[rs] <- rs.string
      out$term[rs] <- gsub("r_(.*)\\.(.*)\\.", "\\1", out$term[rs])
    }

    # did we really had random effects?

    if (obj_has_name(out, "random.effect") &&
        all(sjmisc::is_empty(out$random.effect, first.only = FALSE)))
      out <- dplyr::select(out, -.data$random.effect)
  }


  # categorical model?

  if (inherits(x, "brmsfit") && faminfo$is_categorical) {

    # terms of categorical models are prefixed with "mu"

    if (length(string_starts_with(pattern = "b_mu", x = out$term)) == nrow(out)) {
      out$term <- substr(out$term, 5, max(nchar(out$term)))
      # create "response-level" variable
      out <- add_cols(out, response.level = "", .after = -1)
      out$response.level <- gsub("(.*)\\_(.*)", "\\1", out$term)
      out$term <- gsub("(.*)\\_(.*)", "\\2", out$term)
    }
  }


  # multivariate-response model?

  if (inherits(x, "brmsfit") && faminfo$is_multivariate) {

    # get response variables

    responses <- stats::formula(x)$responses

    # also clean prepared data frame
    resp.sigma1 <- string_starts_with(pattern = "sigma_", x = out$term)
    resp.sigma2 <- string_starts_with(pattern = "b_sigma_", x = out$term)

    resp.sigma <- c(resp.sigma1, resp.sigma2)

    if (!sjmisc::is_empty(resp.sigma))
      out <- dplyr::slice(out, !! -resp.sigma)


    # create "response-level" variable

    out <- add_cols(out, response = "", .after = -1)

    # check if multivariate response model also has random effects
    # we need to clean names for the random effects as well here

    if (obj_has_name(out, "random.effect")) {
      re <- which(!sjmisc::is_empty(sjmisc::trim(out$random.effect), first.only = FALSE))
    } else {
      re <- NULL
    }

    # copy name of response into new character variable
    # and remove response name from term name

    for (i in responses) {
      m <- string_contains(pattern = i, x = out$term)
      out$response[intersect(which(out$response == ""), m)] <- i
      out$term <- gsub(sprintf("b_%s_", i), "", out$term, fixed = TRUE)
      out$term <- gsub(sprintf("s_%s_", i), "", out$term, fixed = TRUE)

      if (!sjmisc::is_empty(re)) {
        out$random.effect[re] <- gsub(sprintf("__%s", i), "", out$random.effect[re], fixed = TRUE)
        out$term[re] <- gsub(sprintf("__%s", i), "", out$term[re], fixed = TRUE)
      }
    }

  }


  if (inherits(x, "stanmvreg")) {

    # get response variables

    responses <- resp_var(x)
    resp.names <- names(responses)


    # create "response-level" variable

    out <- add_cols(out, response = "", .after = -1)


    # copy name of response into new character variable
    # and remove response name from term name

    for (i in 1:length(responses)) {
      pattern <- paste0(resp.names[i], "|")
      m <- string_starts_with(pattern = pattern, x = out$term)
      out$response[intersect(which(out$response == ""), m)] <- responses[i]
      out$term <- gsub(pattern, "", out$term, fixed = TRUE)
    }

  }


  class(out) <- c("tidy_stan", class(out))

  attr(out, "digits") <- digits
  attr(out, "model_name") <- deparse(substitute(x))
  attr(out, "prob") <- prob

  if (inherits(x, "brmsfit"))
    attr(out, "formula") <- as.character(stats::formula(x))[1]
  else
    attr(out, "formula") <- deparse(stats::formula(x))

  # round values
  purrr::modify_if(out, is.numeric, ~ round(.x, digits = digits))
}


#' @importFrom purrr map_dbl
#' @importFrom dplyr slice
#' @importFrom sjmisc is_empty
remove_effects_from_stan <- function(out, type, is.brms) {

  # brmsfit-objects also include sd and cor for mixed
  # effecs models, so remove these here

  if (is.brms) out <- brms_clean(out)


  # remove certain terms like log-posterior etc. from output

  keep <- seq_len(length(out$term))

  for (.x in c("^(?!lp__)", "^(?!log-posterior)", "^(?!mean_PPD)")) {
    keep <- intersect(keep, grep(.x, out$term, perl = TRUE))
  }

  out <- dplyr::slice(out, !! keep)


  # if user wants all terms, return data here

  if (type == "all") return(out)


  # check if we need to remove random or fixed effects
  # therefor, find random effect parts first

  alt.term.names <- make.names(out$term)

  re <- string_starts_with(pattern = "b[", x = out$term)
  re.s <- string_starts_with(pattern = "Sigma[", x = out$term)
  re.i1 <- intersect(
    string_starts_with(pattern = "r_", x = out$term),
    string_ends_with(pattern = ".", x = out$term)
  )
  re.i2 <- intersect(
    string_starts_with(pattern = "r_", x = alt.term.names),
    string_ends_with(pattern = ".", x = alt.term.names)
  )

  removers <- unique(c(re, re.s, re.i1, re.i2))

  if (!sjmisc::is_empty(removers)) {
    if (type == "fixed") {
      # remove all random effects
      out <- dplyr::slice(out, !! -removers)
    } else if (type == "random") {
      # remove all fixed effects
      out <- dplyr::slice(out, !! removers)
    }
  }


  out
}


#' @importFrom dplyr slice
#' @importFrom sjmisc is_empty
brms_clean <- function(out) {

  # brmsfit-objects also include sd and cor for mixed
  # effecs models, so remove these here

  if (obj_has_name(out, "term")) {
    re.sd <- string_starts_with(pattern = "sd_", x = out$term)
    re.cor <- string_starts_with(pattern = "cor_", x = out$term)
    re.s <- string_starts_with(pattern = "Sigma[", x = out$term)
    lp <- string_starts_with(pattern = "lp__", x = out$term)
    priors <- string_starts_with(pattern = "prior_", x = out$term)

    removers <- unique(c(re.sd, re.cor, re.s, lp, priors))

    if (!sjmisc::is_empty(removers))
      out <- dplyr::slice(out, !! -removers)
  }


  # we may have transformed data frame, where columns
  # need to be removed

  re.sd <- string_starts_with(pattern = "sd_", x = colnames(out))
  re.cor <- string_starts_with(pattern = "cor_", x = colnames(out))
  re.s <- string_starts_with(pattern = "Sigma[", x = colnames(out))
  lp <- string_starts_with(pattern = "lp__", x = colnames(out))
  priors <- string_starts_with(pattern = "prior_", x = colnames(out))

  removers <- unique(c(re.sd, re.cor, re.s, lp, priors))

  if (!sjmisc::is_empty(removers))
    out <- dplyr::select(out, !! -removers)

  out
}


#' @importFrom purrr map_dbl
n_of_samples <- function(x) {
  if (inherits(x, "brmsfit") && requireNamespace("brms", quietly = TRUE)) {
    brms::nsamples(x, incl_warmup = FALSE)
  } else {
    purrr::map_dbl(x$stanfit@stan_args, ~ .x$iter - .x$warmup) %>% sum()
  }
}


n_of_chains <- function(x) {
  if (inherits(x, "brmsfit"))
    length(x$fit@stan_args)
  else
    length(x$stanfit@stan_args)
}
