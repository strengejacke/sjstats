# Help-functions


is_merMod <- function(fit) {
  return(any(class(fit) %in% c("lmerMod", "glmerMod", "nlmerMod", "merModLmerTest")))
}

#' @importFrom sjmisc str_contains
get_glm_family <- function(fit) {
  c.f <- class(fit)
  # ------------------------
  # do we have glm? if so, get link family. make exceptions
  # for specific models that don't have family function
  # ------------------------
  if (any(c.f %in% c("lme", "plm"))) {
    fitfam <- ""
    logit_link <- FALSE
  } else {
    fitfam <- stats::family(fit)$family
    logit_link <- stats::family(fit)$link == "logit"
  }
  # --------------------------------------------------------
  # create logical for family
  # --------------------------------------------------------
  binom_fam <- fitfam %in% c("binomial", "quasibinomial")
  poisson_fam <- fitfam %in% c("poisson", "quasipoisson") ||
    sjmisc::str_contains(fitfam, "negative binomial", ignore.case = T)
  return(list(is_bin = binom_fam, is_pois = poisson_fam, is_logit = logit_link))
}
