# Help-functions


is_merMod <- function(fit) {
  return(inherits(fit, c("lmerMod", "glmerMod", "nlmerMod", "merModLmerTest")))
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

# return names of objects passed as ellipses argument
dot_names <- function(dots) unname(unlist(lapply(dots, as.character)))


#' @importFrom tidyr nest
#' @importFrom dplyr select_ filter
#' @importFrom stats complete.cases
get_grouped_data <- function(x) {
  # nest data frame
  grps <- tidyr::nest(x)

  # remove NA category
  cc <- grps %>%
    dplyr::select_("-data") %>%
    stats::complete.cases()
  # select only complete cases
  grps <- grps %>% dplyr::filter(cc)

  # arrange data
  if (length(attr(x, "vars", exact = T)) == 1)
    reihe <- order(grps[[1]])
  else
    reihe <- order(grps[[1]], grps[[2]])
  grps <- grps[reihe, ]

  grps
}
