# Help-functions

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' @importFrom sjmisc typical_value
#' @export
sjmisc::typical_value


data_frame <- function(...) {
  x <- data.frame(..., stringsAsFactors = FALSE)
  rownames(x) <- NULL
  x
}


is_merMod <- function(fit) {
  inherits(fit, c("lmerMod", "glmerMod", "nlmerMod", "merModLmerTest"))
}


is_stan_model <- function(fit) {
  inherits(fit, c("stanreg", "stanfit", "brmsfit"))
}


#' @importFrom sjmisc str_contains
get_glm_family <- function(fit) {
  c.f <- class(fit)

  # do we have glm? if so, get link family. make exceptions
  # for specific models that don't have family function
  if (any(c.f %in% c("lme", "plm"))) {
    fitfam <- ""
    logit_link <- FALSE
  } else {
    fitfam <- stats::family(fit)$family
    logit_link <- stats::family(fit)$link == "logit"
  }

  # create logical for family
  binom_fam <- fitfam %in% c("binomial", "quasibinomial")
  poisson_fam <- fitfam %in% c("poisson", "quasipoisson") ||
    sjmisc::str_contains(fitfam, "negative binomial", ignore.case = T)

  list(is_bin = binom_fam, is_pois = poisson_fam, is_logit = logit_link)
}

# return names of objects passed as ellipses argument
dot_names <- function(dots) unname(unlist(lapply(dots, as.character)))


#' @importFrom tidyr nest
#' @importFrom dplyr select filter group_vars
#' @importFrom stats complete.cases
#' @importFrom rlang .data
get_grouped_data <- function(x) {
  # nest data frame
  grps <- tidyr::nest(x)

  # remove NA category
  cc <- grps %>%
    dplyr::select(-.data$data) %>%
    stats::complete.cases()

  # select only complete cases
  grps <- grps %>% dplyr::filter(!! cc)

  # arrange data
  if (length(dplyr::group_vars(x)) == 1)
    reihe <- order(grps[[1]])
  else
    reihe <- order(grps[[1]], grps[[2]])
  grps <- grps[reihe, ]

  grps
}
