# Help-functions

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
    sjmisc::str_contains(fitfam, "negative binomial", ignore.case = TRUE)

  list(is_bin = binom_fam, is_pois = poisson_fam, is_logit = logit_link)
}

# return names of objects passed as ellipses argument
dot_names <- function(dots) unname(unlist(lapply(dots, as.character)))


#' @importFrom tidyr nest
#' @importFrom dplyr select filter group_vars
#' @importFrom stats complete.cases
#' @importFrom rlang .data
get_grouped_data <- function(x) {
  # retain observations that are complete wrt grouping vars, then nest
  grps <- x %>%
    dplyr::group_modify(~ dplyr::filter(.x, stats::complete.cases(.y))) %>%
    tidyr::nest()

  # arrange data
  if (length(dplyr::group_vars(x)) == 1)
    reihe <- order(grps[[1]])
  else
    reihe <- order(grps[[1]], grps[[2]])
  grps <- grps[reihe, ]

  grps
}


.compact_character <- function(x) {
  x[!sapply(x, function(i) is.null(i) || nchar(i) == 0 || is.na(i) || any(i == "NULL", na.rm = TRUE))]
}


.format_symbols <- function(x) {
  if (.unicode_symbols()) {
    x <- gsub("Delta", "\u0394", x, ignore.case = TRUE)
    x <- gsub("Phi", "\u03D5", x, ignore.case = TRUE)
    x <- gsub("Eta", "\u03B7", x, ignore.case = TRUE)
    x <- gsub("Epsilon", "\u03b5", x, ignore.case = TRUE)
    x <- gsub("Omega", "\u03b5", x, ignore.case = TRUE)
    x <- gsub("R2", "R\u00b2", x, ignore.case = TRUE)
    x <- gsub("Chi2", "\u03C7\u00b2", x, ignore.case = TRUE)
    x <- gsub("Chi-squared", "\u03C7\u00b2", x, ignore.case = TRUE)
    x <- gsub("Chi", "\u03C7", x, ignore.case = TRUE)
    x <- gsub("Sigma", "\u03C3", x, ignore.case = TRUE)
    x <- gsub("Rho", "\u03C1", x, ignore.case = TRUE)
    x <- gsub("Mu", "\u03BC", x, ignore.case = TRUE)
    x <- gsub("Theta", "\u03B8", x, ignore.case = TRUE)
  }
  x
}


.unicode_symbols <- function() {
  win_os <- tryCatch(
    {
      si <- Sys.info()
      if (is.null(si["sysname"])) {
        FALSE
      } else {
        si["sysname"] == "Windows" || startsWith(R.version$os, "mingw")
      }
    },
    error = function(e) {
      TRUE
    }
  )
  l10n_info()[["UTF-8"]] && ((win_os && getRversion() >= "4.2") || (!win_os && getRversion() >= "4.0"))
}
