#' @importFrom purrr map_chr
#' @importFrom stats formula
#' @importFrom sjmisc is_empty trim
#' @rdname pred_vars
#' @export
resp_var <- function(x, combine = TRUE) {
  if (inherits(x, "brmsfit")) {
    if (is.null(stats::formula(x)$responses)) {
      rv <- deparse(stats::formula(x)$formula[[2L]])
      # check for brms Additional Response Information
      if (!sjmisc::is_empty(string_contains("|", rv))) {
        r1 <- sjmisc::trim(sub("(.*)\\|(.*)", "\\1", rv))
        r2 <- sjmisc::trim(sub("(.*)\\|(.*)\\(([^,)]*).*", "\\3", rv))
        rv <- c(r1, r2)
      }
    } else {
      rv <- purrr::map_chr(
        x$formula$forms,
        ~ stats::formula(.x)[[2L]] %>% all.vars()
      )
      # stats::formula(x)$responses
    }
  } else if (inherits(x, "stanmvreg")) {
    rv <- purrr::map_chr(stats::formula(x), ~ deparse(.x[[2L]]))
  } else if (inherits(x, "felm")) {
    rv <- x$lhs
  } else if (inherits(x, "clm2")) {
    rv <- all.vars(attr(x$location, "terms", exact = TRUE)[[2L]])
  } else if (inherits(x, "aovlist")) {
    rv <- all.vars(attr(x, "terms")[[2L]])
  } else if (inherits(x, "gam") && is.list(stats::formula(x))) {
    rv <- deparse(stats::formula(x)[[1]][[2L]])
  } else
    rv <- deparse(stats::formula(x)[[2L]])

  if (!combine && grepl("cbind\\((.*)\\)", rv) && !inherits(x, "brmsfit")) {
    rv <- sub("cbind\\(([^,].*)([\\)].*)" ,"\\1", rv) %>%
      strsplit(split  = ",", fixed = TRUE) %>%
      unlist() %>%
      sjmisc::trim()

    if (!sjmisc::is_empty(string_contains("-", rv[2])))
      rv[2] <- sjmisc::trim(sub("(.*)(\\-)(.*)", "\\1", rv[2]))
  }

  rv
}

