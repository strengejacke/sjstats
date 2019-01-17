#' @rdname pred_vars
#' @importFrom purrr map_chr
#' @importFrom lme4 findbars
#' @importFrom stats formula
#' @importFrom sjmisc trim
#' @export
re_grp_var <- function(x) {
  tryCatch({
    if (inherits(x, "brmsfit"))
      f <- stats::formula(x)[[1]]
    if (inherits(x, "MixMod"))
      return(x$id_name)
    else
      f <- stats::formula(x)

    re <- purrr::map_chr(lme4::findbars(f), deparse)
    sjmisc::trim(substring(re, regexpr(pattern = "\\|", re) + 1))
  },
  error = function(x) { NULL }
  )
}


#' @rdname pred_vars
#' @export
grp_var <- function(x) {
  re_grp_var(x)
}
