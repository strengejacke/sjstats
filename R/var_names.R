#' @rdname pred_vars
#' @importFrom purrr map_chr
#' @export
var_names <- function(x) {
  if (is.character(x))
    get_vn_helper(x)
  else
    colnames(model_frame(x))
}


#' @importFrom sjmisc is_empty trim
#' @importFrom purrr map_chr
get_vn_helper <- function(x) {

  # return if x is empty
  if (sjmisc::is_empty(x)) return("")

  # for gam-smoothers/loess, remove s()- and lo()-function in column name
  # for survival, remove strata(), and so on...
  pattern <- c(
    "as.factor", "factor", "offset", "log-log", "log", "lag", "diff", "lo", "bs", "ns",
    "t2", "te", "ti", "tt", "mi", "mo", "gp", "pspline", "poly", "strata", "scale",
    "interaction", "s", "I"
  )

  # do we have a "log()" pattern here? if yes, get capture region
  # which matches the "cleaned" variable name
  purrr::map_chr(1:length(x), function(i) {
    for (j in 1:length(pattern)) {
      if (pattern[j] == "offset") {
        x[i] <- sjmisc::trim(unique(sub("^offset\\(([^-+ )]*).*", "\\1", x[i])))
      } else if (pattern[j] == "I") {
        x[i] <- sjmisc::trim(unique(sub("I\\((\\w*).*", "\\1", x[i])))
      } else if (pattern[j] == "log-log") {
        x[i] <- sjmisc::trim(unique(sub("^log\\(log\\(([^,)]*)).*", "\\1", x[i])))
      } else {
        p <- paste0("^", pattern[j], "\\(([^,)]*).*")
        x[i] <- unique(sub(p, "\\1", x[i]))
      }
    }
    # for coxme-models, remove random-effect things...
    sjmisc::trim(sub("^(.*)\\|(.*)", "\\2", x[i]))
    # x[i]
  })
}
