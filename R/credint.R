#' @export
#' @rdname hdi
cred_int <- function(x, ...) {
  UseMethod("cred_int")
}


#' @rdname hdi
#' @export
cred_int.stanreg <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  cred_int_worker(x = x, prob = prob, trans = trans, type = type)
}


#' @rdname hdi
#' @export
cred_int.brmsfit <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all"), ...) {
  # check arguments
  type <- match.arg(type)

  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  cred_int_worker(x = x, prob = prob, trans = trans, type = type)
}


#' @export
cred_int.stanfit <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  cred_int_worker(x = x, prob = prob, trans = trans, type = type)
}


#' @export
cred_int.data.frame <- function(x, prob = .9, trans = NULL, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  cred_int_worker(x = x, prob = prob, trans = trans, type = type)
}


#' @export
cred_int.default <- function(x, prob = .9, trans = NULL, ...) {
  cred_int_helper(x, prob, trans)
}


#' @importFrom purrr map_df map2_df
#' @importFrom sjmisc rotate_df
cred_int_worker <- function(x, prob, trans, type) {
  dat <- purrr::map(
    prob,
    function(i) {
      tmp <- as.data.frame(x)
      out <- purrr::map2_df(tmp, colnames(tmp), ~ cred_int_helper(x = .x, prob = i, trans = trans, cn = .y)) %>%
        sjmisc::rotate_df() %>%
        rownames_as_column()

      colnames(out) <- c("term", "ci.low", "ci.high")
      out

    }) %>%
    dplyr::bind_cols()

  dat <- dplyr::select(dat, 1, string_starts_with(pattern = "ci.", x = colnames(dat)))



  # for multiple HDIs, fix column names

  if (length(prob) > 1) {
    suffix <- prob %>%
      purrr::map(~ rep(.x, 2)) %>%
      purrr::flatten_dbl()

    colnames(dat)[2:ncol(dat)] <-
      sprintf(
        "%s_%s",
        rep(c("ci.low", "ci.high"), length(prob)),
        as.character(suffix)
      )
  }

  attr(dat, "prob") <- prob

  if (is_stan_model(x)) {
    # check if we need to remove random or fixed effects
    dat <- remove_effects_from_stan(dat, type, is.brms = inherits(x, "brmsfit"))
  }

  class(dat) <- c("sj_credint", class(dat))
  dat
}


# based on Kruschke 2015, pp727f
#' @importFrom stats quantile
cred_int_helper <- function(x, prob, trans, cn = NULL) {
  x <- stats::quantile(x, probs = c((1 - prob) / 2, (1 + prob) / 2))
  if (!is.null(trans) && (is.null(cn) || !grepl("^simo_mo", cn))) {
    trans <- match.fun(trans)
    x <- trans(x)
  }
  x
}
