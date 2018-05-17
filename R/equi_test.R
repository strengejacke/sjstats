#' @rdname hdi
#' @export
equi_test <- function(x, rope, eff_size) {
  UseMethod("equi_test")
}


#' @export
equi_test.default <- function(x, rope, eff_size) {
  equi_test_worker(x = x, rope = rope, eff_size = eff_size, fm = NULL)
}


#' @export
equi_test.stanreg <- function(x, rope, eff_size) {
  dat <- equi_test_worker(x = x, rope = rope, eff_size = eff_size, fm = model_family(x))

  # check if we need to remove random or fixed effects
  dat <- remove_effects_from_stan(dat, type = "fixed", is.brms = FALSE)

  class(dat) <- c("sj_equi_test", class(dat))
  dat
}


#' @export
equi_test.brmsfit <- function(x, rope, eff_size) {
  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  dat <- equi_test_worker(x = x, rope = rope, eff_size = eff_size, fm = model_family(x))

  # check if we need to remove random or fixed effects
  dat <- remove_effects_from_stan(dat, type = "fixed", is.brms = TRUE)

  class(dat) <- c("sj_equi_test", class(dat))
  dat
}


#' @export
equi_test.stanfit <- function(x, rope, eff_size) {
  dat <- equi_test_worker(x = x, rope = rope, eff_size = eff_size, fm = model_family(x))

  # check if we need to remove random or fixed effects
  dat <- remove_effects_from_stan(dat, type = "fixed", is.brms = FALSE)

  class(dat) <- c("sj_equi_test", class(dat))
  dat
}


#' @export
equi_test.data.frame <- function(x, rope, eff_size) {
  dat <- equi_test_worker(x = x, rope = rope, eff_size = eff_size, fm = NULL)
  class(dat) <- c("sj_equi_test", class(dat))
  dat
}


#' @importFrom purrr map_df
#' @importFrom sjmisc add_columns var_rename is_empty
#' @importFrom tibble add_column
#' @importFrom dplyr case_when select pull
#' @importFrom stats sd
equi_test_worker <- function(x, rope, eff_size, fm) {

  if (fm$is_multivariate)
    stop("Multivariate response models not supported yet.", call. = F)


  dat <- as.data.frame(x)
  if (missing(eff_size)) eff_size <- .1

  if (!is.null(fm)) {
    if (missing(rope)) {
      if (fm$is_linear)
        eff_range <- stats::sd(resp_val(x)) * eff_size
      else
        eff_range <- stats::sd(dat[[1]]) * eff_size

      rope <- c(-1, 1) * eff_range
    }
  }

  .hdi <- hdi(x = x, prob = .95, trans = NULL, type = "fixed")
  .rope <- rope(x = x, rope = rope, trans = NULL, type = "fixed")
  .neff <- nrow(dat)

  result <- dplyr::case_when(
    .hdi$hdi.low > rope[2] ~ "reject",
    .hdi$hdi.high < rope[1] ~ "reject",
    .hdi$hdi.low >= rope[1] & .hdi$hdi.high <= rope[2] ~ "accept",
    TRUE ~ "undecided"
  )

  # for convenience reasons, also add proportion of values outside rope
  dat <- .hdi %>%
    dplyr::select(-1) %>%
    sjmisc::add_columns(.rope) %>%
    dplyr::select(-3) %>%
    tibble::add_column(decision = result, .after = 1) %>%
    sjmisc::var_rename(rope = "inside.rope")

  # indicate parameters with critical number of effective samples

  critical <- NULL
  if (inherits(x, c("stanfit", "stanreg", "brmsfit"))) {
    critical <- which(n_eff(x, type = "fixed") %>% dplyr::pull(-1) < .9)
    if (!sjmisc::is_empty(critical))
      dat$term[critical] <- sprintf("%s (*)", dat$term[critical])
  }

  attr(dat, "rope") <- rope
  attr(dat, "eff_size") <- eff_size
  attr(dat, "nsamples") <- .neff
  attr(dat, "critical") <- !sjmisc::is_empty(critical)

  dat
}
