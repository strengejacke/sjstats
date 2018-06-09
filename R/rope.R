#' @rdname hdi
#' @export
rope <- function(x, rope, ...) {
  UseMethod("rope")
}


#' @export
rope.default <- function(x, rope, trans = NULL, type = c("fixed", "random", "all"), ...) {
  rope_helper(x, rope, trans)
}


#' @rdname hdi
#' @export
rope.stanreg <- function(x, rope, trans = NULL, type = c("fixed", "random", "all"), ...) {
  type <- match.arg(type)
  rope_worker(x = x, rope = rope, trans = trans, type = type)
}


#' @rdname hdi
#' @export
rope.brmsfit <- function(x, rope, trans = NULL, type = c("fixed", "random", "all"), ...) {
  # check arguments
  type <- match.arg(type)

  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  rope_worker(x = x, rope = rope, trans = trans, type = type)
}


#' @export
rope.stanfit <- function(x, rope, trans = NULL, type = c("fixed", "random", "all"), ...) {
  # check arguments
  type <- match.arg(type)
  rope_worker(x = x, rope = rope, trans = trans, type = type)
}


#' @export
rope.data.frame <- function(x, rope, trans = NULL, type = c("fixed", "random", "all"), ...) {
  # check arguments
  type <- match.arg(type)
  rope_worker(x = x, rope = rope, trans = trans, type = type)
}


#' @importFrom purrr map_df
#' @importFrom sjmisc rotate_df
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr mutate
#' @importFrom rlang .data
rope_worker <- function(x, rope, trans, type) {
  # get posterior data
  dat <- x %>%
    as.data.frame() %>%
    purrr::map_df(~ rope_helper(.x, rope, trans)) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column()

  colnames(dat) <- c("term", "rope")

  # for convenience reasons, also add proportion of values outside rope
  dat <- dplyr::mutate(dat, outside.rope = 100 - .data$rope)

  if (is_stan_model(x)) {
    # check if we need to remove random or fixed effects
    dat <- remove_effects_from_stan(dat, type, is.brms = inherits(x, "brmsfit"))
  }

  class(dat) <- c("sj_rope", class(dat))
  dat
}


#' @importFrom dplyr between
rope_helper <- function(x, rope, trans) {
  # stop if argument is not correct
  if (length(rope) != 2)
    stop("Argument `rope` needs to be a vector of length two.", call. = F)

  # switch values, if lower bound is larger than upper bound
  if (rope[1] > rope[2]) {
    tmp <- rope[2]
    rope[2] <- rope[1]
    rope[1] <- tmp
  }

  # check if we have correct function
  if (!is.null(trans)) {
    trans <- match.fun(trans)
    x <- trans(x)
  }

  # sort values, to compute rope
  x <- sort(x)
  r <- dplyr::between(x, rope[1], rope[2])

  # compute proportion of values within boundaries
  round(100 * sum(r) / length(x), 3)
}

