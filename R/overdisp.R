#' @title Deprecated functions
#' @name overdisp
#' @description A list of deprecated functions.
#'
#' @param x An object.
#' @param ... Currently not used.
#'
#' @return Nothing.
#'
#' @importFrom performance check_overdispersion check_zeroinflation check_convergence check_singularity item_reliability item_split_half cronbachs_alpha item_difficulty item_intercor principal_components
#' @export
overdisp <- function(x, ...) {
  .Deprecated("performance::check_overdispersion()")
  performance::check_overdispersion(x)
}

#' @rdname overdisp
#' @export
zero_count <- function(x, ...) {
  .Deprecated("performance::check_zeroinflation()")
  performance::check_zeroinflation(x)
}

#' @rdname overdisp
#' @export
converge_ok <- function(x, ...) {
  .Deprecated("performance::check_convergence()")
  performance::check_convergence(x)
}

#' @rdname overdisp
#' @export
is_singular <- function(x, ...) {
  .Deprecated("performance::check_singularity()")
  performance::check_singularity(x)
}

#' @rdname overdisp
#' @export
reliab_test <- function(x, ...) {
  .Deprecated("performance::item_reliability()")
  performance::item_reliability(x)
}

#' @rdname overdisp
#' @export
split_half <- function(x, ...) {
  .Deprecated("performance::item_split_half()")
  performance::item_split_half(x)
}

#' @rdname overdisp
#' @export
cronb <- function(x, ...) {
  .Deprecated("performance::cronbachs_alpha()")
  performance::cronbachs_alpha(x)
}

#' @rdname overdisp
#' @export
difficulty <- function(x, ...) {
  .Deprecated("performance::item_difficulty()")
  performance::item_difficulty(x)
}

#' @rdname overdisp
#' @export
mic <- function(x, ...) {
  .Deprecated("performance::item_intercor()")
  performance::item_intercor(x)
}

#' @rdname overdisp
#' @export
pca <- function(x, ...) {
  .Deprecated("performance::principal_components()")
  performance::principal_components(x)
}

#' @rdname overdisp
#' @export
pca_rotate <- function(x, ...) {
  .Deprecated("performance::principal_components()")
  performance::principal_components(x, rotation = "varimax")
}


#' @rdname overdisp
#' @export
n_eff <- function(x, ...) {
  .Deprecated("bayestestR::effective_sample()")
  bayestestR::effective_sample(x, ...)
}


#' @rdname overdisp
#' @export
mcse <- function(x, ...) {
  .Deprecated("bayestestR::mcse()")
  bayestestR::mcse(x, ...)
}
