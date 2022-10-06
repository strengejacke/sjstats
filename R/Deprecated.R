#' @title Deprecated functions
#' @name r2
#' @description A list of deprecated functions.
#'
#' @param x An object.
#' @param ... Currently not used.
#'
#' @return Nothing.
#'
#' @export
r2 <- function(x) {
  .Defunct("performance::r2()")
  performance::r2(x)
}


#' @rdname r2
#' @export
cohens_f <- function(x, ...) {
  .Defunct("effectsize::cohens_f()")
  effectsize::cohens_f(x)
}


#' @rdname r2
#' @export
eta_sq <- function(x, ...) {
  .Defunct("effectsize::eta_squared()")
  effectsize::eta_squared(x)
}


#' @rdname r2
#' @export
epsilon_sq <- function(x, ...) {
  .Defunct("effectsize::epsilon_squared()")
  effectsize::epsilon_squared(x)
}


#' @rdname r2
#' @export
omega_sq <- function(x, ...) {
  .Defunct("effectsize::omega_sqared()")
  effectsize::omega_squared(x)
}


#' @rdname r2
#' @export
scale_weights <- function(x, ...) {
  .Defunct("datawizard::rescale_weights()")
  datawizard::rescale_weights(x, ...)
}


#' @rdname r2
#' @export
robust <- function(x, ...) {
  .Defunct("parameters::standard_error()")
  parameters::standard_error(x, ...)
}


#' @rdname r2
#' @export
icc <- function(x) {
  .Defunct("performance::icc()")
  performance::icc(x)
}

#' @rdname r2
#' @export
p_value <- function(x, ...) {
  .Defunct("parameters::p_value()")
  parameters::p_value(x)
}


#' @rdname r2
#' @export
se <- function(x, ...) {
  .Defunct("parameters::standard_error()")
  parameters::standard_error(x)
}
