#' @title Deprecated functions
#' @name r2
#' @description A list of deprecated functions.
#'
#' @param x An object.
#' @param ... Currently not used.
#'
#' @return Nothing.
#'
#' @importFrom performance r2
#' @export
r2 <- function(x) {
  .Deprecated("performance::r2()")
  performance::r2(x)
}


#' @importFrom performance icc
#' @rdname r2
#' @export
icc <- function(x) {
  .Deprecated("performance::icc()")
  performance::icc(x)
}


#' @importFrom parameters p_value
#' @rdname r2
#' @export
p_value <- function(x, ...) {
  .Deprecated("parameters::p_value()")
  parameters::p_value(x)
}


#' @importFrom parameters standard_error
#' @rdname r2
#' @export
se <- function(x, ...) {
  .Deprecated("parameters::standard_error()")
  parameters::standard_error(x)
}


#' @importFrom effectsize cohens_f
#' @rdname r2
#' @export
cohens_f <- function(x, ...) {
  .Deprecated("effectsize::cohens_f()")
  effectsize::cohens_f(x)
}






#' @importFrom effectsize eta_squared
#' @rdname r2
#' @export
eta_sq <- function(x, ...) {
  .Deprecated("effectsize::eta_squared()")
  effectsize::eta_squared(x)
}


#' @importFrom effectsize epsilon_squared
#' @rdname r2
#' @export
epsilon_sq <- function(x, ...) {
  .Deprecated("effectsize::epsilon_squared()")
  effectsize::epsilon_squared(x)
}


#' @importFrom effectsize omega_squared
#' @rdname r2
#' @export
omega_sq <- function(x, ...) {
  .Deprecated("effectsize::omega_sqared()")
  effectsize::omega_squared(x)
}


#' @importFrom parameters rescale_weights
#' @rdname r2
#' @export
scale_weights <- function(x, ...) {
  .Deprecated("parameters::rescale_weights()")
  parameters::rescale_weights(x, ...)
}


#' @importFrom parameters model_parameters
#' @rdname r2
#' @export
tidy_stan <- function(x, ...) {
  .Deprecated("parameters::model_parameters()")
  parameters::model_parameters(x, ...)
}
