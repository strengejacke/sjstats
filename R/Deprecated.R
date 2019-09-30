#' @title Deprecated functions
#' @name overdisp
#' @description A list of deprecated functions.
#'
#' @param x An object.
#' @param ... Currently not used.
#'
#' @return Nothing.
#'
#' @importFrom performance check_overdispersion
#' @export
overdisp <- function(x, ...) {
  .Deprecated("performance::check_overdispersion()")
  performance::check_overdispersion(x)
}



#' @importFrom performance check_zeroinflation
#' @rdname overdisp
#' @export
zero_count <- function(x, ...) {
  .Deprecated("performance::check_zeroinflation()")
  performance::check_zeroinflation(x)
}


#' @importFrom parameters principal_components
#' @rdname overdisp
#' @export
pca <- function(x, ...) {
  .Deprecated("parameters::principal_components()")
  parameters::principal_components(x)
}


#' @importFrom parameters principal_components
#' @rdname overdisp
#' @export
pca_rotate <- function(x, ...) {
  .Deprecated("parameters::principal_components()")
  parameters::principal_components(x, rotation = "varimax")
}


#' @importFrom performance r2
#' @rdname overdisp
#' @export
r2 <- function(x) {
  .Deprecated("performance::r2()")
  performance::r2(x)
}


#' @importFrom performance icc
#' @rdname overdisp
#' @export
icc <- function(x) {
  .Deprecated("performance::icc()")
  performance::icc(x)
}


#' @importFrom parameters p_value
#' @rdname overdisp
#' @export
p_value <- function(x, ...) {
  .Deprecated("parameters::p_value()")
  parameters::p_value(x)
}


#' @importFrom parameters se
#' @rdname overdisp
#' @export
se <- function(x, ...) {
  .Deprecated("parameters::se()")
  parameters::se(x)
}


#' @importFrom parameters parameters_standardize
#' @rdname overdisp
#' @export
std_beta <- function(x, ...) {
  .Deprecated("parameters::parameters_standardize()")
  parameters::parameters_standardize(x, method = "smart")
}
