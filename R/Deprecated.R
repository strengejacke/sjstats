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


#' @importFrom parameters p_value standard_error
#' @rdname overdisp
#' @export
p_value <- function(x, ...) {
  # .Deprecated("parameters::p_value()")
  message("'p_value()' is deprecated. Use 'parameters::p_value()' instead.")
  pv <- parameters::p_value(x)
  se <- parameters::standard_error(x)
  out <- merge(pv, se, by = "Parameter")
  colnames(out) <- c("term", "p.value", "std.error")
  out
}


#' @importFrom parameters standard_error
#' @rdname overdisp
#' @export
se <- function(x, ...) {
  .Deprecated("parameters::standard_error()")
  parameters::standard_error(x)
}
