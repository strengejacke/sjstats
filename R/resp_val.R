#' @rdname pred_vars
#' @importFrom insight get_response
#' @export
resp_val <- function(x) {
  .Deprecated("insight::get_response()")
  insight::get_response(x)
}
