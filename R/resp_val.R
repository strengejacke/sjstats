#' @rdname pred_vars
#' @importFrom insight get_response
#' @export
resp_val <- function(x) {
  insight::get_response(x)
}
