#' @rdname pred_vars
#' @importFrom insight find_response
#' @export
resp_var <- function(x, combine = TRUE) {
  insight::find_response(x = x, combine = combine)
}
