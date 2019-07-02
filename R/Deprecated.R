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


#' @importFrom performance check_convergence
#' @rdname overdisp
#' @export
converge_ok <- function(x, ...) {
  .Deprecated("performance::check_convergence()")
  performance::check_convergence(x)
}


#' @importFrom performance check_singularity
#' @rdname overdisp
#' @export
is_singular <- function(x, ...) {
  .Deprecated("performance::check_singularity()")
  performance::check_singularity(x)
}


#' @importFrom performance item_reliability
#' @rdname overdisp
#' @export
reliab_test <- function(x, ...) {
  .Deprecated("performance::item_reliability()")
  performance::item_reliability(x)
}


#' @importFrom performance item_split_half
#' @rdname overdisp
#' @export
split_half <- function(x, ...) {
  .Deprecated("performance::item_split_half()")
  performance::item_split_half(x)
}


#' @importFrom performance cronbachs_alpha
#' @rdname overdisp
#' @export
cronb <- function(x, ...) {
  .Deprecated("performance::cronbachs_alpha()")
  performance::cronbachs_alpha(x)
}


#' @importFrom performance item_difficulty
#' @rdname overdisp
#' @export
difficulty <- function(x, ...) {
  .Deprecated("performance::item_difficulty()")
  performance::item_difficulty(x)
}


#' @importFrom performance item_intercor
#' @rdname overdisp
#' @export
mic <- function(x, ...) {
  .Deprecated("performance::item_intercor()")
  performance::item_intercor(x)
}


#' @importFrom performance principal_components
#' @rdname overdisp
#' @export
pca <- function(x, ...) {
  .Deprecated("performance::principal_components()")
  performance::principal_components(x)
}


#' @importFrom performance principal_components
#' @rdname overdisp
#' @export
pca_rotate <- function(x, ...) {
  .Deprecated("performance::principal_components()")
  performance::principal_components(x, rotation = "varimax")
}


#' @importFrom bayestestR effective_sample
#' @rdname overdisp
#' @export
n_eff <- function(x, ...) {
  .Deprecated("bayestestR::effective_sample()")
  bayestestR::effective_sample(x, ...)
}


#' @importFrom bayestestR mcse
#' @rdname overdisp
#' @export
mcse <- function(x, ...) {
  .Deprecated("bayestestR::mcse()")
  bayestestR::mcse(x, ...)
}



#' @importFrom insight find_predictors
#' @rdname overdisp
#' @export
pred_vars <- function(x, ...) {
  .Deprecated("insight::find_predictors()")
  insight::find_predictors(x, ...)
}



#' @importFrom insight model_info
#' @rdname overdisp
#' @export
model_family <- function(x, ...) {
  .Deprecated("insight::model_info()")
  insight::model_info(x, ...)
}


#' @importFrom insight get_data
#' @rdname overdisp
#' @export
model_frame <- function(x, ...) {
  .Deprecated("insight::get_data()")
  insight::get_data(x, ...)
}



#' @rdname overdisp
#' @importFrom insight get_response
#' @export
resp_val <- function(x, ...) {
  .Deprecated("insight::get_response()")
  insight::get_response(x, ...)
}



#' @rdname overdisp
#' @importFrom insight find_response
#' @export
resp_var <- function(x, ...) {
  .Deprecated("insight::find_response()")
  insight::find_response(x, ...)
}



#' @importFrom insight find_random
#' @rdname overdisp
#' @export
grp_var <- function(x) {
  .Deprecated("insight::find_random()")
  insight::find_random(x)
}



#' @importFrom insight find_random
#' @rdname overdisp
#' @export
re_grp_var <- function(x) {
  .Deprecated("insight::find_random()")
  insight::find_random(x)
}



#' @importFrom insight clean_names
#' @rdname overdisp
#' @export
var_names <- function(x) {
  .Deprecated("insight::clean_names()")
  insight::clean_names(x)
}


#' @importFrom performance performance_accuracy
#' @rdname overdisp
#' @export
pred_accuracy <- function(x, ...) {
  .Deprecated("performance::performance_accuracy()")
  performance::performance_accuracy(x, ...)
}


#' @importFrom performance check_outliers
#' @rdname overdisp
#' @export
outliers <- function(x) {
  .Deprecated("performance::check_outliers()")
  performance::check_outliers(x)
}


#' @importFrom performance check_heteroscedasticity
#' @rdname overdisp
#' @export
heteroskedastic <- function(x) {
  .Deprecated("performance::check_heteroscedasticity()")
  performance::check_heteroscedasticity(x)
}


#' @importFrom performance check_autocorrelation
#' @rdname overdisp
#' @export
autocorrelation <- function(x) {
  .Deprecated("performance::check_autocorrelation()")
  performance::check_autocorrelation(x)
}


#' @importFrom performance check_normality
#' @rdname overdisp
#' @export
normality <- function(x) {
  .Deprecated("performance::check_normality()")
  performance::check_normality(x)
}


#' @importFrom performance check_collinearity
#' @rdname overdisp
#' @export
multicollin <- function(x) {
  .Deprecated("performance::check_collinearity()")
  performance::check_collinearity(x)
}


#' @importFrom performance check_model
#' @rdname overdisp
#' @export
check_assumptions <- function(x) {
  .Deprecated("performance::check_model()")
  performance::check_model(x)
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


#' @importFrom bayestestR hdi
#' @rdname overdisp
#' @export
hdi <- function(x) {
  .Deprecated("bayestestR::hdi()")
  bayestestR::hdi(x)
}


#' @importFrom bayestestR rope
#' @rdname overdisp
#' @export
rope <- function(x) {
  .Deprecated("bayestestR::rope()")
  bayestestR::rope(x)
}
