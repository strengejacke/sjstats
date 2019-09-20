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



#' @rdname samplesize_mixed
#' @export
smpsize_lmm <- function(eff.size, df.n = NULL, power = .8, sig.level = .05, k, n, icc = 0.05) {
  .Deprecated("samplesize_mixed()")
  samplesize_mixed(eff.size, df.n, power, sig.level, k = k, n, icc)
}


#' @rdname design_effect
#' @export
deff <- function(n, icc = 0.05) {
  .Deprecated("design_effect()")
  design_effect(n, icc)
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


#' @importFrom performance cronbachs_alpha
#' @rdname overdisp
#' @export
cronb <- function(x, ...) {
  .Deprecated("performance::cronbachs_alpha()")
  performance::cronbachs_alpha(x)
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



#' @importFrom performance performance_accuracy
#' @rdname overdisp
#' @export
pred_accuracy <- function(x, ...) {
  .Deprecated("performance::performance_accuracy()")
  performance::performance_accuracy(x, ...)
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
