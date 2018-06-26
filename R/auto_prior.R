#' @importFrom stats sd
auto_prior <- function(formula, data, gaussian, locations = NULL) {

  if (!requireNamespace("brms", quietly = TRUE))
    stop("Package `brms` required.", call. = FALSE)

  scale.b <- 2.5
  scale.y <- 10

  pred <- pred_vars(formula)
  y.name <- resp_var(formula)
  y <- data[[y.name]]

  if (gaussian) {
    scale.factor <- stats::sd(y, na.rm = TRUE)
    scale.b <- scale.b * scale.factor
    scale.y <- scale.y * scale.factor
  }

  if (!is.null(locations))
    location.y <- locations[1]
  else
    location.y <- 0

  priors <- brms::set_prior(
    sprintf("normal(%f, %f)", round(location.y, 2), round(scale.y, 2)),
    class = "Intercept"
  )

  is.fac <- NULL
  term.names <- NULL
  scale.pred <- NULL

  # we need to check which predictors are categorical and then "mimic"
  # their coefficient name as it is represented in the model (i.e. variable
  # name + category name)

  for (i in pred) {
    f <- data[[i]]
    if (is.factor(f)) {
      i <- sprintf("%s%s", i, levels(f)[2:nlevels(f)])
      is.fac <- c(is.fac, rep(TRUE, nlevels(f) - 1))
      scale.pred <- c(scale.pred, rep(scale.b, nlevels(f) - 1))
    } else {
      is.fac <- c(is.fac, FALSE)
      scale.pred <- c(scale.pred, scale.b / stats::sd(f, na.rm = TRUE))
    }
    term.names <- c(term.names, i)
  }

  for (i in 1:length(term.names)) {

    if (!is.null(locations) && length(locations) >= (i + 1))
      location.b <- locations[i + 1]
    else
      location.b <- 0

    priors <- priors + brms::set_prior(
      sprintf("normal(%f, %f)", round(location.b, 2), round(scale.pred[i]), 2),
      class = "b",
      coef = term.names[i]
    )
  }

  priors
}
