#' @title Create default priors for brms-models
#' @name auto_prior
#'
#' @description This function creates default priors for brms-regression
#'   models, based on the same automatic prior-scale adjustment as in
#'   \pkg{rstanarm}.
#'
#' @param formula A formula describing the model, which just needs to contain
#'   the model terms, but no notation of interaction, splines etc. Usually,
#'   you want only those predictors in the formula, for which automatic
#'   priors should be generated. Add informative priors afterwards to the
#'   returned \code{brmsprior}-object.
#' @param data The data that will be used to fit the model.
#' @param gaussian Logical, if the outcome is gaussian or not.
#' @param locations A numeric vector with location values for the priors. If
#'   \code{locations = NULL}, \code{0} is used as location parameter.
#'
#' @return A \code{brmsprior}-object.
#'
#' @details \code{auto_prior()} is a small, convenient function to create
#'   some default priors for brms-models with automatically adjusted prior
#'   scales, in a similar way like \pkg{rstanarm} does. The default scale for
#'   the intercept is 10, for coefficients 2.5. If the outcome is gaussian,
#'   both scales are multiplied with \code{sd(y)}. Then, for categorical
#'   variables, nothing more is changed. For numeric variables, the scales
#'   are divided by the standard deviation of the related variable.
#'   \cr \cr
#'   All prior distributions are \emph{normal} distributions. \code{auto_prior()}
#'   is intended to quickly create default priors with feasible scales. If
#'   more precise definitions of priors is necessary, this needs to be done
#'   directly with brms-functions like \code{set_prior()}.
#'
#' @note As \code{auto_prior()} also sets priors on the intercept, the model
#'   formula used in \code{brms::brm()} must be rewritten to something like
#'   \code{y ~ 0 + intercept ...}, see \code{\link[brms]{set_prior}}.
#'
#' @examples
#' library(sjmisc)
#' data(efc)
#' efc$c172code <- as.factor(efc$c172code)
#' efc$c161sex <- to_label(efc$c161sex)
#'
#' mf <- formula(neg_c_7 ~ c161sex + c160age + c172code)
#'
#' if (requireNamespace("brms", quietly = TRUE))
#'   auto_prior(mf, efc, TRUE)
#'
#' ## compare to
#' # library(rstanarm)
#' # m <- stan_glm(mf, data = efc, chains = 2, iter = 200)
#' # ps <- prior_summary(m)
#' # ps$prior_intercept$adjusted_scale
#' # ps$prior$adjusted_scale
#'
#' ## usage
#' # ap <- auto_prior(mf, efc, TRUE)
#' # brm(mf, data = efc, priors = ap)
#'
#' # add informative priors
#' mf <- formula(neg_c_7 ~ c161sex + c172code)
#'
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   auto_prior(mf, efc, TRUE) +
#'     brms::prior(normal(.1554, 40), class = "b", coef = "c160age")
#' }
#'
#' # example with binary response
#' efc$neg_c_7d <- ifelse(efc$neg_c_7 < median(efc$neg_c_7, na.rm = TRUE), 0, 1)
#' mf <- formula(neg_c_7d ~ c161sex + c160age + c172code + e17age)
#'
#' if (requireNamespace("brms", quietly = TRUE))
#'   auto_prior(mf, efc, FALSE)
#'
#' @importFrom stats sd na.omit
#' @importFrom dplyr select n_distinct
#' @importFrom insight find_predictors find_response
#' @export
auto_prior <- function(formula, data, gaussian, locations = NULL) {

  if (!requireNamespace("brms", quietly = TRUE))
    stop("Package `brms` required.", call. = FALSE)

  scale.b <- 2.5
  scale.y <- 10

  pred <- insight::find_predictors(formula, effects = "all", flatten = TRUE)
  y.name <- insight::find_response(formula, combine = TRUE)

  cols <- c(y.name, pred)

  data <- data %>%
    dplyr::select(!! cols) %>%
    stats::na.omit() %>%
    as.data.frame()

  y <- data[[y.name]]

  # check if response is binary
  if (missing(gaussian) && dplyr::n_distinct(y, na.rm = TRUE) == 2) gaussian <- FALSE

  if (isTRUE(gaussian) && dplyr::n_distinct(y, na.rm = TRUE) == 2)
    warning("Priors were calculated based on assumption that the response is Gaussian, however it seems to be binary.", call. = F)


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
    sprintf("normal(%s, %s)", round(location.y, 2), round(scale.y, 2)),
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
      sprintf("normal(%s, %s)", round(location.b, 2), round(scale.pred[i], 2)),
      class = "b",
      coef = term.names[i]
    )
  }

  priors
}
