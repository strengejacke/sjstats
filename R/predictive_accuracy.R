#' @title Accuracy of predictions from model fit
#' @name pred_accuracy
#'
#' @seealso \code{\link{cv_error}}
#'
#' @description This function calculates the predictive accuracy of linear
#'              or logistic regression models.
#'
#' @param fit Fitted model object of class \code{lm} or \code{glm}, the latter
#'          being a logistic regression model (binary response).
#' @param k The number of folds for the kfold-crossvalidation.
#' @param method Character string, indicating whether crossvalidation
#'          (\code{method = "cv"}) or bootstrapping (\code{method = "boot"})
#'          is used to compute the accuracy values.
#'
#' @inheritParams bootstrap
#'
#' @return A list with two values: The \code{accuracy} of the model predictions, i.e.
#'         the proportion of accurately predicted values from the model and
#'         its standard error, \code{std.error}.
#'
#' @details For linar models, the accuracy is the correlation coefficient
#'          between the actual and the predicted value of the outcome. For
#'          logistic regression models, the accuracy corresponds to the
#'          AUC-value, calculated with the \code{\link[pROC]{auc}}-function.
#'          \cr \cr
#'          The accuracy is the mean value of multiple correlation resp.
#'          AUC-values, which are either computed with crossvalidation
#'          or nonparametric bootstrapping (see argument \code{method}).
#'          The standard error is the standard deviation of the computed
#'          correlation resp. AUC-values.
#'
#' @examples
#' data(efc)
#' fit <- lm(neg_c_7 ~ barthtot + c161sex, data = efc)
#'
#' # accuracy for linear model, with crossvalidation
#' pred_accuracy(efc, fit)
#'
#' # accuracy for linear model, with bootstrapping
#' pred_accuracy(efc, fit, method = "boot", n = 100)
#'
#' # accuracy for logistic regression, with crossvalidation
#' efc$services <- sjmisc::dicho(efc$tot_sc_e, dich.by = 0, as.num = TRUE)
#' fit <- glm(services ~ neg_c_7 + c161sex + e42dep,
#'            data = efc, family = binomial(link = "logit"))
#' pred_accuracy(efc, fit)
#'
#' @importFrom modelr crossv_kfold add_predictions
#' @importFrom purrr map map2 map2_dbl
#' @importFrom stats lm cor glm predict predict.glm model.frame formula binomial
#' @importFrom dplyr mutate
#' @export
pred_accuracy <- function(data, fit, method = c("cv", "boot"), k = 5, n = 1000) {

  method <- match.arg(method)

  # check package availability
  if (!requireNamespace("pROC", quietly = TRUE) && inherits(fit, "glm")) {
    stop("Package `pROC` needed for this function to work. Please install it.", call. = F)
  }

  # get formula from model fit
  formula <- stats::formula(fit)

  # get name of response
  resp.name <- resp_var(fit)

  # accuracy for linear models
  if (inherits(fit, "lm") && !inherits(fit, "glm")) {

    measure <- "Correlation between observed and predicted"

    # check if bootstrapping or cross validation is requested
    if (method == "boot") {

      # accuracy linear models with bootstrapping
      cv <- data %>%
        bootstrap(n) %>%
        dplyr::mutate(
          models = purrr::map(.data$strap, ~ stats::lm(formula, data = .x)),
          predictions = purrr::map(.data$models, ~ stats::predict(.x, type = "response")),
          response = purrr::map(.data$models, ~ resp_val(.x)),
          accuracy = purrr::map2_dbl(.data$predictions, .data$response, ~ stats::cor(.x, .y, use = "pairwise.complete.obs"))
        )
    } else {

      # accuracy linear models with cross validation
      cv <- modelr::crossv_kfold(data, k = k) %>%
        dplyr::mutate(
          models = purrr::map(.data$train, ~ stats::lm(formula, data = .x)),
          predictions = purrr::map2(.data$test, .data$models, ~ modelr::add_predictions(.x, .y)[["pred"]]),
          response = purrr::map(.data$test, ~ as.data.frame(.x)[[resp.name]]),
          accuracy = purrr::map2_dbl(.data$predictions, .data$response, ~ stats::cor(.x, .y, use = "pairwise.complete.obs"))
        )
    }

  } else if (inherits(fit, "glm") && stats::family(fit) == "binomial") {

    measure <- "Area under Curve"

    # check if bootstrapping or cross validation is requested
    if (method == "boot") {

      # accuracy linear models with bootstrapping
      cv <- data %>%
        bootstrap(n) %>%
        dplyr::mutate(
          models = purrr::map(.data$strap, ~ stats::glm(formula, data = .x, family = stats::binomial(link = "logit"))),
          accuracy = purrr::map_dbl(.data$models, ~ pROC::auc(pROC::roc(response = resp_val(.x), predictor = stats::predict.glm(.x, stats::model.frame(.x)))))
        )

    } else {

      # accuracy logistic regression models with cross validation
      cv <- modelr::crossv_kfold(data, k = k) %>%
        dplyr::mutate(
          models = purrr::map(.data$train, ~ stats::glm(formula, data = .x, family = stats::binomial(link = "logit"))),
          predictions = purrr::map2(.data$models, .data$test, ~ stats::predict.glm(.x, newdata = .y)),
          response = purrr::map(.data$test, ~ as.data.frame(.x)[[resp.name]]),
          accuracy = purrr::map2_dbl(.data$response, .data$predictions, ~ pROC::auc(pROC::roc(response = .x, predictor = .y)))
        )
    }
  }

  # return mean value of accuracy
  structure(
    class = c("sjstats_pred_accuracy", "list"),
    list(
      accuracy = mean(cv$accuracy),
      std.error = sd(cv$accuracy),
      stat = measure
    )
  )
}
