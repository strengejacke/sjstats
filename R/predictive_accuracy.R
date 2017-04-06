utils::globalVariables(c("predictions", "response", "train", "test"))

#' @title Accuracy of predictions from model fit
#' @name pred_accuracy
#'
#' @description This function calculates the predictive accuracy of (generalized)
#'              linear models.
#'
#' @param fit Fitted model object of class \code{lm} or \code{glm}.
#' @param k The number of folds for the kfold-crossvalidation.
#' @param method Character string, indicating whether crossvalidation
#'          (\code{method = "cv"}) or bootstrapping (\code{method = "boot"})
#'          is used to compute the accuracy values.
#'
#' @inheritParams bootstrap
#'
#' @return The accuracy of the model predictions, i.e. the proportion of
#'         accurately predicted values from the model.
#'
#' @details For linar models, the accuracy is the correlation coefficient
#'          between the actual and the predicted value of the outcome. For
#'          logistic regression models, the accuracy corresponds to the
#'          AUC-value, calculated with the \code{\link[pROC]{auc}}-function.
#'          \cr \cr
#'          The accuracy is the mean value of multiple correlation resp.
#'          AUC-values, which are either computed with crossvalidation
#'          or nonparametric bootstrapping (see argument \code{method}).
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

    # check if bootstrapping or cross validation is requested
    if (method == "boot") {

      # accuracy linear models with bootstrapping
      cv <- data %>%
        bootstrap(n) %>%
        dplyr::mutate(
          models = purrr::map(strap, ~ stats::lm(formula, data = .x)),
          predictions = purrr::map(models, ~ stats::predict(.x, type = "response")),
          response = purrr::map(models, ~ resp_val(.x)),
          accuracy = purrr::map2_dbl(predictions, response, ~ stats::cor(.x, .y, use = "pairwise.complete.obs"))
        )
    } else {

      # accuracy linear models with cross validation
      cv <- modelr::crossv_kfold(data, k = k) %>%
        dplyr::mutate(
          models = purrr::map(train, ~ stats::lm(formula, data = .x)),
          predictions = purrr::map2(test, models, ~ modelr::add_predictions(.x, .y)[["pred"]]),
          response = purrr::map(test, ~ as.data.frame(.x)[[resp.name]]),
          accuracy = purrr::map2_dbl(predictions, response, ~ stats::cor(.x, .y, use = "pairwise.complete.obs"))
        )
    }

  } else if (inherits(fit, "glm") && stats::family(fit) == "binomial") {

    # check if bootstrapping or cross validation is requested
    if (method == "boot") {

      # accuracy linear models with bootstrapping
      cv <- data %>%
        bootstrap(n) %>%
        dplyr::mutate(
          models = purrr::map(strap, ~ stats::lm(formula, data = .x)),
          accuracy = purrr::map_dbl(models, ~ pROC::auc(pROC::roc(response = resp_val(.x), predictor = stats::predict.glm(.x, stats::model.frame(.x)))))
        )

    } else {

      # accuracy logistic regression models with cross validation
      cv <- modelr::crossv_kfold(data, k = k) %>%
        dplyr::mutate(
          models = purrr::map(train, ~ stats::glm(formula, data = ., family = stats::binomial(link = "logit"))),
          accuracy = purrr::map_dbl(models, ~ pROC::auc(pROC::roc(response = resp_val(.x), predictor = stats::predict.glm(.x, stats::model.frame(.x)))))
        )
    }
  }

  # return mean value of accuracy
  mean(cv$accuracy)
}
