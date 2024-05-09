#' @title Test and training error from model cross-validation
#' @name cv_error
#'
#' @description \code{cv_error()} computes the root mean squared error from a model fitted
#'          to kfold cross-validated test-training-data. \code{cv_compare()}
#'          does the same, for multiple formulas at once (by calling \code{cv_error()}
#'          for each formula).
#'
#' @param data A data frame.
#' @param formula The formula to fit the linear model for the test and training data.
#' @param formulas A list of formulas, to fit linear models for the test and training data.
#' @param k The number of folds for the kfold-crossvalidation.
#'
#' @return A data frame with the root mean squared errors for the training and test data.
#'
#' @details \code{cv_error()} first generates cross-validated test-training pairs, using
#'          \code{\link[modelr]{crossv_kfold}} and then fits a linear model, which
#'          is described in \code{formula}, to the training data. Then, predictions
#'          for the test data are computed, based on the trained models.
#'          The \emph{training error} is the mean value of the \code{\link{rmse}} for
#'          all \emph{trained} models; the \emph{test error} is the rmse based on all
#'          residuals from the test data.
#'
#' @examples
#' data(efc)
#' cv_error(efc, neg_c_7 ~ barthtot + c161sex)
#'
#' cv_compare(efc, formulas = list(
#'   neg_c_7 ~ barthtot + c161sex,
#'   neg_c_7 ~ barthtot + c161sex + e42dep,
#'   neg_c_7 ~ barthtot + c12hour
#' ))
#'
#' @export
cv_error <- function(data, formula, k = 5) {
  insight::check_if_installed("datawizard")

  # response
  resp <- insight::find_response(formula)

  # compute cross validation data
  cv_data <- lapply(k, function(i) {
    datawizard::data_partition(data, proportion = 0.8)
  })
  # get train and test datasets
  train_data <- lapply(cv_data, function(cvdat) cvdat[[1]])
  test_data <- lapply(cv_data, function(cvdat) cvdat[[2]])

  # fit models on datasets
  trained_models <- lapply(train_data, function(x) stats::lm(formula, data = x))
  test_models <- lapply(test_data, function(x) stats::lm(formula, data = x))

  # RMSE
  train_error <- mean(vapply(trained_models, performance::rmse, numeric(1)), na.rm = TRUE)
  test_error <- mean(vapply(test_models, performance::rmse, numeric(1)), na.rm = TRUE)

  data_frame(
    model = deparse(formula),
    train.error = round(train_error, 4),
    test.error = round(test_error, 4)
  )
}



#' @rdname cv_error
#' @export
cv_compare <- function(data, formulas, k = 5) {
  do.call(rbind, lapply(formulas, function(f) cv_error(data, formula = f, k = k)))
}
