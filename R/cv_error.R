utils::globalVariables(c("rmse.train", "predicted", "residuals"))

#' @title Test and training error from model cross-validation
#' @name cv_error
#'
#' @description \code{cv_error} computes the root mean squared error from a model fitted
#'          to kfold cross-validated test-training-data. \code{cv_compare}
#'          does the same, for multiple formulas at once (by calling \code{cv_error}
#'          for each formula).
#'
#' @param data A data frame.
#' @param formula The formula to fit the linear model for the test and training data.
#' @param formulas A list of formulas, to fit linear models for the test and training data.
#'
#' @inheritParams pred_accuracy
#'
#' @return A tibble with the root mean squared errors for the training and test data.
#'
#' @details \code{cv_error} first generates cross-validated test-training pairs, using
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
#' @importFrom modelr crossv_kfold
#' @importFrom dplyr mutate summarise
#' @importFrom purrr map map2 map_dbl map_df
#' @importFrom broom augment
#' @importFrom tidyr unnest
#' @importFrom tibble tibble
#' @export
cv_error <- function(data, formula, k = 5) {

  # compute cross validation data
  cv_data <- data %>%
    modelr::crossv_kfold(k = k) %>%
    dplyr::mutate(
      trained.models = purrr::map(.data$train, ~ stats::lm(formula, data = .x)),
      predicted = purrr::map2(.data$trained.models, .data$test, ~ broom::augment(.x, newdata = .y)),
      residuals = purrr::map(.data$predicted, ~.x[[resp_var(formula)]] - .x$.fitted),
      rmse.train = purrr::map_dbl(.data$trained.models, ~ sjstats::rmse(.x))
    )


  # Training error
  train.error <- dplyr::summarise(cv_data, train.error = mean(rmse.train, na.rm = TRUE))


  # Test error
  test.error <- cv_data %>%
    tidyr::unnest(predicted, residuals) %>%
    dplyr::summarise(test.error = sqrt(mean(residuals ^ 2, na.rm = TRUE)))

  tibble::tibble(
    model = deparse(formula),
    train.error = round(train.error[[1]], 4),
    test.error = round(test.error[[1]], 4)
  )
}



#' @rdname cv_error
#' @export
cv_compare <- function(data, formulas, k = 5) {
  purrr::map_df(formulas, ~ cv_error(data, formula = .x, k = k))
}
