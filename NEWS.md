# sjstats 0.8.0

## New functions

* `svy()` to compute robust standard errors for weighted models, adjusting the residual degrees of freedom to simulate sampling weights.
* `zero_count()` to check whether a poisson-model is over- or underfitting zero-counts in the outcome.
* `pred_accuracy()` to calculate accuracy of predictions from model fit.
* `outliers()` to detect outliers in (generalized) linear models.
* `heteroskedastic()` to check linear models for (non-)constant error variance.
* `autocorrelation()` to check linear models for auto-correlated residuals.
* `normality()` to check whether residuals in linear models are normally distributed or not.
* `multicollin()` to check predictors in a model for multicollinearity.
* `check_assumptions()` to run a set of model assumption checks.

## Changes to functions

* `prop()` no longer works within dplyr's `summarise()` function. Instead, when now used with grouped data frames, a summary of proportions is directly returned as tibble.
* `se()` now computes adjusted standard errors for generalized linear (mixed) models, using the Taylor series-based delta method. 


# sjstats 0.7.1

## General

* Package depends on R-version >= 3.3.

## Changes to functions

* `prop()` gets a `digits`-argument to round the return value to a specific number of decimal places.

# sjstats 0.7.0

## General

* Largely revised the documentation.

## New functions

* `prop()` to calculate proportion of values in a vector.
* `mse()` to calculate the mean square error for models.
* `robust()` to calculate robust standard errors and confidence intervals for regression models, returned as tidy data frame.

# sjstats 0.6.0

## New functions

* `split_half()` to compute the split-half-reliability of tests or questionnaires.
* `sd_pop()` and `var_pop()` to compute population variance and population standard deviation.

## Changes to functions

* `se()` now also computes the standard error from estimates (regression coefficients) and p-values.

# sjstats 0.5.0

## New functions

* Added S3-`print`-method for `mwu()`-function.
* `get_model_pval()` to return a tidy data frame (tibble) of model term names, p-values and standard errors from various regression model types.
* `se_ybar()` to compute standard error of sample mean for mixed models, considering the effect of clustering on the standard error.
* `std()` and `center()` to standardize and center variables, supporting the pipe-operator.

## Changes to functions

* `se()` now also computes the standard error for intraclass correlation coefficients, as returned by the `icc()`-function.
* `std_beta()` now always returns a tidy data frame (tibble) with model term names, standardized estimate, standard error and confidence intervals.
* `r2()` now also computes alternative omega-squared-statistics, if null model is given.
