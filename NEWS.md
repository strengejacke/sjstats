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
