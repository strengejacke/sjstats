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

# sjstats 0.4.0

## New functions
* `inequ_trend()` to calculate proportional change of absolute and relative inequalities between two status groups for a vector of given prevalence rates.


## Changes to functions

* `bootstrap()` is now much more memory efficient due to use of pointers (thanks to [Hadley Wickham](https://twitter.com/hadleywickham) for the hint).
* `boot_ci()`, `boot_se()` and `boot_p()` now accept multiple variables as input.
* `resp_val()` now also applies to models fitted with `nlme::lme()`.
