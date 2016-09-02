# sjstats 0.5.0

## New functions

* Added S3-`print`-method for `mwu()`-function.
* `get_model_pval()` to return a tidy data frame (tibble) of model term names, p-values and standard errors from various regression model types.

## Changes to functions

* `se` now also computes the standard error for intraclass correlation coefficients, as returned by the `icc()`-function.

# sjstats 0.4.0

## New functions
* `inequ_trend` to calculate proportional change of absolute and relative inequalities between two status groups for a vector of given prevalence rates.


## Changes to functions

* `bootstrap` is now much more memory efficient due to use of pointers (thanks to [Hadley Wickham](https://twitter.com/hadleywickham) for the hint).
* `boot_ci`, `boot_se` and `boot_p` now accept multiple variables as input.
* `resp_val` now also applies to models fitted with `nlme::lme`.
