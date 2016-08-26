# sjstats 0.4.0

## New functions
* `inequ_trend` to calculate proportional change of absolute and relative inequalities between two status groups for a vector of given prevalence rates.


## Changes to functions

* `bootstrap` is now much more memory efficient due to use of pointers (thanks to [Hadley Wickham](https://twitter.com/hadleywickham) for the hint).
* `boot_ci`, `boot_se` and `boot_p` now accept multiple variables as input.
* `resp_val` now also applies to models fitted with `nlme::lme`.


# sjstats 0.3.0

## General

* Removed nom-necessary checks for package-availability.

## New functions

* `bootstrap` to generate bootstrap replicates of data frames.
* `boot_ci` to compute confidence intervals from bootstrapped values.
* `pred_vars` to get the names of predictor variables from fitted models.
* `resp_var` to get the name of the response variable from fitted models.
* `resp_val` to get the values of the response vector from fitted models.
