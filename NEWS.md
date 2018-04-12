# sjstats 0.14.3

## General

* Remove _tidyverse_ from suggested packages, as requested by maintainers.

## New functions

* `model_family()` to get model-information about family and link-functions. This function is intended to be "generic" and work with many different model objects, because not all packages provide a `family()` function.

## Bug fixes

* Fix issue with `p_value()` for unconditional mixed models.
* Fix typo in `xtab_statistics()`.
* Fix issue with wrong calculation of Nagelkerke's r-squared value in `r2()`.
* Fix issue for factors with character leves in `typical_value()`, when argument `fun` for factors was set to `mode`.
* Don't show prior-samples in `hdi()`, `tidy_stan()` etc. for _brmsfit_-objects.

# sjstats 0.14.2

## General

* Revise examples, vignettes and package description to make sure all used packages are available for CRAN checks on operating systems.

# sjstats 0.14.2

## New functions

* `residuals.svyglm.nb()` as S3-generic `residuals()` method for objects fitted with `svyglm.nb()`.

## Changes to functions

* `icc()` gets a `posterior`-argument, to compute ICC-values from `brmsfit`-objects, for the whole posterior distribution.
* `icc()` now gives a warning when computed for random-slope-intercept models, to warn user about probably inappropriate inference.
* `r2()` now computes Bayesian version of R-squared for `stanreg` and `brmsfit` objects.
* Argument `prob` in `hdi()` now accepts a vector of scalars to compute HDIs for multiple probability tresholds at once.
* Argument `probs` in `tidy_stan()` was renamed into `prob`, to be consistent with `hdi()`.
* `mwu()` gets an `out`-argument, to print output to console, or as HTML table in the viewer or web browser.
* `scale_weights()` now also works if weights have missing values.
* `hdi()` and `rope()` get `data.frame`-methods.
* `omega_sq()` and `eta_sq()` get a `ci.lvl`-argument to compute confidence intervals for the effect size statistics.
* `omega_sq()`, `eta_sq()` and `cohens_f()` now always return a data frame with at least two columns: term name and effect size. Confidence intervals are added as additional columns, if the `ci.lvl`-argument is `TRUE`.
* `omega_sq()` gets a `partial`-argument to compute partial omega-squared.
* `omega_sq()`, `eta_sq()`, `cohens_f()` and `anova_stats()` now support `anova.rms`-objects from the *rms*-package.

## Bug fixes

* Fix unnecessary warning for tibbles in `mic()`.
* Make sure that `model_frame()` does not return duplicated column names.
* Fix issue in `tidy_stan()` with incorrect *n_eff* statistics for _sigma_ parameter in mixed models.
* Fix issue in `tidy_stan()`, which did not work when `probs` was of length greater than 2.
* Fix issue in `icc()` with _brmsfit_-models, which was broken probably due to internal changes in _brms_.

# sjstats 0.14.1

## General

* Remove unused imports.
* Cross refences from `dplyr::select_helpers` were updated to `tidyselect::select_helpers`.

## Changes to functions

* `var_names()` now also cleans variable names from variables modeled with the `mi()` function (multiple imputation on the fly in *brms*).
* `reliab_test()` gets an `out`-argument, to print output to console, or as HTML table in the viewer or web browser.

## Bug fixes

* Fix issues with `mcse()`, `n_eff()` and `tidy_stan()` with more complex _brmsfit_-models.
* Fix issue in `typical_value()` to prevent error for R-oldrel-Windows.
* `model_frame()` now returns response values from models, which are in matrix form (bound with `cbind()`), as is.
* Fixed issues in `grpmean()`, where values instead of value labels were printed if some categories were not present in the data.

# sjstats 0.14.0

## General

* Beautiful colored output for `grpmean()` and `mwu()`.

## New functions

* `mcse()` to compute the Monte Carlo standard error for `stanreg`- and `brmsfit`-models.
* `n_eff()` to compute the effective sample size for `stanreg`- and `brmsfit`-models.

## Changes to functions

* `grpmean()` now uses `contrasts()` from package *emmeans* to compute p-values, which correclty indicate whether the sub-group mean is significantly different from the total mean.
* `grpmean()` gets an `out`-argument, to print output to console, or as HTML table in the viewer or web browser.
* `tidy_stan()` now includes information on the Monte Carlo standard error.
* `model_frame()`, `p_value()` and `link_inverse()` now support Zelig-relogit-models.
* `typical_value()` gets an explicit `weight.by`-argument.

## Bug fixes

* `model_frame()` did not work properly for variables that were standardized with `scale()`.
* In certain cases, `weight.by`-argument did not work in `grpmean()`.

# sjstats 0.13.0

## General

* Remove deprecated `get_model_pval()`.
* Revised documentation for `overdisp()`.

## New functions

* `scale_weights()` to rescale design weights for multilevel models.
* `pca()` and `pca_rotate()` to create tidy summaries of principal component analyses or rotated loadings matrices from PCA.
* `gmd()` to compute Gini's mean difference.
* `is_prime()` to check whether a number is a prime number or not.

## Changes to functions

* `link_inverse()` now supports `brmsfit`, `multinom` and `clm`-models.
* `p_value()` now supports `polr` and `multinom`-models.
* `zero_count()` gets a `tolerance`-argument, to accept models with a ratio within a certain range of 1.
* `var_names()` now also cleans variable names from variables modelled with the `offset()`, `lag()` or `diff()` function.
* `icc()`, `re_var()` and `get_re_var()` now support `brmsfit`-objects (models fitted with the *brms*-package).
* For `fun = "weighted.mean"`, `typical_value()` now checks if vector of weights is of same length as `x`.
* The print-method for `grpmean()` now also prints the overall p-value from the model.

## Bug fixes

* `resp_val()`, `cv_error()` and `pred_accuracy()` did not work for formulas with transforming function for response terms, e.g. `log(response)`.


# sjstats 0.12.0

## General

* Fixed examples, to resolve issues with CRAN package checks.
* More model objects supported in `p_value()`.

## New functions

* `model_frame()` to get the model frame from model objects, also of those models that don't have a S3-generic model.frame-function.
* `var_names()` to get cleaned variable names from model objects.
* `link_inverse()` to get the inverse link function from model objects.

## Changes to functions

* The `fun`-argument in `typical_value()` can now also be a named vector, to apply different functions for numeric and categorical variables.

## Bug fixes

* Fixed issue with specific model formulas in `pred_vars()`.
* Fixed issue with specific model objects in `resp_val()`.
* Fixed issue with nested models in `re_var()`.

# sjstats 0.11.2

## New functions

* `tidy_stan()` to return a tidy summary of Stan-models.

## Changes to functions

* `hdi()` and `rope()` now also work for `brmsfit`-models, from package *brms*.
* `hdi()` and `rope()` now have a `type`-argument, to return fixed, random or all effects for mixed effects models.

# sjstats 0.11.1

## Changes to functions

* `typical_value()` gets a "zero"-option for the `fun`-argument.
* Changes to `icc()`, which used `stats::sigma()` and thus required R-version 3.3 or higher. Now should depend on R 3.2 again.
* `se()` now also supports `stanreg` and `stanfit` objects.
* `hdi()` now also supports `stanfit`-objects.
* `std_beta()` gets a `ci.lvl`-argument, to specify the level of the calculated confidence interval for standardized coefficients.
* `get_model_pval()` is now deprecated. Please use `p_value()` instead.

## New functions

* `rope()` to calculate the region of practical equivalence for MCMC samples.

# sjstats 0.11.0

## General

* Added vignettes for various functions.
* Fixed issue with latest tidyr-update on CRAN.

## New functions

* `grpmean()` to compute mean values by groups (One-way Anova).
* `hdi()` to compute high density intervals (HDI) for MCMC samples.
* `find_beta()` and `find_beta2()` to find the shape parameters of a Beta distribution.
* `find_normal()` and `find_cauchy()` to find the parameters of a normal or cauchy distribution.

# sjstats 0.10.3

## New functions

* `typical_value()`, to return the typical value of a variable.
* `eta_sq()`, `cohens_f()` and `omega_sq()` to compute (partial) eta-squared or omega-squared statistics, or Cohen's F for anova tables.
* `anova_stats()` to compute a complete model summary, including (partial) eta-squared, omega-squared and Cohen's F statistics for anova tables, returned as tidy data frame.
* `svy_md()` as convenient shortcut to compute the median for variables from survey designs.
* `is_singular()` to check a model fit for singularity in case of post-fitting convergence warnings.

## Changes to functions

* Computation of `r2()` for `glm`-objects is now based on log-Likelihood methods and also accounts for count models.
* Better `print()`-method for `overdisp()`.
* `print()`-method for `svyglm.nb()` now also prints the dispersion parameter Theta.
* `overdisp()` now supports `glmmTMB`-objects.
* `boot_ci()` also displays CI based on sample quantiles.

## Bug fixes

* `std_beta()` did not work for models with only one predictor.

# sjstats 0.10.2

## Changes to functions

* `icc()`, `re_var()` and `get_re_var()` now support `glmmTMB`-objects.
* `pred_accuracy()` now also reports the standard error of accuracy, and gets a print-method.

## Bug fixes

* `pred_accuracy()` with cross-validation-method did not correctly account for the generated test data.
* Fixed issue with calculation in `smpsize_lmm()` and `se_ybar()`.

# sjstats 0.10.1

## General

* Revised imports: Labelled data functions from package *sjmisc* have been moved to package *sjlabelled*.

## New functions

* `boot_est()` to return the estimate from bootstrap replicates.

## Changes to functions

* The `print()`-method for `svyglm.nb()`-objects now also prints confidence intervals.

## Bug fixes

* `se()` did not work for `icc()`-objects, when the mixed model had more than one random effect term.

# sjstats 0.10.0

## New functions

* `cv_error()` and `cv_compare()` to compute the root mean squared error for test and training data from cross-validation.
* `props()` to calculate proportions in a vector, supporting multiple logical statements.
* `or_to_rr()` to convert odds ratio estimates into risk ratio estimates.
* `mn()`, `md()` and `sm()` to calculate mean, median or sum of a vector, but using `na.rm = TRUE` as default.
* S3-generics for `svyglm.nb`-models: `family()`, `print()`, `formula()`, `model.frame()` and `predict()`.

## Bug fixes

* Fixed error in computation of `mse()`.

# sjstats 0.9.0

## General

* Functions `std()` and `center()` were removed and are now in the [sjmisc-package](https://cran.r-project.org/package=sjmisc).

## New functions

* `svyglm.nb()` to compute survey-weighted negative binomial regressions. 
* `xtab_statistics()` to compute various measures of assiciation for contingency tables.
* Added S3-`model.frame()`-function for `gee`-models.

## Changes to functions

* `se()` gets a `type`-argument, which applies to generalized linear mixed models. You can now choose to compute either standard errors with delta-method approximation for fixed effects only, or standard errors for joint random and fixed effects.


## Bug fixes

* `prop()` did not work for non-labelled data frames when used with grouped data frames.

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
