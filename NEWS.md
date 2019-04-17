# sjstats 0.17.5

## New functions

* `epsilon_sq()`, to compute epsilon-squared effect-size.

## Deprecated and defunct

Following functions are deprecated:

* `link_inverse()`, please use `insight::link_inverse()`
* `model_family()`, please use `insight::model_info()`
* `model_frame()`, please use `insight::get_data()`
* `pred_vars()`, please use `insight::find_predictors()`
* `re_grp_var()`, please use `insight::find_random()`
* `grp_var()`, please use `insight::find_random()`
* `resp_val()`, please use `insight::get_response()`
* `resp_var()`, please use `insight::find_response()`
* `var_names()`, please use `insight::clean_names()`
* `overdisp()`, please use `performance::check_overdispersion()`
* `zero_count()`, please use `performance::check_zeroinflation()`
* `converge_ok()`, please use `performance::check_convergence()`
* `is_singular()`, please use `performance::check_singularity()`
* `reliab_test()`, please use `performance::item_reliability()`
* `split_half()`, please use `performance::item_split_half()`
* `cronb()`, please use `performance::cronbachs_alpha()`
* `difficulty()`, please use `performance::item_difficulty()`
* `mic()`, please use `performance::item_intercor()`
* `pca()`, please use `performance::principal_components()`
* `pca_rotate()`, please use `performance::principal_components()`
* `r2()`, please use `performance::r2()`
* `icc()`, please use `performance::icc()`
* `rmse()`, please use `performance::rmse()`
* `rse()`, please use `performance::rse()`
* `mse()`, please use `performance::mse()`
* `hdi()`, please use `bayestestR::hdi()`
* `cred_int()`, please use `bayestestR::ci()`
* `rope()`, please use `bayestestR::rope()`
* `equi_test()`, please use `bayestestR::equivalence_test()`

## Changes to functions

* Anova-stats functions (like `eta_sq()`) get a `method`-argument to define the method for computing confidence intervals from bootstrapping.

## Bug fixes

* In some situations, `smpsize_lmm()` could result in negative sample-size recommendations. This was fixed, and a warning is now shown indicating that the parameters for the power-calculation should be modified.

# sjstats 0.17.4

## General

* Following models/objects are now supported by model-information functions like `model_family()`, `link_inverse()` or `model_frame()`: `MixMod` (package **GLMMadaptive**), **MCMCglmm**, `mlogit` and `gmnl`.
* Reduce package dependencies.

## New functions

* `cred_int()`, to compute uncertainty intervals of Bayesian models. Mimics the behaviour and style of `hdi()` and is thus a convenient complement to functions like `posterior_interval()`.

## Changes to functions

* `equi_test()` now finds better defaults for models with binomial outcome (like logistic regression models).
* `r2()` for mixed models now also should work properly for mixed models fitted with **rstanarm**.
* `anova_stats()` and alike (e.g. `eta_sq()`) now all preserve original term names.
* `model_family()` now returns `$is_count = TRUE`, when model is a count-model, and `$is_beta = TRUE` for models with beta-family.
* `pred_vars()` checks that return value has only unique values.
* `pred_vars()` gets a `zi`-argument to return the variables from a model's zero-inflation-formula.

## Bug fixes

* Fix minor issues in `wtd_sd()` and `wtd_mean()` when weight was `NULL` (which usually shoudln't be the case anyway).
* Fix potential issue with `deparse()`, cutting off very long formulas in various functions.
* Fix encoding issues in help-files.

# sjstats 0.17.3

## General

* Export `dplyr::n()`, to meet forthcoming changes in dplyr 0.8.0.

## Changes to functions

* `boot_ci()` gets a `ci.lvl`-argument.
* The `rotation`-argument in `pca_rotate()` now supports all rotations from `psych::principal()`.
* `pred_vars()` gets a `fe.only`-argument to return only fixed effects terms from mixed models, and a `disp`-argument to return the variables from a model's dispersion-formula.
* `icc()` for Bayesian models gets a `adjusted`-argument, to calculate adjusted and conditional ICC (however, only for Gaussian models).
* For `icc()` for non-Gaussian Bayes-models, a message is printed that recommends setting argument `ppd` to `TRUE`.
* `resp_val()` and `resp_var()` now also work for **brms**-models with additional response information (like `trial()` in formula).
* `resp_var()` gets a `combine`-argument, to return either the name of the matrix-column or the original variable names for matrix-columns.
* `model_frame()` now also returns the original variables for matrix-column-variables.
* `model_frame()` now also returns the variable from the dispersion-formula of **glmmTMB**-models.
* `model_family()` and `link_inverse()` now supports **glmmPQL**, **felm** and **lm_robust**-models.
* `anova_stats()` and alike (`omeqa_sq()` etc.) now support gam-models from package **gam**.
* `p_value()` now supports objects of class `svyolr`.

## Bug fixes

* Fix issue with `se()` and `get_re_var()` for objects returned by `icc()`.
* Fix issue with `icc()` for Stan-models.
* `var_names()` did not clear terms with log-log transformation, e.g. `log(log(y))`.
* Fix issue in `model_frame()` for models with splines with only one column.

# sjstats 0.17.2

## General

* Revised help-files for `r2()` and `icc()`, also by adding more references.

## New functions

* `re_grp_var()` to find group factors of random effects in mixed models.

## Changes to functions

* `omega_sq()` and `eta_sq()` give more informative messages when using non-supported objects.
* `r2()` and `icc()` give more informative warnings and messages.
* `tidy_stan()` supports printing simplex parameters of monotonic effects of **brms** models.
* `grpmean()` and `mwu()` get a `file` and `encoding` argument, to save the HTML output as file.

## Bug fixes

* `model_frame()` now correctly names the offset-columns for terms provided as `offset`-argument (i.e. for models where the offset was not specified inside the formula).
* Fixed issue with `weights`-argument in `grpmean()` when variable name was passed as character vector.
* Fixed issue with `r2()` for **glmmTMB** models with `ar1` random effects structure.

# sjstats 0.17.1

## New functions

* `wtd_chisqtest()` to compute a weighted Chi-squared test.
* `wtd_median()` to compute the weighted median of variables.
* `wtd_cor()` to compute weighted correlation coefficients of variables.

## Changes to functions

* `mediation()` can now cope with models from different families, e.g. if the moderator or outcome is binary, while the treatment-effect is continuous.
* `model_frame()`, `link_inverse()`, `pred_vars()`, `resp_var()`, `resp_val()`, `r2()` and `model_family()` now support `clm2`-objects from package **ordinal**.
* `anova_stats()` gives a more informative message for non-supported models or ANOVA-options.

## Bug fixes

* Fixed issue with `model_family()` and `link_inverse()` for models fitted with `pscl::hurdle()` or `pscl::zeroinfl()`.
* Fixed issue with wrong title in `grpmean()` for grouped data frames, when grouping variable was an unlabelled factor.
* Fix issue with `model_frame()` for **coxph**-models with polynomial or spline-terms.
* Fix issue with `mediation()` for logical variables.

# sjstats 0.17.0

## General

* Reduce package dependencies.

## New functions

* `wtd_ttest()` to compute a weighted t-test.
* `wtd_mwu()` to compute a weighted Mann-Whitney-U or Kruskal-Wallis test.

## Changes to functions

* `robust()` was revised, getting more arguments to specify different types of covariance-matrix estimation, and handling these more flexible.
* Improved `print()`-method for `tidy_stan()` for _brmsfit_-objects with categorical-families.
* `se()` now also computes standard errors for relative frequencies (proportions) of a vector.
* `r2()` now also computes r-squared values for _glmmTMB_-models from `genpois`-families.
* `r2()` gives more precise warnings for non-supported model-families.
* `xtab_statistics()` gets a `weights`-argument, to compute measures of association for contingency tables for weighted data.
* The `statistics`-argument in `xtab_statistics()` gets a `"fisher"`-option, to force Fisher's Exact Test to be used.
* Improved variance calculation in `icc()` for generalized linear mixed models with Poisson or negative binomial families.
* `icc()` gets an `adjusted`-argument, to calculate the adjusted and conditional ICC for mixed models.
* To get consistent argument names accross functions, argument `weight.by` is now deprecated and renamed into `weights`.

## Bug fixes

* Fix issues with effect size computation for repeated-measure Anova when using bootstrapping to compute confidence intervals.
* `grpmean()` now also adjusts the `n`-columm for weighted data.
* `icc()`, `re_var()` and  `get_re_var()` now correctly compute the random-effect-variances for models with multiple random slopes per random effect term (e.g., `(1 + rs1 + rs2 | grp)`).
* Fix issues in `tidy_stan()`, `mcse()`, `hdi()` and `n_eff()` for `stan_polr()`-models.
* Plotting `equi_test()` did not work for intercept-only models.

# sjstats 0.16.0

## General

* The S3-generics for functions like `hdi()`, `rope()`, `equi_test()` etc. are now more generic, and function usage for each supported object is now included in the documentation.
* Following functions are now S3-generic: `icc()`, `r2()`, `p_value()`, `se()`, and `std_beta()`.
* Added `print()`-methods for some more functions, for a clearer output.
* Revised `r2()` for mixed models (packages **lme4**, **glmmTMB**). The r-squared value should be much more precise now, and reports the marginal and conditional r-squared values.
* Reduced package dependencies and removed _apaTables_ and _MBESS_ from suggested packages
* `stanmvreg`-models are now supported by many functions.

## New functions

* `binned_resid()` to plot binned residuals for logistic regression models.
* `error_rate()` to compute model quality for logistic regression models.
* `auto_prior()` to quickly create automatically adjusted priors for brms-models.
* `difficulty()` to compute the item difficulty.

## Changes to functions

* `icc()` gets a `ppd`-argument for Stan-models (*brmsfit* and *stanreg*), which performs a variance decomposition based on the posterior predictive distribution. This is the recommended way for non-Gaussian models.
* For Stan-models (*brmsfit* and *stanreg*), `icc()` now also computes the HDI for the ICC and random-effect variances. Use the `prob`-argument to specify the limits of this interval.
* `link_inverse()` and `model_family()` now support _clmm_-models (package *ordinal*) and _glmRob_ and _lmRob_-models (package *robust*).
* `model_family()` gets a `multi.resp`-argument, to return a list of family-informations for multivariate-response models (of class `brmsfit` or `stanmvreg`).
* `link_inverse()` gets a `multi.resp`-argument, to return a list of link-inverse-functions for multivariate-response models (of class `brmsfit` or `stanmvreg`).
* `p_value()` now supports _rlm_-models (package *MASS*).
* `check_assumptions()` for single models with `as.logical = FALSE` now has a nice print-method.
* `eta_sq()` and `omega_sq()` now also work for repeated-measure Anovas, i.e. Anova with error term (requires broom > 0.4.5).

## Bug fixes

* `model_frame()` and `var_names()` now correctly cleans nested patterns like `offset(log(x + 10))` from column names.
* `model_frame()` now returns proper column names from _gamm4_ models.
* `model_frame()` did not work when the model frame had spline-terms and weights.
* Fix issue in `robust()` when `exponentiate = TRUE` and `conf.int = FALSE`.
* `reliab_test()` returned an error when the provided data frame has less than three columns, instead of returning `NULL`.

# sjstats 0.15.0

## General

* Added new Vignette _Statistics for Bayesian Models_.

## New functions

* `equi_test()` to test if parameter values in Bayesian estimation should be accepted or rejected.
* `mediation()` to print a summary of a mediation analysis from multivariate response models fitted with _brms_.

## Changes to functions

* `link_inverse()` now also returns the link-inverse function for cumulative-family _brms_-models.
* `model_family()` now also returns an `is_ordinal`-element with information if the model is ordinal resp. a cumulative link model.
* Functions that access model information (like `model_family()`) now better support `vglm`-models (package _VGAM_).
* `r2()` now also calculates the standard error for _brms_ or _stanreg_ models.
* `r2()` gets a `loo`-argument to calculate LOO-adjusted rsquared values for _brms_ or _stanreg_ models. This measure comes conceptionally closer to an adjusted r-squared measure.
* Effect sizes (`anova_stats()`, `eta_sq()` etc.) are now also computed for mixed models.
* To avoid confusion, `n_eff()` now computes the number of effective samples, and no longer its ratio in relation to the total number of samples.
* The column name for the ratio of the number of effective samples in `tidy_stan()` is now named *neff_ratio*, to avoid confusion.

## Bug fixes

* Fixed issue in `se()` for `icc()`-objects, where random effect term could not be found.
* Fixed issue in `se()` for `merMod`-objects.
* Fixed issue in `p_value()` for mixed models with KR-approximation, which is now more accurate.

# sjstats 0.14.3

## General

* Remove _tidyverse_ from suggested packages, as requested by maintainers.

## Breaking Changes

* `mwu()` now requires a data frame as first argument, followed by the names of the two variables to perform the Mann-Whitney-U-Test on.

## Changes to functions

* `tidy_stan()` was improved especially for more complex multilevel models.
* Make `tidy_stan()` for large `brmsfit`-objects (esp. with random effects) more efficient.
* Better `print()`-method for `tidy_stan()`, `hdi()`, `rope()`, `icc()` and some other functions.
* `link_inverse()` now also should return the link-inverse function for most (or some or all?) custom families of _brms_-models.
* The `weight.by`-arguments in `grpmean()` and `mwu()` now should be a variable name from a variable in `x`, and no longer a separate vector.

## New functions

* `model_family()` to get model-information about family and link-functions. This function is intended to be "generic" and work with many different model objects, because not all packages provide a `family()` function.

## Bug fixes

* Fix issue with `omega_sq()`, `eta_sq()` etc. when confidence intervals were computed with bootstrapping and the model-formula contained function calls like `scale()` or `as.factor()`.
* Fix issue with `p_value()` for unconditional mixed models.
* Fix typo in `xtab_statistics()`.
* Fix issue with wrong calculation of Nagelkerke's r-squared value in `r2()`.
* Fix issue for factors with character leves in `typical_value()`, when argument `fun` for factors was set to `mode`.
* Don't show prior-samples in `hdi()`, `tidy_stan()` etc. for _brmsfit_-objects.
* Fixed issues in `model_frame()`with spline-terms when missing values were removed due to casewise deletion.

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
