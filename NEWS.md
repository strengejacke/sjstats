# sjstats 0.17.8

## Deprecated and defunct

_sjstats_ is being re-structured, and many functions are re-implemented in new packages that are part of a new project called **easystats**.

Therefore, following functions are now deprecated:

* `cohens_f()`, please use `effectsize::cohens_f()`.
* `std_beta()`, please use `effectsize::standardize_parameters()`.
* `tidy_stan()`, please use `parameters::model_parameters()`.
* `scale_weights()`, please use `parameters::rescale_weights()`.
* `robust()`, please use `parameters::standard_error_robust()`.

## General

* Functions for weighted statistics with prefix `wtd_*()` have been renamed to `weighted_*()`.
* `svy_md()` was renamed to `survey_median()`.
* `mannwhitney()` is an alias for `mwu()`.
* `means_by_group()` is an alias for `grpmean()`.

# sjstats 0.17.7

## Deprecated and defunct

_sjstats_ is being re-structured, and many functions are re-implemented in new packages that are part of a new project called **easystats**. The aim of **easystats** is to provide a unifying and consistent framework to tame, discipline and harness the scary R statistics and their pesky models.

Therefore, following functions are now deprecated:

* `p_value()`, please use `parameters::p_value()`
* `se()`, please use `parameters::standard_error()`

## General

* Revise some functions to cope with the forthcoming _insight_ update.

# sjstats 0.17.6

## General

* Minor revisions to meet the changes in the forthcoming update from *tidyr*.
* `design_effect()` is an alias for `deff()`.
* `samplesize_mixed()` is an alias for `smpsize_lmm()`.
* `crosstable_statistics()` is an alias for `xtab_statistics()`.

## New functions

* `svyglm.zip()` to fit zero-inflated Poisson models for survey-designs.

## Changes to functions

* `phi()` and `cramer()` can now compute confidence intervals.
* `tidy_stan()` removes prior parameters from output.
* `tidy_stan()` now also prints the probability of direction.

## Bug fixes

* Fix bug with wrong computation in `odds_to_rr()`.

# sjstats 0.17.5

## New functions

* `epsilon_sq()`, to compute epsilon-squared effect-size.

## Deprecated and defunct

_sjstats_ is being re-structured, and many functions are re-implemented in new packages that are part of a new project called **easystats**. The aim of **easystats** is to provide a unifying and consistent framework to tame, discipline and harness the scary R statistics and their pesky models.

Therefore, following functions are now deprecated:

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
* `predictive_accurarcy()`, please use `performance::performance_accuracy()`
* `cronb()`, please use `performance::cronbachs_alpha()`
* `difficulty()`, please use `performance::item_difficulty()`
* `mic()`, please use `performance::item_intercor()`
* `pca()`, please use `parameters::principal_components()`
* `pca_rotate()`, please use `parameters::principal_components()`
* `r2()`, please use `performance::r2()`
* `icc()`, please use `performance::icc()`
* `rmse()`, please use `performance::rmse()`
* `rse()`, please use `performance::rse()`
* `mse()`, please use `performance::mse()`
* `hdi()`, please use `bayestestR::hdi()`
* `cred_int()`, please use `bayestestR::ci()`
* `rope()`, please use `bayestestR::rope()`
* `n_eff()`, please use `bayestestR::effective_sample()`
* `equi_test()`, please use `bayestestR::equivalence_test()`
* `multicollin()`, please use `performance::check_collinearity()`
* `normality()`, please use `performance::check_normality()`
* `autocorrelation()`, please use `performance::check_autocorrelation()`
* `heteroskedastic()`, please use `performance::check_heteroscedasticity()`
* `outliers()`, please use `performance::check_outliers()`

## Changes to functions

* Anova-stats functions (like `eta_sq()`) get a `method`-argument to define the method for computing confidence intervals from bootstrapping.

## Bug fixes

* In some situations, `smpsize_lmm()` could result in negative sample-size recommendations. This was fixed, and a warning is now shown indicating that the parameters for the power-calculation should be modified.
* Fixed issue with wrong calculated effect size `r` in `mwu()` if group-factor contained more than two groups.

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
