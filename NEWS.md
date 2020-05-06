# sjstats 0.18.0

## General

* Effect size computation functions (like `eta_sq()`) now internally call the related functions from the *effectsize* package.
* Remove packages from "Suggest" that have been removed from CRAN.

# sjstats 0.17.9

## Bug fixes

* Fixed documentation for `chisq_gof()`.
* Fixed issue in `anova_stats()` with incorrect effect sizes for certain Anova types (that included an intercept).

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
* To get consistent argument names across functions, argument `weight.by` is now deprecated and renamed into `weights`.

## Bug fixes

* Fix issues with effect size computation for repeated-measure Anova when using bootstrapping to compute confidence intervals.
* `grpmean()` now also adjusts the `n`-columm for weighted data.
* `icc()`, `re_var()` and  `get_re_var()` now correctly compute the random-effect-variances for models with multiple random slopes per random effect term (e.g., `(1 + rs1 + rs2 | grp)`).
* Fix issues in `tidy_stan()`, `mcse()`, `hdi()` and `n_eff()` for `stan_polr()`-models.
* Plotting `equi_test()` did not work for intercept-only models.
