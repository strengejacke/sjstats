# sjstats - Collection of Convenient Functions for Common Statistical Computations <img src="man/figures/logo.png" align="right" />

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1284472.svg)](https://doi.org/10.5281/zenodo.1284472)

Collection of convenient functions for common statistical computations, which are not directly provided by R's base or stats packages. 

This package aims at providing, **first**, shortcuts for statistical measures, which otherwise could only be calculated with additional effort (like Cramer's V, Phi, or effict size statistics like Eta or Omega squared), or for which currently no functions available.

**Second**, another focus lies on weighted variants of common statistical measures and tests like weighted standard error, mean, t-test, correlation, and more.

The comprised tools include:

* Especially for mixed models: design effect, sample size calculation
* Especially for Bayesian models: mediation analysis
* For anova-tables: Eta-squared, Partial Eta-squared, Omega-squared, Partial Omega-squared and Epsilon-squared statistics
* Weighted statistics and tests for: mean, median, standard error, standard deviation, correlation, Chi-squared test, t-test, Mann-Whitney-U-test

## Documentation

Please visit [https://strengejacke.github.io/sjstats/](https://strengejacke.github.io/sjstats/) for documentation and vignettes.

## Installation

### Latest development build

To install the latest development snapshot (see latest changes below), type following commands into the R console:

```r
library(devtools)
devtools::install_github("strengejacke/sjstats")
```
### Officiale, stable release

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sjstats)](https://cran.r-project.org/package=sjstats)
&#160;&#160;
[![downloads](http://cranlogs.r-pkg.org/badges/sjstats)](http://cranlogs.r-pkg.org/)
&#160;&#160;
[![total](http://cranlogs.r-pkg.org/badges/grand-total/sjstats)](http://cranlogs.r-pkg.org/)

To install the latest stable release from CRAN, type following command into the R console:

```r
install.packages("sjstats")
```

## Citation

In case you want / have to cite my package, please use `citation('sjstats')` for citation information. 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1284472.svg)](https://doi.org/10.5281/zenodo.1284472)
