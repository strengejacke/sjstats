# sjstats - Collection of Convenient Functions for Common Statistical Computations <img src="man/figures/logo.png" align="right" />

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1284472.svg)](https://doi.org/10.5281/zenodo.1284472)

Collection of convenient functions for common statistical computations, which are not directly provided by R's base or stats packages. 

This package aims at providing, **first**, shortcuts for statistical measures, which otherwise could only be calculated with additional effort (like Cramer's V, Phi, or effict size statistics like Eta or Omega squared), or for which currently no functions are available.

**Second**, another focus lies on implementations of common statistical significance tests with a consistent syntax, like t-test, Mann-Whitney test, Chi-squared test, and more. These functions are designed to be more user-friendly and also support weights, i.e. weighted statistics can be calculated.

**Finally**, the package includes miscellaneous functions that are either not yet available in R (like `svyglm.nb()` or `svyglm.zip()` to calculate negative binomial or zero-inflated poisson models for survey data) or are just convenient for daily work (like functions for bootstrapping, or ANOVA summary tables).

The comprised tools include:

* Especially for mixed models: design effect, sample size calculation
* Significance tests: Correlation, Chi-squared test, t-test, Mann-Whitney-U test, Wilcoxon rank sum test, Kruskal-Wallis test.

Note that most functions that formerly were available in this package have been moved to the [**easystats** project](http://easystats.github.io/easystats).

## Documentation

Please visit [https://strengejacke.github.io/sjstats/](https://strengejacke.github.io/sjstats/) for documentation and vignettes.

## Installation

### Latest development build

To install the latest development snapshot (see latest changes below), type following commands into the R console:

```r
library(remotes)
remotes::install_github("strengejacke/sjstats")
```
### Officiale, stable release

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sjstats)](https://cran.r-project.org/package=sjstats)

To install the latest stable release from CRAN, type following command into the R console:

```r
install.packages("sjstats")
```

## Citation

In case you want / have to cite my package, please use `citation('sjstats')` for citation information. 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1284472.svg)](https://doi.org/10.5281/zenodo.1284472)
