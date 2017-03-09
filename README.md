sjstats - Collection of Convenient Functions for Common Statistical Computations
------------------------------------------------------------------------------
Collection of convenient functions for common statistical computations, which are not directly provided by R's base or stats packages. 

This package aims at providing, **first**, shortcuts for statistical measures, which otherwise could only be calculated with additional effort (like standard errors or root mean squared errors).

**Second**, these shortcut functions are generic (if appropriate), and can be applied not only to vectors, but also to other objects as well (e.g., the Coefficient of Variation can be computed for vectors, linear models, or linear mixed models; the `r2()`-function returns the r-squared value for _lm_, _glm_, _merMod_ or _lme_ objects).

Most functions of this package are designed as _summary functions_, i.e. they do not transform the input vector; rather, they return a summary, which is sometimes a vector and sometimes a [tidy data frame](https://cran.r-project.org/package=broom/vignettes/broom.html). The focus of most functions lies on summary statistics or fit measures for regression models, including generalized linear models and mixed effects models. However, some of the functions deal with other statistical measures, like Cronbach's Alpha, Cramer's V, Phi etc.

The comprised tools include:

* For regression and mixed models: Coefficient of Variation, Root Mean Squared Error, Residual Standard Error, Coefficient of Discrimination, R-squared and pseudo-R-squared values, standardized beta values
* Especially for mixed models: Design effect, ICC, sample size calculation, convergence and overdispersion tests

Other statistics:

* Cramer's V, Cronbach's Alpha, Mean Inter-Item-Correlation, Mann-Whitney-U-Test, Item-scale reliability tests


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
