sjstats - Statistical Functions for Regression Models
------------------------------------------------------------------------------
**sjstats** aims at providing convenient functions for common statistical computations, which are not provided by R's base or stats packages. This package comprises tools that help computing various statistics for simple regression or mixed models. These include:

* For regression and mixed models: Coefficient of Variation, Root Mean Squared Error, Residual Standard Error, Coefficient of Discrimination, R-squared and pseudo-R-squared values, standardized beta values
* Especially for mixed models: Design effect, ICC, sample size calculation, convergence and overdispersion tests

Some functions are also for simple statistics, like:

* Cramer's V, Cronbach's Alpha, Mean Inter-Item-Correlation, Mann-Whitney-U-Test, Item-scale reliability tests

## Installation

### Latest development build

To install the latest development snapshot (see latest changes below), type following commands into the R console:

```r
library(devtools)
devtools::install_github("sjPlot/sjstats")
```

### Officiale, stable release


## Citation

In case you want / have to cite my package, please use `citation('sjstats')` for citation information. 
