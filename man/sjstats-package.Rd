\encoding{UTF-8}
\name{sjstats-package}
\alias{sjstats-package}
\alias{sjstats}
\docType{package}
\title{Collection of Convenient Functions for Common Statistical Computations}
\description{Collection of convenient functions for common statistical computations, which are not directly provided by R's base or stats packages. This package aims at providing, first, shortcuts for statistical measures, which otherwise could only be calculated with additional effort (like standard errors or root mean squared errors). Second, these shortcut functions are generic (if appropriate), and can be applied not only to vectors, but also to other objects as well (e.g., the Coefficient of Variation can be computed for vectors, linear models, or linear mixed models; the \code{r2}-function returns the r-squared value for lm, glm, merMod or lme objects). The focus of most functions lies on summary statistics or fit measures for regression models, including generalized linear models and mixed effects models. However, some of the functions deal with other statistical measures, like Cronbach's Alpha, Cramer's V, Phi etc.

The comprised tools include:
\itemize{
  \item For regression and mixed models: Coefficient of Variation, Root Mean Squared Error, Residual Standard Error, Coefficient of Discrimination, R-squared and pseudo-R-squared values, standardized beta values
  \item Especially for mixed models: Design effect, ICC, sample size calculation, convergence and overdispersion tests
}
Other statistics:
\itemize{
  \item Cramer's V, Cronbach's Alpha, Mean Inter-Item-Correlation, Mann-Whitney-U-Test, Item-scale reliability tests
}
}
