\encoding{UTF-8}
\name{sjstats-package}
\alias{sjstats-package}
\alias{sjstats}
\docType{package}
\title{Collection of Convenient Functions for Common Statistical Computations}
\description{Collection of convenient functions for common statistical computations, which are not directly provided by R's base or stats packages.

This package aims at providing, first, shortcuts for statistical measures, which otherwise could only be calculated with additional effort (like standard errors or root mean squared errors).

Second, these shortcut functions are generic (if appropriate), and can be applied not only to vectors, but also to other objects as well (e.g., the Coefficient of Variation can be computed for vectors, linear models, or linear mixed models; the \code{r2()}-function returns the r-squared value for \code{lm}, \code{glm}, \code{merMod}, \code{glmmTMB}, or \code{lme} and other objects).

Most functions of this package are designed as \emph{summary functions}, i.e. they do not transform the input vector; rather, they return a summary, which is sometimes a vector and sometimes a \href{https://cran.r-project.org/package=broom/vignettes/broom.html}{tidy data frame}. he focus of most functions lies on summary statistics or fit measures for regression models, including generalized linear models, mixed effects models or Bayesian models. However, some of the functions deal with other statistical measures, like Cronbach's Alpha, Cramer's V, Phi etc.

The comprised tools include:
\itemize{
  \item For regression and mixed models: Coefficient of Variation, Root Mean Squared Error, Residual Standard Error, Coefficient of Discrimination, R-squared and pseudo-R-squared values, standardized beta values
  \item Especially for mixed models: Design effect, ICC, sample size calculation and convergence tests
  \item Especially for Bayesian models: Highest Density Interval, region of practical equivalence (rope), Monte Carlo Standard Errors, ratio of number of effective samples, mediation analysis, Test for Practical Equivalence
  \item Fit and accuracy measures for regression models: Overdispersion tests, accuracy of predictions, test/training-error comparisons, error rate and binned residual plots for logistic regression models
  \item For anova-tables: Eta-squared, Partial Eta-squared, Omega-squared and Partial Omega-squared statistics
}

Furthermore, \strong{sjstats} has functions to access information from model objects, which either support more model objects than their \strong{stats} counterparts, or provide easy access to model attributes, like:

\itemize{
  \item \code{model_frame()} to get the model frame
  \item \code{model_family()} to get information about the model family, link functions etc.
  \item \code{link_inverse()} to get the link-inverse function
  \item \code{pred_vars()} and \code{resp_var()} to get the names of either the dependent or independent variables, or
  \item \code{var_names()} to get the "cleaned" variables names from a model object (cleaned means, things like \code{s()} or \code{log()} are removed from the returned character vector with variable names.)
}

Other statistics:
\itemize{
  \item Cramer's V, Cronbach's Alpha, Mean Inter-Item-Correlation, Mann-Whitney-U-Test, Item-scale reliability tests
}
}
