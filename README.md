
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- <!-- Global site tag (gtag.js) - Google Analytics -->

–>
<!-- <script async src="https://www.googletagmanager.com/gtag/js?id=G-6XKC5E4PWL"></script> -->
<!-- <script> --> <!--   window.dataLayer = window.dataLayer || []; -->
<!--   function gtag(){dataLayer.push(arguments);} -->
<!--   gtag('js', new Date()); -->

<!--   gtag('config', 'G-6XKC5E4PWL'); -->
<!-- </script> -->

# LAWBL: Latent (Variable) Analysis With Bayesian Learning

[![Project Status: Active ? The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/LAWBL)](https://cran.r-project.org/package=LAWBL)
[![](https://cranlogs.r-pkg.org/badges/LAWBL?color=brightgreen)](https://cran.r-project.org/package=LAWBL)
[![](http://cranlogs.r-pkg.org/badges/grand-total/LAWBL?color=green)](https://cran.r-project.org/package=LAWBL)

## How to cite the package

Chen, J. (2022). LAWBL: Latent (variable) analysis with Bayesian
learning (R package version 1.5.0). Retrieved from
<https://CRAN.R-project.org/package=LAWBL>

## Introduction

**LAWBL** represents a partially exploratory-confirmatory approach to
model latent variables based on Bayesian learning. Built on the power of
statistical learning, it can address psychometric challenges such as
parameter specification, local dependence, and factor extraction. Built
on the scalability and flexibility of Bayesian inference and resampling
techniques, it can accommodate modeling frameworks such as factor
analysis, item response theory, cognitive diagnosis modeling and causal
or explanatory modeling. The package can also handle different response
formats or a mix of them, with or without missingness.

## Features

-   Partially CFA (PCFA) for continuous data: regularization of loading
    specification and local dependence; PCFA with local independence
    (PCFA-LI); CFA with local dependence (CFA-LD)
-   Generalized PCFA (GPCFA) for continuous, categorical, or mixed-type
    data, with or without missingness; GPCFA with local independence
    (GPCFA-LI); Generalized CFA with local dependence (GCFA-LD)
-   Partially confirmatory item response model (PCIRM) for continuous
    and dichotomous data with intercept terms; PCIRM-LI; CIRM-LD
-   Bayesian regularized EFA (BREFA): factor extraction and parameter
    estimation in one step; Fully and partially EFA: unknown number of
    factors without or with partial knowledge
-   Estimation using different Bayesian learning methods and MCMC
    algorithms
-   Simulating data based on all aforementioned models
-   Plotting trace, density or Gelman-Rubin diagnostics based on
    eigenvalue
-   Summary of all parameters with both point and interval estimates

Please refer to the [online
tutorials](https://jinsong-chen.github.io/LAWBL/articles/LAWBL.html) for
more details.

## Installation

1.  Install the stable version from [CRAN](https://CRAN.R-project.org)
    with:

``` r
install.packages("LAWBL")
```

2.  Install the `devtools` package (if necessary), and install the
    development version from the Github.

``` r
# install.packages("devtools")
devtools::install_github("Jinsong-Chen/LAWBL")
```

## References

Chen, J. (2020). A partially confirmatory approach to the
multidimensional item response theory with the Bayesian Lasso.
*Psychometrika*. 85(3), 738-774. DOI: 10.1007/s11336-020-09724-3.

Chen, J., Guo, Z., Zhang, L., & Pan, J. (2021). A partially confirmatory
approach to scale development with the Bayesian Lasso. *Psychological
Methods*. 26(2), 210–235. DOI: 10.1037/met0000293.

Chen, J. (2021). A generalized partially confirmatory factor analysis
framework with mixed Bayesian Lasso methods. *Multivariate Behavioral
Research*. DOI: 10.1080/00273171.2021.1925520.

Chen, J. (2021). A Bayesian regularized approach to exploratory factor
analysis in one step. *Structural Equation Modeling: A Multidisciplinary
Journal*. DOI: 10.1080/10705511.2020.1854763.

Chen, J. (2022). Partially confirmatory approach to factor analysis with
Bayesian learning: A LAWBL tutorial. *Structural Equation Modeling: A
Multidisciplinary Journal*. DOI: 10.1080/00273171.2021.1925520.

Chen, J. (In Press). Fully and partially exploratory factor analysis
with bi-level Bayesian regularization. *Behavior Research Methods*.
