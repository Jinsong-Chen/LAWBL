
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LAWBL: Latent (variable) Analysis with Bayesian Learning

[![Project Status: Active ? The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/LAWBL)](https://cran.r-project.org/package=LAWBL)

## Installation

1)  Install the stable version from [CRAN](https://CRAN.R-project.org)
    with:

<!-- end list -->

``` r
install.packages("LAWBL")
```

2)  Install the `devtools` package (if necessary), and install the
    development version from the Github.

<!-- end list -->

``` r
# install.packages("devtools")
devtools::install_github("Jinsong-Chen/LAWBL")
```

## What can this package do?

The long-term goal of **LAWBL** is to provide an analytical framework
for modeling latent variables with different Bayesian learning methods.

Currently, this package includes the Partially Confirmatory Factor
Analysis (PCFA), a partially confirmatory approach covering a wide range
of the exploratory-confirmatory continuum in factor analytic models. The
PCFA (Chen, Guo, Zhang, & Pan, 2020) is only for continuous data, while
the generalized PCFA (GPCFA) covers both continuous, categorical, or
mixed-type data. There are two major model variants with different
constraints for identification. One assumes local independence (LI) with
a more exploratory tendency, which can be also called the E-step. The
other allows local dependence (LD) with a more confirmatory tendency,
which can be also called the C-step. Parameters are obtained by sampling
from the posterior distributions with the Markov chain Monte Carlo
(MCMC) techniques. Different Bayesian Lasso methods are used to
regularize the loading pattern and local dependence.

For examples of how to use the package, see vignettes or
[here](https://jinsong-chen.github.io/LAWBL/articles/pcfa-examples.html)
for PCFA with continuous data;
[here](https://jinsong-chen.github.io/LAWBL/articles/gpcfa-examples.html)
for GPCFA with categorical data.

## How to use this package?

  - To estimate the E-step (when only a few loadings can be specified,
    e.g., 2 per factor), use *m \<- pcfa(dat=dat,Q=Q,LD=F)*
  - To estimate the C-step (with one specified loading per item), use *m
    \<- pcfa(dat=dat,Q=Q,LD=T)*
  - To summarize basic information after estimation, use *summary(m)*
  - To summarize significant loadings in pattern/Q-matrix format, use
    *summary(m,what=‘qlambda’)*
  - To summarize factorial eigenvalues, use *summary(m,what=‘eigen’)*
  - To summarize significant LD terms, use *summary(m,what=‘offpsx’)*
  - To plot eigenvalues’ trace, use *plot\_eigen(m)*
  - To plot eigenvalues’ density, use *plot\_eigen(m, what=‘density’)*
  - To plot eigenvalues’ adjusted PSRF, use *plot\_eigen(m,
    what=‘APSR’)*

## Reference

Chen, J., Guo, Z., Zhang, L., & Pan, J. (2020). A partially confirmatory
approach to scale development with the Bayesian Lasso. *Psychological
Methods*. Advance online publication.
<http://dx.doi.org/10.1037/met0000293>
