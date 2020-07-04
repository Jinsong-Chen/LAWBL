
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LAWBL: Latent (variable) Analysis with Bayesian Learning

The long-term goal of **LAWBL** is to provide a analytical framework for
modeling latent variables with different Bayesian learning methods.

## Installation

1)  Install the released version from [CRAN](https://CRAN.R-project.org)
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

Currently, this package includes the Partially Confirmatory Factor
Analysis (PCFA), a partially confirmatory approach covering a wide range
of the exploratory-confirmatory continuum in factor analytic models
(Chen, Guo, Zhang, & Pan, 2020). There are two major model variants with
different constraints for identification. One assumes local independence
(LI) with a more exploratory tendency, which can be also called the
E-step. The other allows local dependence (LD) with a more confirmatory
tendency, which can be also called the C-step. Parameters are obtained
by sampling from the posterior distributions with the Markov chain Monte
Carlo (MCMC) techniques. Different Bayesian Lasso methods are used to
regularize the loading pattern and local dependence.

Although only continuous data are supported currently, inclusion of
mixed-type data is on schedule. More Bayesian learning approaches will
be also included in future releases of this package.

For an example of how to use the package, see vignettes

## How to use this package?

  - To estimate the E-step (when only a few loadings can be specified,
    e.g., 2 per factor), use **m \<- pcfa(dat=dat,Q=Q,LD=F)**
  - To estimate the C-step (with one specified loading per item), use
    **m \<- pcfa(dat=dat,Q=Q,LD=T)**
  - To summarize basic information after estimation, use **summary(m)**
  - To summarize significant loadings in pattern/Q-matrix format, use
    **summary(m,what=‘qlambda’)**
  - To summarize factorial eigenvalues, use **summary(m,what=‘eigen’)**
  - To summarize significant LD terms, use **summary(m,what=‘offpsx’)**
  - To plot eigenvalues’ trace, use **plot\_eigen(m)**
  - To plot eigenvalues’ density, use **plot\_eigen(m, what=‘density’)**
  - To plot eigenvalues’ adjusted PSRF, use **plot\_eigen(m,
    what=‘APSR’)**

## Reference

Chen, J., Guo, Z., Zhang, L., & Pan, J. (In Press). A partially
confirmatory approach to scale development with the Bayesian Lasso.
*Psychological Methods*.
