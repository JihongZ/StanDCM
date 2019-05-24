
<!-- README.md is generated from README.Rmd. Please edit that file -->

# StanDCM Package for Using Stan to Estimate Diagnostic Classification Models



## Learning resources

  - Check [Jiang and Carter’s (2019) article](https://doi.org/10.3758/s13428-018-1069-9) on Using Hamiltonian Monte Carlo to estimate the log-linear cognitive diagnosis model via Stan

## Features of the package

  - Estimating log-linear cognitive diagnosis model (LCDM) and a variety of widely-used models subsumed
    by the LCDM, including the DINA model, DINO model, additive-CDM (A-CDM), linear logistic model (LLM), reduced reparametrized unified model (RRUM), multiple-strategy DINA model
    for dichotomous responses
  - Specifying customized prior distributions for parameters
  - Computing can be achieved in parallel environments
  - Estimating models within the LCDM model framework using
    user-specified design matrix
  - Estimating rating scale DCM for ordinal and nominal responses

## Installation

To install this package from source:

1)  Users may need to install the
    [rstan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) in order to execute the functions of StanDCM package.

2)  Windows users should avoid using space when installing rstan.

3)  After installing rstan package, users can use the lines beblow to install StanDCM package.

<!-- end list -->

``` r
# install.packages("devtools")
devtools::install_github("JihongZ/StanDCM")
```

The parametric version of DCM R package named GDINA can be found in R CRAN at
[here](https://CRAN.R-project.org/package=GDINA)
