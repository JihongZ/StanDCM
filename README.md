
# StanDCM <a href='https://www.jihongzhang.org/rpackages/'><img src='man/figures/logo.png' align="right" height="139" /></a>

[![Travis build
status](https://travis-ci.org/JihongZ/StanDCM.svg?branch=master)](https://travis-ci.org/JihongZ/StanDCM)
[![R build
status](https://github.com/JihongZ/StanDCM/workflows/R-CMD-check/badge.svg)](https://github.com/JihongZ/StanDCM/actions)

## Overview

**StanDCM** is a helpful tool of estimating Diagnostic Classification
Models (DCM) via Stan

## Learning resources

  - Check [Jiang and Carter’s (2019)
    article](https://doi.org/10.3758/s13428-018-1069-9) on Using
    Hamiltonian Monte Carlo to estimate the log-linear cognitive
    diagnosis model via Stan

## Features of the package

  - Estimating log-linear cognitive diagnosis model (LCDM) and a variety
    of widely-used models subsumed by the LCDM, including the DINA
    model, DINO model, additive-CDM (A-CDM), linear logistic model
    (LLM), reduced reparametrized unified model (RRUM),
    multiple-strategy DINA model for dichotomous responses
  - Specifying customized prior distributions for parameters
  - Computing can be achieved in parallel environments
  - Estimating models within the LCDM model framework using
    user-specified design matrix
  - Estimating rating scale DCM for ordinal and nominal responses

## Installation

To install this package from source:

1)  Users may need to install the
    [rstan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
    in order to execute the functions of StanDCM package.

2)  Windows users should avoid using space when installing rstan.

3)  After installing rstan package, users can use the lines beblow to
    install StanDCM package.

<!-- end list -->

``` r
# the development version from GitHub
# install.packages("devtools")
devtools::install_github("JihongZ/StanDCM")
```

The parametric version of DCM R package named GDINA can be found in R
CRAN at [here](https://CRAN.R-project.org/package=GDINA)

## Usage

You can fit a LCDM model:

``` r
mod.LCDM <- StanLCDM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
summary(mod.LCDM)
```

or you can fit a DINA model by simply using:

``` r
mod.DINA <- StanDINA.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
```

Posterior predictive model checking (PPMC) could be also conducted:

``` r
StanDCM.ppmc(mod.LCDM, respMatrix, n.sim = 500, n.burnin = 100)
```
