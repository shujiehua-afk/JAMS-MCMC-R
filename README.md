# Adaptive MCMC for Multimodal Targets (R Implementation)

## Overview

This repository contains an R implementation of several MCMC components designed for sampling from multimodal target distributions.
The codebase includes:

* A core JAMS-style adaptive MCMC framework
* Tools for local moves, jump moves, empirical covariance adaptation
* Demonstration notebooks for 2D, high-dimensional Gaussian mixtures, banana-shaped distributions, and parallel tempering experiments

This repository is intended both as a research record and a growing personal portfolio.

---

## Repository Structure

```
/src
    MCMC_func.R               # Main JAMS-style MCMC implementation
    mcmc2.R                   # Additional tools: banana-shaped distribution, PT utilities
/examples
    MCMC_real_example.Rmd     # 2D/3D multimodal demo, ESS evaluation, diagnostics
    mcmc2.Rmd                 # Banana-shaped + PT demonstrations

README.md
```

---

## Files Description

### `/src/MCMC_func.R`

Core implementation including:

* Local proposal kernel
* Jump proposal kernel
* Empirical covariance update
* Adaptive parameter update
* Deterministic jump mechanism
* Burn-in components

### `/src/mcmc2.R`

Non-JAMS utilities:

* Banana-shaped target demos
* Parallel tempering experiments
* Adaptive covariance experiments

### `/examples/MCMC_real_example.Rmd`

Reproducible R Markdown demonstration:

* 2D multimodal examples
* ESS evaluation
* Covariance / weight adaptation checks

### `/examples/mcmc2.Rmd`

* Banana-shaped distributions
* High-dimensional experiments
* Parallel tempering verification

---

## Usage

### 1. Install Dependencies

```r
install.packages(c("mvtnorm", "ggplot2", "MASS"))
```

### 2. Run Core MCMC

```r
source("src/MCMC_func.R")
out <- JAMS_MCMC(target, init, n_iter = 50000)
```

### 3. Reproduce Examples

Open the files in `/examples` and knit them.

---

## Goals of This Repository

1. Maintain a clean, incremental research log on multimodal adaptive MCMC.
2. Provide reproducible experiments for multimodal sampling challenges.
3. Serve as a personal technical portfolio for MCMC algorithm development.
4. Gradually expand to more general adaptive and tempered MCMC methods.

---

## Future Work

* Add unit tests (testthat)
* Add visualization utilities for mode transitions
* Provide vignettes for JAMS core logic
* Extend to higher-dimensional PT and AIS
