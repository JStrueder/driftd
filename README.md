
<!-- README.md is generated from README.Rmd. Please edit that file -->

# driftd

<!-- badges: start -->
<!-- badges: end -->

driftd provides a streamlined workflow for fitting Drift‑Diffusion
Models (DDMs) to trial‑level RT and binary response data, performing
participant‑level maximum‑likelihood fits (via `nlminb`/`rtdists`) and
simulation‑based posterior predictive checks (PPCs).

## Installation

From CRAN (when available):

``` r
install.packages("driftd")
```

Development version from GitHub:

``` r
devtools::install_github("JStrueder/driftd")
```

## Quick Example: Fitting a DDM and Running a Posterior Predictive Check

This example demonstrates the typical workflow. Fitting can be slow, so
the fitting step is shown with `eval = FALSE`. Use the included example
dataset to follow along.

``` r
library(driftd)

# load example data included with the package
data(example_data)

# Quick inspect
head(example_data)

# Fit the model (this may take time; run interactively)
results <- driftd(
  data = example_data,
  sub_id = "ID",
  response = "R",
  rt = "RT",
  vary_z = c("val"),
  vary_v = c("prob", "val"),
  fix_st0 = TRUE,
  fix_sz = TRUE,
  fix_sv = TRUE,
  print_starting_vals = FALSE
)

# Examine results
results$individual_estimates
results$standard_errors   # asymptotic standard errors (see ?driftd)
results$summary

# Posterior predictive check (simulate 1 dataset)
ppc <- posterior_predict(results, n_sim = 1)

# Plot PPC summary
plot(ppc, type = "summary")

# Plot the comparison of observed vs. simulated statistics for each participant
plot(ppc, type = "participant")

# Plot the comparison of the full RT distributions for a specific condition
plot(ppc, type = "density", condition = "pos_0.55")
```

## License and citation

This package is licensed under the MIT license (see `LICENSE`). Please
cite the package when used in research:

``` r
citation("driftd")
```

## Contact and issues

- Maintainer: J. D. Strueder (<jstrueder@iese.edu>)
