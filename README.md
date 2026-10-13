
<!-- README.md is generated from README.Rmd. Please edit that file -->

# driftd

<!-- badges: start -->
[![R-CMD-check](https://github.com/JStrueder/driftd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JStrueder/driftd/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/JStrueder/driftd/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/JStrueder/driftd/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

driftd is a wrapper around the
[rtdists](https://CRAN.R-project.org/package=rtdists) package (Singmann,
Brown, Gretton, & Heathcote) that makes it easier to fit the diffusion
decision model to response time data and to analyze the results. You
specify which parameters vary across which conditions; driftd builds the
likelihood, fits it to each participant by maximum likelihood, and
provides standard errors, analysis-ready parameter tables, model fit
checks, model comparison, parameter and model recovery, and simulation.
All diffusion model computations are done by rtdists.

## Installation

From CRAN (when available):

``` r
install.packages("driftd")
```

Development version from GitHub:

``` r
devtools::install_github("JStrueder/driftd")
```

To also install the vignette, a complete worked analysis from the data
to reporting:

``` r
devtools::install_github("JStrueder/driftd", build_vignettes = TRUE)
vignette("driftd", package = "driftd")
```

## Quick example

The example data, `speed_acc`, are from a lexical decision experiment
with speed and accuracy instructions (Wagenmakers et al., 2008), as
distributed with rtdists. Fitting all 17 participants takes about 15
minutes per model, so the code is not run here; the vignette shows the
output.

``` r
library(driftd)
data(speed_acc)

# boundary separation by instruction, drift rate by stimulus type
fit <- driftd(speed_acc, sub_id = "id", response = "response", rt = "rt",
              vary_a = "condition", vary_v = "stimulus", fix_st0 = FALSE)
fit$summary
estimates(fit)                     # one row per participant and parameter, with SEs
estimates(fit, format = "cells")   # one row per participant and condition

# model fit
pred <- predict(fit, seed = 1)
plot(pred)                         # RT quantiles
plot(pred, type = "choice")        # response probabilities
summary(pred)                      # participants with the most misfit

# model comparison: does the instruction also affect the drift rate?
fit_v <- driftd(speed_acc, sub_id = "id", response = "response", rt = "rt",
                vary_a = "condition", vary_v = c("stimulus", "condition"), fix_st0 = FALSE)
cmp <- compare_models(a = fit, a_v = fit_v)
plot(cmp)

# parameter and model recovery
plot(recover_parameters(fit, seed = 1))
plot(recover_models(a = fit, a_v = fit_v, seed = 1))

# simulate data from the fitted model
sims <- simulate(fit, nsim = 10, seed = 1)
```

## Credits

driftd is built on [rtdists](https://CRAN.R-project.org/package=rtdists)
(Singmann, Brown, Gretton, & Heathcote), which computes all diffusion
model densities, distribution functions and random numbers, based on the
fast-dm code of Andreas and Jochen Voss. If you use driftd, please cite
rtdists as well; `citation("driftd")` gives both references.

The `speed_acc` example data were collected by Wagenmakers, Ratcliff,
Gomez, and McKoon (2008, Experiment 1) and are taken from the rtdists
package.

Wagenmakers, E.-J., Ratcliff, R., Gomez, P., & McKoon, G. (2008). A
diffusion model account of criterion shifts in the lexical decision
task. *Journal of Memory and Language, 58*(1), 140-159.
<https://doi.org/10.1016/j.jml.2007.04.006>

## License and citation

driftd is licensed under the GNU General Public License, version 3 or
later (see `LICENSE.md`). Please cite the package, and rtdists, when
used in research:

``` r
citation("driftd")
```

## Contact and issues

- Maintainer: J. D. Strueder (<jstrueder@iese.edu>)
