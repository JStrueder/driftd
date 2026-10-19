# driftd 0.1.0

First public release.

* `driftd()` fits the diffusion decision model to each participant by maximum likelihood, with densities from rtdists. You specify which parameters vary across which conditions; starting values come from EZ-diffusion, and alternative starts are tried when a fit fails or does not replicate.
* Wald standard errors and covariance matrices from the numerical Hessian.
* Analysis-ready parameter tables, and model fit checks that compare observed and predicted response probabilities and RT quantiles, with prediction intervals.
* Model comparison with `compare_models()`.
* Parameter and model recovery, and simulation from the fitted model with `simulate()`.
* Optional outlier mixture (`p_outlier`, off by default), as in HDDM: a fixed share of trials is modelled as uniform contaminants, which protects the estimates against fast guesses and lapses.
* Example data `speed_acc` (Wagenmakers et al., 2008, Experiment 1, as distributed with rtdists).
