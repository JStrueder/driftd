# Shared helpers for the driftd tests. testthat loads files starting with "helper"
# before running the tests.

# Run driftd() without its console output (progress bar, summary table, messages).
quiet_driftd <- function(...) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(driftd(...)))
  out
}

# Run driftd() without its console output, but let messages through (to test them).
driftd_no_output <- function(...) {
  out <- NULL
  utils::capture.output(out <- driftd(...))
  out
}

# Expect every element of `actual` to be within `tol` (absolute) of `expected`.
expect_near <- function(actual, expected, tol) {
  diff <- max(abs(unname(unlist(actual)) - unname(unlist(expected))))
  expect_true(diff <= tol, label = sprintf("max absolute difference %.3g <= %g", diff, tol))
}

# Simulate trial-level data from the diffusion model with rtdists::rdiffusion().
# `v` is a named vector with one drift rate per condition; z is relative (0-1).
simulate_ddm <- function(n_per_cond, a, t0, z, v, sv = 0, ids = 1, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(ids, function(id) {
    do.call(rbind, lapply(names(v), function(cond) {
      s <- rtdists::rdiffusion(n_per_cond, a = a, v = v[[cond]], t0 = t0, z = z * a, sv = sv)
      data.frame(id = id, rt = s$rt, resp = as.integer(s$response == "upper"), cond = cond)
    }))
  }))
}
