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

# Run driftd() without its console output and return the messages it gave.
driftd_messages <- function(...) {
  msgs <- character()
  withCallingHandlers(
    utils::capture.output(driftd(...)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  msgs
}

# Small data set used by the edge-case tests: three participants, two drift conditions.
edge_data <- function() {
  simulate_ddm(150, a = 1.2, t0 = 0.3, z = 0.5, v = c(easy = 1.5, hard = 0.5), ids = 1:3, seed = 1)
}

# Parameter columns of a fit's individual estimates, without the ID column.
estimate_values <- function(fit) {
  est <- fit$individual_estimates
  est[, names(est) != fit$settings$columns$sub_id]
}

# Run every function that works on a fit, and build every plot, so that a design or model
# specification that breaks any of them fails the test.
expect_downstream_works <- function(fit) {
  pred <- predict(fit, n_sim = 20, seed = 1)
  expect_s3_class(pred, "driftd_prediction")
  expect_s3_class(aggregate(pred), "data.frame")
  expect_s3_class(summary(pred), "data.frame")
  for (type in c("rt", "choice", "scatter", "participants", "qp")) {
    expect_no_error(suppressMessages(ggplot2::ggplot_build(plot(pred, type = type))))
  }
  expect_s3_class(simulate(fit, seed = 1), "data.frame")
  expect_s3_class(estimates(fit), "data.frame")
  expect_s3_class(estimates(fit, format = "cells"), "data.frame")
  expect_output(fit_summary <- summary(fit))
  expect_s3_class(fit_summary, "data.frame")
  expect_output(print(fit))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit))
}
