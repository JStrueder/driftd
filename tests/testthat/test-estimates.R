d <- subset(example_data, ID %in% c(1, 2, 3))
fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_z = "val", vary_v = c("prob", "val"))
n_par <- 17   # a, t0, 2 z, 10 v, st0, sz, sv

test_that("long format has one row per participant and parameter", {
  est <- estimates(fit)
  expect_equal(names(est), c("ID", "parameter", "val", "prob", "name", "estimate", "se", "fixed", "converge"))
  expect_equal(nrow(est), 3 * n_par)
  expect_equal(sort(unique(est$parameter)), c("a", "st0", "sv", "sz", "t0", "v", "z"))
  # condition values come from the condition columns, with their original types
  expect_type(est$prob, "double")
  row <- est[est$ID == 2 & est$name == "v_0.45_pos", ]
  expect_equal(row$parameter, "v")
  expect_equal(row$val, "pos")
  expect_equal(row$prob, 0.45)
  expect_equal(row$estimate, fit$individual_estimates$v_0.45_pos[2])
  expect_equal(row$se, fit$standard_errors$v_0.45_pos[2])
  # parameters that do not vary across a column have NA there
  expect_true(all(is.na(est$val[est$parameter %in% c("a", "t0")])))
  expect_true(all(is.na(est$prob[est$parameter == "z"])))
  # fixed parameters are marked and have no SE
  expect_equal(unique(est$parameter[est$fixed]), c("st0", "sz", "sv"))
  expect_true(all(is.na(est$se[est$fixed])))
  expect_equal(nrow(estimates(fit, include_fixed = FALSE)), 3 * (n_par - 3))
})

test_that("cells format has one row per participant and condition cell", {
  cells <- estimates(fit, format = "cells")
  expect_equal(names(cells), c("ID", "val", "prob", "n_trials", "a", "t0", "z", "v", "st0", "sz", "sv", "converge"))
  expect_equal(nrow(cells), 3 * 10)
  row <- cells[cells$ID == 3 & cells$val == "neg" & cells$prob == 0.55, ]
  e <- fit$individual_estimates[fit$individual_estimates$ID == 3, ]
  expect_equal(row$v, e$v_0.55_neg)
  expect_equal(row$z, e$z_neg)
  expect_equal(row$a, e$a)
  expect_equal(row$n_trials, sum(d$ID == 3 & d$val == "neg" & d$prob == 0.55))
  # sorted by participant, then by condition values
  expect_equal(cells$val[1:5], rep("neg", 5))
  expect_equal(cells$prob[1:5], c(0.4, 0.45, 0.5, 0.55, 0.6))
})

test_that("the cells format can be analysed directly, e.g. with a repeated-measures ANOVA", {
  cells <- estimates(fit, format = "cells")
  cells$ID <- factor(cells$ID)
  cells$prob <- factor(cells$prob)
  res <- summary(stats::aov(v ~ prob * val + Error(ID / (prob * val)), data = cells))
  expect_equal(res[["Error: ID:prob"]][[1]]["prob", "Df"], 4)
  expect_equal(res[["Error: ID:val"]][[1]]["val", "Df"], 1)
})

test_that("condition values are correct for levels containing '_' and for designs without conditions", {
  d_us <- subset(example_data, ID == 1)
  d_us$val <- paste0(d_us$val, "_cue")
  fit_us <- quiet_driftd(d_us, sub_id = "ID", response = "R", rt = "RT", vary_v = "val", compute_se = FALSE)
  est <- estimates(fit_us)
  expect_equal(est$val[est$parameter == "v"], c("neg_cue", "pos_cue"))
  expect_equal(estimates(fit_us, "cells")$v, est$estimate[est$parameter == "v"])

  fit0 <- quiet_driftd(subset(example_data, ID %in% c(1, 2)), sub_id = "ID", response = "R", rt = "RT",
                       compute_se = FALSE)
  cells0 <- estimates(fit0, format = "cells")
  expect_equal(nrow(cells0), 2)
  expect_equal(names(cells0), c("ID", "n_trials", "a", "t0", "z", "v", "st0", "sz", "sv", "converge"))
})

test_that("the summary table describes parameters only", {
  expect_equal(names(fit$summary), c("Parameter", "Mean", "SD", "Fixed"))
  expect_false(any(c("converge", "LLE", "AIC", "BIC") %in% fit$summary$Parameter))
  expect_equal(fit$summary$Parameter[fit$summary$Fixed], c("st0", "sz", "sv"))
})

test_that("results from older driftd versions are rejected clearly", {
  old <- fit
  old$settings$parameters <- NULL
  expect_error(estimates(old), "older version of driftd")
})
