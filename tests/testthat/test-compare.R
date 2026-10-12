# Data simulated with a drift rate that differs between conditions
d_effect <- simulate_ddm(300, a = 1.4, t0 = 0.3, z = 0.5, v = c(easy = 2.0, hard = 0.6), ids = 1:4, seed = 5)
fit_const <- quiet_driftd(d_effect, sub_id = "id", response = "resp", rt = "rt", compute_se = FALSE)
fit_vary <- quiet_driftd(d_effect, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond", compute_se = FALSE)
fit_vary_a <- quiet_driftd(d_effect, sub_id = "id", response = "resp", rt = "rt", vary_a = "cond", compute_se = FALSE)

test_that("the model that generated the data is preferred", {
  cmp <- compare_models(constant = fit_const, varying = fit_vary)
  expect_s3_class(cmp, "driftd_comparison")
  expect_equal(cmp$table$model[which.min(cmp$table$BIC)], "varying")
  expect_equal(cmp$table$model[which.min(cmp$table$AIC)], "varying")
  expect_equal(cmp$table$n_best_BIC, c(0, 4))
  expect_lt(cmp$lrt$p, 1e-6)
  expect_equal(cmp$lrt$n_significant, 4)
})

test_that("information criteria are computed per participant and summed", {
  cmp <- compare_models(fit_const, fit_vary)
  per <- cmp$participants
  est <- fit_vary$individual_estimates
  n <- as.vector(table(d_effect$id))
  k <- 5   # a, t0, z, v_easy, v_hard
  rows <- per[per$model == "fit_vary", ]
  expect_equal(rows$k, rep(k, 4))
  expect_equal(rows$n_trials, n)
  expect_equal(rows$logLik, est$LLE)
  expect_equal(rows$AIC, 2 * k - 2 * est$LLE)
  expect_equal(rows$BIC, k * log(n) - 2 * est$LLE)
  # the BIC stored by driftd() agrees
  expect_equal(rows$BIC, est$BIC)
  tab <- cmp$table
  expect_equal(tab$logLik[2], sum(est$LLE))
  expect_equal(tab$BIC[2], sum(rows$BIC))
  expect_equal(sum(tab$w_AIC), 1)
  expect_equal(sum(tab$w_BIC), 1)
  expect_equal(sum(tab$n_best_BIC), 4)
  expect_true(all(per$dBIC >= 0))
})

test_that("likelihood-ratio tests are run for nested models only", {
  cmp <- compare_models(constant = fit_const, varying = fit_vary, boundary = fit_vary_a)
  # constant is nested in both other models; varying and boundary are not nested in each other
  expect_equal(nrow(cmp$lrt), 2)
  expect_equal(cmp$lrt$restricted, c("constant", "constant"))
  expect_setequal(cmp$lrt$full, c("varying", "boundary"))
  row <- cmp$lrt[cmp$lrt$full == "varying", ]
  lr <- 2 * (fit_vary$individual_estimates$LLE - fit_const$individual_estimates$LLE)
  expect_equal(row$chisq, sum(lr))
  expect_equal(row$df, 4)
  expect_equal(row$p, stats::pchisq(sum(lr), 4, lower.tail = FALSE))
  expect_false(row$boundary)
  expect_null(compare_models(fit_vary, fit_vary_a)$lrt)
})

test_that("nesting is detected from the model specifications", {
  fit_sv <- quiet_driftd(subset(d_effect, id == 1), sub_id = "id", response = "resp", rt = "rt",
                         vary_v = "cond", fix_sv = FALSE, compute_se = FALSE, max_retries = 0)
  fit_v1 <- quiet_driftd(subset(d_effect, id == 1), sub_id = "id", response = "resp", rt = "rt",
                         vary_v = "cond", compute_se = FALSE)
  fit_z <- quiet_driftd(subset(d_effect, id == 1), sub_id = "id", response = "resp", rt = "rt",
                        vary_v = "cond", fix_z = TRUE, z_fixed = 0.4, compute_se = FALSE)
  fit_z2 <- quiet_driftd(subset(d_effect, id == 1), sub_id = "id", response = "resp", rt = "rt",
                         vary_v = "cond", fix_z = TRUE, z_fixed = 0.6, compute_se = FALSE)
  expect_true(driftd:::nesting(fit_v1, fit_sv))          # sv fixed at 0 vs estimated: boundary
  expect_false(driftd:::nesting(fit_z, fit_v1))          # z fixed at 0.4 vs estimated
  expect_true(is.na(driftd:::nesting(fit_sv, fit_v1)))   # more complex model is not nested
  expect_true(is.na(driftd:::nesting(fit_z, fit_z2)))    # fixed at different values
  expect_true(is.na(driftd:::nesting(fit_vary, fit_vary_a)))
})

test_that("an optimization problem in the full model is flagged", {
  broken <- fit_vary
  broken$individual_estimates$LLE[2] <- fit_const$individual_estimates$LLE[2] - 5
  cmp <- compare_models(fit_const, broken)
  expect_equal(cmp$lrt$n_negative, 1)
  expect_output(print(cmp), "full model fits worse")
})

test_that("models fitted to different data, and invalid input, are rejected", {
  fit_other <- quiet_driftd(subset(d_effect, id != 4), sub_id = "id", response = "resp", rt = "rt",
                            compute_se = FALSE)
  expect_error(compare_models(fit_const, fit_other), "not fitted to the same data")
  expect_error(compare_models(fit_const), "at least two fits")
  expect_error(compare_models(fit_const, list()), "driftd_results")
  expect_error(compare_models(fit_const, fit_vary, names = c("a", "a")), "unique")
  expect_equal(compare_models(fit_const, fit_vary, names = c("m1", "m2"))$table$model, c("m1", "m2"))
})

test_that("print and plot methods work", {
  cmp <- compare_models(constant = fit_const, varying = fit_vary)
  expect_output(print(cmp), "Lowest BIC: varying")
  expect_output(print(cmp), "Likelihood-ratio tests")
  for (criterion in c("BIC", "AIC")) {
    g <- plot(cmp, criterion = criterion)
    expect_s3_class(g, "ggplot")
    expect_no_error(ggplot2::ggplot_build(g))
  }
})
