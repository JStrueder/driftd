d_rec <- simulate_ddm(150, a = 1.4, t0 = 0.3, z = 0.5, v = c(easy = 2.0, hard = 0.6), ids = 1:4, seed = 8)
fit_rec <- quiet_driftd(d_rec, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond", compute_se = FALSE)

test_that("driftd() stores the model specification so it can be refitted", {
  spec <- fit_rec$settings$model
  expect_equal(spec$vary_v, "cond")
  expect_true(spec$fix_sv)
  expect_equal(spec$sv_fixed, 0)
  refit <- driftd:::refit_model(fit_rec, fit_rec$data)
  expect_equal(refit$individual_estimates$LLE, fit_rec$individual_estimates$LLE, tolerance = 1e-8)
})

test_that("recover_parameters() compares true and recovered values of every estimated parameter", {
  rec <- recover_parameters(fit_rec, seed = 1, verbose = FALSE)
  expect_s3_class(rec, "driftd_recovery")
  expect_equal(names(rec$data), c("rep", "id", "parameter", "cond", "name", "true", "recovered", "converge"))
  expect_equal(nrow(rec$data), 4 * 5)   # 4 participants x (a, t0, z, v_easy, v_hard)
  truth <- estimates(fit_rec, include_fixed = FALSE)
  expect_equal(rec$data$true, truth$estimate)
  # recovered values are close to the true ones
  expect_lt(mean(abs(rec$data$recovered - rec$data$true)), 0.15)
  expect_equal(rec$summary$name, c("a", "t0", "z", "v_easy", "v_hard"))
  expect_equal(names(rec$summary), c("name", "parameter", "n", "r", "bias", "rmse"))
  expect_equal(rec$summary_by_type$parameter, c("a", "t0", "z", "v"))
  expect_equal(rec$summary_by_type$n[rec$summary_by_type$parameter == "v"], 8)
  expect_equal(dim(rec$error_correlation), c(5, 5))
  expect_equal(unname(diag(rec$error_correlation)), rep(1, 5))
})

test_that("recover_parameters() supports replications, seeds and standard errors", {
  rec2 <- recover_parameters(fit_rec, n_rep = 2, seed = 3, compute_se = TRUE, verbose = FALSE)
  expect_equal(nrow(rec2$data), 2 * 4 * 5)
  expect_equal(sort(unique(rec2$data$rep)), 1:2)
  expect_true("coverage" %in% names(rec2$summary))
  expect_true(all(rec2$summary$coverage >= 0 & rec2$summary$coverage <= 1))
  expect_equal(recover_parameters(fit_rec, seed = 3, verbose = FALSE)$data,
               recover_parameters(fit_rec, seed = 3, verbose = FALSE)$data)
})

test_that("fixed parameters stay fixed in the refits", {
  fit_z <- quiet_driftd(d_rec, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond",
                        fix_z = TRUE, z_fixed = 0.45, compute_se = FALSE)
  rec <- recover_parameters(fit_z, seed = 1, verbose = FALSE)
  expect_false("z" %in% rec$data$parameter)
  expect_equal(fit_z$settings$model$z_fixed, 0.45)
})

test_that("recover_models() recovers the generating model", {
  fit_const <- quiet_driftd(d_rec, sub_id = "id", response = "resp", rt = "rt", compute_se = FALSE)
  mr <- recover_models(constant = fit_const, varying = fit_rec, seed = 2, verbose = FALSE)
  expect_s3_class(mr, "driftd_model_recovery")
  for (crit in c("AIC", "BIC")) {
    m <- mr$confusion[[crit]]
    expect_equal(dimnames(m), list(generating = c("constant", "varying"), selected = c("constant", "varying")))
    expect_equal(unname(rowSums(m)), c(1, 1))
    expect_equal(unname(rowSums(mr$confusion_group[[crit]])), c(1, 1))
  }
  # data generated with a strong drift effect are recognised as such for every participant
  expect_equal(mr$confusion$BIC["varying", "varying"], 1)
  expect_equal(mr$confusion_group$BIC["varying", "varying"], 1)
  expect_equal(names(mr$comparisons), c("constant_rep1", "varying_rep1"))
  expect_error(recover_models(fit_rec), "at least two fits")
})

test_that("print and plot methods work for recovery results", {
  rec <- recover_parameters(fit_rec, seed = 1, verbose = FALSE)
  expect_output(print(rec), "Parameter recovery")
  for (by in c("parameter", "name")) expect_no_error(ggplot2::ggplot_build(plot(rec, by = by)))

  fit_const <- quiet_driftd(d_rec, sub_id = "id", response = "resp", rt = "rt", compute_se = FALSE)
  mr <- recover_models(constant = fit_const, varying = fit_rec, seed = 2, verbose = FALSE)
  expect_output(print(mr), "Model recovery")
  expect_no_error(ggplot2::ggplot_build(plot(mr)))
  expect_no_error(ggplot2::ggplot_build(plot(mr, criterion = "AIC", level = "group")))
})

test_that("results from older driftd versions are rejected clearly", {
  old <- fit_rec
  old$settings$model <- NULL
  expect_error(recover_parameters(old), "older version of driftd")
})
