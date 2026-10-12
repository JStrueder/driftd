d <- subset(example_data, ID %in% c(1, 2))
fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_z = "val", vary_v = c("prob", "val"),
                    compute_se = FALSE)

test_that("the results object has the documented structure", {
  expect_s3_class(fit, "driftd_results")
  expect_true(all(c("individual_estimates", "summary", "settings", "data", "convergence",
                    "fit_info", "standard_errors", "vcov") %in% names(fit)))
  expect_equal(nrow(fit$individual_estimates), 2)
  expect_equal(fit$settings$n_participants, 2)
  expect_equal(fit$settings$n_conditions, 10)
  # each condition cell maps to exactly one parameter of each estimated type
  expect_equal(nrow(fit$settings$cell_parameters), 10)
  expect_false(anyNA(fit$settings$cell_parameters[, c("a", "t0", "z", "v")]))
})

test_that("print() shows the estimates and returns the object invisibly", {
  out <- utils::capture.output(vis <- withVisible(print(fit)))
  expect_false(vis$visible)
  expect_identical(vis$value, fit)
  expect_true(any(grepl("v_0.6_pos", out)))
  expect_false(any(out == "NULL"))
})

test_that("summary() and individual_estimates() return the stored tables", {
  utils::capture.output(s <- summary(fit))
  expect_identical(s, fit$summary)
  expect_identical(individual_estimates(fit), fit$individual_estimates)
  expect_error(individual_estimates(list()), "driftd_results")
})

test_that("posterior_predict() simulates one dataset per observed trial", {
  set.seed(1)
  utils::capture.output(ppc <- posterior_predict(fit, n_sim = 2))
  expect_s3_class(ppc, "driftd_ppc")
  expect_equal(nrow(ppc$simulated_data), 2 * nrow(d))
  expect_equal(sort(unique(ppc$simulated_data$n_sim)), 1:2)
  # same number of simulated trials as observed trials in every participant x condition cell
  obs <- table(ppc$observed_data$ID, ppc$observed_data$cond)
  sim <- table(ppc$simulated_data$ID, ppc$simulated_data$cond) / 2
  expect_equal(unclass(sim), unclass(obs))
  expect_error(posterior_predict(list()), "driftd_results")
})

test_that("plot methods run without error", {
  set.seed(1)
  utils::capture.output(ppc <- posterior_predict(fit))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, type = "convergence"))
  expect_no_error(plot(ppc, type = "summary"))
  expect_no_error(utils::capture.output(plot(ppc, type = "participant")))
  expect_no_error(plot(ppc, type = "density", condition = "pos_0.55"))
  expect_error(plot(fit, type = "nope"), "Unknown plot type")
})
