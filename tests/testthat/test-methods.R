d <- subset(example_data, ID %in% c(1, 2))
fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_z = "val", vary_v = c("prob", "val"),
                    compute_se = FALSE)

test_that("the results object has the documented structure", {
  expect_s3_class(fit, "driftd_results")
  expect_true(all(c("individual_estimates", "summary", "settings", "data", "convergence",
                    "fit_info", "standard_errors", "vcov") %in% names(fit)))
  # the processed data keep the original condition columns
  expect_equal(names(fit$data), c("ID", ".cell", "R", "RT", "val", "prob"))
  expect_equal(fit$settings$condition_columns, c("val", "prob"))
  expect_equal(nrow(fit$settings$cells), 10)
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

test_that("parameter plots run without error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, type = "convergence"))
  expect_error(plot(fit, type = "nope"), "Unknown plot type")
})
