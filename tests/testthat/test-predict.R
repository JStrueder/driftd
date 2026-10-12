d <- subset(example_data, ID %in% c(1, 2))
fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_z = "val", vary_v = c("prob", "val"),
                    compute_se = FALSE)
pred <- predict(fit, n_sim = 100, seed = 1)

test_that("predict() returns one row per participant, cell, response and statistic", {
  expect_s3_class(pred, "driftd_prediction")
  expect_equal(nrow(pred), 2 * 10 * 2 * 6)
  expect_equal(names(pred), c("ID", ".cell", "val", "prob", "response", "statistic", "quantile", "n_cell",
                              "n_response", "observed", "predicted", "lower", "upper"))
  # the two response probabilities of a cell add up to 1
  p <- pred[pred$statistic == "probability", ]
  expect_equal(as.vector(tapply(p$predicted, list(p$ID, p$.cell), sum)), rep(1, 20), tolerance = 1e-6)
  expect_equal(as.vector(tapply(p$observed, list(p$ID, p$.cell), sum)), rep(1, 20))
})

test_that("predictions are the exact rtdists probabilities and quantiles", {
  e <- fit$individual_estimates[fit$individual_estimates$ID == 1, ]
  pars <- list(a = e$a, v = e$v_0.4_neg, t0 = e$t0, z = e$z_neg * e$a)
  rows <- pred[pred$ID == 1 & pred$.cell == "neg_0.4" & pred$response == "lower", ]
  p_lower <- rtdists::pdiffusion(Inf, "lower", a = pars$a, v = pars$v, t0 = pars$t0, z = pars$z)
  q_lower <- rtdists::qdiffusion(c(0.1, 0.3, 0.5, 0.7, 0.9) * p_lower, "lower",
                                 a = pars$a, v = pars$v, t0 = pars$t0, z = pars$z)
  expect_equal(rows$predicted, c(p_lower, q_lower), tolerance = 1e-6)
})

test_that("observed statistics are computed per cell, with quantiles only from min_n responses", {
  cell <- d[d$ID == 1 & d$prob == 0.4 & d$val == "neg", ]
  rows <- pred[pred$ID == 1 & pred$.cell == "neg_0.4", ]
  lower_rt <- cell$RT[cell$R == 0]
  expect_equal(rows$n_cell[1], nrow(cell))
  expect_equal(rows$observed[rows$response == "lower" & rows$statistic == "probability"], mean(cell$R == 0))
  expect_equal(rows$observed[rows$response == "lower" & rows$statistic == "quantile"],
               unname(stats::quantile(lower_rt, c(0.1, 0.3, 0.5, 0.7, 0.9))))
  # cells with fewer than min_n responses have no observed quantiles
  few <- pred$statistic == "quantile" & pred$n_response < 5
  expect_true(any(few))
  expect_true(all(is.na(pred$observed[few])))
  expect_true(all(!is.na(pred$observed[pred$statistic == "quantile" & pred$n_response >= 5])))
})

test_that("prediction intervals are ordered, reproducible with a seed, and optional", {
  ok <- !is.na(pred$lower)
  expect_true(all(pred$lower[ok] <= pred$upper[ok]))
  expect_identical(predict(fit, n_sim = 100, seed = 1), pred)
  no_int <- predict(fit, intervals = FALSE)
  expect_true(all(is.na(no_int$lower)))
  expect_equal(no_int$predicted, pred$predicted)
  expect_error(summary(no_int), "needs prediction intervals")
})

test_that("predict(), simulate() and their seeds leave the user's random numbers unchanged", {
  set.seed(42)
  a <- stats::runif(1)
  set.seed(42)
  invisible(predict(fit, n_sim = 20, seed = 7))
  invisible(simulate(fit, seed = 7))
  expect_equal(stats::runif(1), a)
})

test_that("observed statistics fall outside prediction intervals at about the nominal rate", {
  # data simulated from the model itself: about 5% of statistics should fall outside
  sim_d <- simulate_ddm(250, a = 1.4, t0 = 0.3, z = 0.5, v = c(easy = 1.8, hard = 0.6), ids = 1:4, seed = 3)
  sim_fit <- quiet_driftd(sim_d, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond", compute_se = FALSE)
  s <- summary(predict(sim_fit, n_sim = 200, seed = 2))
  expect_lt(sum(s$n_outside) / sum(s$n_statistics), 0.12)
})

test_that("aggregate() averages participant statistics (Vincentizing)", {
  g <- aggregate(pred)
  expect_equal(names(g), c("val", "prob", "response", "statistic", "quantile", "n_participants",
                           "observed", "predicted", "lower", "upper"))
  row <- g[g$val == "pos" & g$prob == 0.6 & g$response == "upper" & g$statistic == "quantile" & g$quantile == 0.5, ]
  src <- pred[pred$val == "pos" & pred$prob == 0.6 & pred$response == "upper" & pred$statistic == "quantile" &
                pred$quantile == 0.5 & !is.na(pred$observed), ]
  expect_equal(row$observed, mean(src$observed))
  expect_equal(row$predicted, mean(src$predicted))
  expect_equal(row$n_participants, nrow(src))
  # averaging over a condition column
  g_val <- aggregate(pred, by = "val")
  expect_equal(nrow(g_val), 2 * 2 * 6)
  expect_error(aggregate(pred, by = "nope"), "must name condition columns")
})

test_that("summary() lists participants from worst to best fit", {
  s <- summary(pred)
  expect_equal(names(s), c("ID", "n_statistics", "n_outside", "prop_outside", "mean_abs_dev"))
  expect_equal(nrow(s), 2)
  expect_true(all(diff(s$prop_outside) <= 0))
})

test_that("simulate() reproduces the design of the observed data", {
  sims <- simulate(fit, nsim = 2, seed = 1)
  expect_equal(names(sims), c("sim", "ID", ".cell", "val", "prob", "R", "RT"))
  expect_equal(nrow(sims), 2 * nrow(d))
  expect_equal(levels(sims$R), c("lower", "upper"))
  obs <- table(fit$data$ID, fit$data$.cell)
  expect_equal(unclass(table(sims$ID[sims$sim == 2], sims$.cell[sims$sim == 2])), unclass(obs))
  expect_identical(simulate(fit, nsim = 2, seed = 1), sims)
})

test_that("parameter uncertainty can be included in the simulations", {
  fit_se <- quiet_driftd(subset(d, ID == 1), sub_id = "ID", response = "R", rt = "RT", vary_v = "val")
  pu <- predict(fit_se, n_sim = 20, parameter_uncertainty = TRUE, seed = 1)
  expect_true(all(is.finite(pu$lower[pu$statistic == "probability"])))
  expect_equal(nrow(simulate(fit_se, nsim = 3, parameter_uncertainty = TRUE, seed = 1)), 3 * nrow(subset(d, ID == 1)))
})

test_that("results from older driftd versions are rejected clearly", {
  old <- fit
  old$settings$columns <- NULL
  expect_error(predict(old), "older version of driftd")
  expect_error(simulate(old), "older version of driftd")
})

test_that("all plot types build for designs with two, one and no condition columns", {
  fit1 <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_v = "val", compute_se = FALSE)
  fit0 <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", compute_se = FALSE)
  for (p in list(pred, predict(fit1, n_sim = 20, seed = 1), predict(fit0, n_sim = 20, seed = 1))) {
    for (type in c("qp", "scatter", "cdf", "participants")) {
      g <- plot(p, type = type)
      expect_s3_class(g, "ggplot")
      expect_no_error(suppressMessages(ggplot2::ggplot_build(g)))
    }
  }
})

test_that("plot options select, facet, color and filter as requested", {
  build <- function(g) suppressMessages(ggplot2::ggplot_build(g))
  expect_no_error(build(plot(pred, facet = "val")))
  expect_no_error(build(plot(pred, type = "scatter", level = "participant", color = "val")))
  expect_no_error(build(plot(pred, type = "scatter", intervals = FALSE)))
  expect_no_error(build(plot(pred, participant = 2)))
  expect_no_error(build(plot(pred, type = "cdf", by = "prob", select = list(prob = c(0.4, 0.6)))))
  expect_message(plot(pred, type = "cdf", max_cells = 4), "Showing the first 4 of 10 cells")
  expect_error(plot(pred, facet = "nope"), "condition columns")
  expect_error(plot(pred, participant = 99), "No participant")
  expect_error(plot(pred, select = list(prob = 0.99)), "No data left")
})
