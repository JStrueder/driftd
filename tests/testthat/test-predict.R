d <- subset(example_data, ID %in% c(1, 2))
fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_z = "val", vary_v = c("prob", "val"),
                    compute_se = FALSE)
pred <- predict(fit, n_sim = 100, seed = 1)

test_that("predict() returns one row per participant, cell, response and statistic", {
  expect_s3_class(pred, "driftd_prediction")
  expect_equal(nrow(pred), 2 * 10 * 2 * 6)
  expect_equal(names(pred), c("ID", ".cell", "val", "prob", "response", "statistic", "quantile", "n_cell",
                              "n_response", "observed", "predicted", "lower", "upper", ".row"))
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

test_that("aggregate() and plot() give the same results on subsets made with base R or dplyr", {
  full <- aggregate(pred)
  pos <- full[full$val == "pos", ]
  rownames(pos) <- NULL
  sort_rows <- function(g) {
    g <- g[order(g$prob, g$response, g$statistic, g$quantile), ]
    rownames(g) <- NULL
    g
  }
  # base R subset, a reordered subset, and dplyr::filter all match the full table
  expect_equal(sort_rows(aggregate(pred[pred$val == "pos", ])), sort_rows(pos))
  shuffled <- pred[rev(which(pred$val == "pos")), ]
  expect_equal(sort_rows(aggregate(shuffled)), sort_rows(pos))
  expect_equal(sort_rows(aggregate(dplyr::filter(pred, val == "pos"))), sort_rows(pos))
  expect_no_error(suppressMessages(ggplot2::ggplot_build(plot(pred[pred$ID == 1, ]))))
  # without the .row column, intervals of averages cannot be computed
  no_row <- pred[, names(pred) != ".row"]
  expect_warning(g <- aggregate(no_row), ".row column is missing")
  expect_true(all(is.na(g$lower)))
  expect_error(aggregate.driftd_prediction(as.data.frame(unclass(pred))), "output of predict")
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
    for (type in c("rt", "choice", "scatter", "participants", "qp")) {
      g <- plot(p, type = type)
      expect_s3_class(g, "ggplot")
      expect_no_error(suppressMessages(ggplot2::ggplot_build(g)))
    }
  }
})

test_that("rt and choice plots assign condition columns to roles by default", {
  roles <- driftd:::assign_roles(aggregate(pred), c("val", "prob"), need_x = TRUE)
  expect_equal(roles, list(x = "prob", color = "val", panels = character(0)))
  roles <- driftd:::assign_roles(aggregate(pred), c("val", "prob"))
  expect_equal(roles, list(x = NULL, color = "val", panels = "prob"))
  # user choices are respected, and every column in `by` gets a role
  roles <- driftd:::assign_roles(aggregate(pred), c("val", "prob"), color = "prob")
  expect_equal(roles, list(x = NULL, color = "prob", panels = "val"))
  expect_error(driftd:::assign_roles(aggregate(pred), c("val", "prob"), color = "val", panels = "val"),
               "different condition columns")

  # the rt plot: one column of panels per response, one row per prob level, colors for val
  b <- ggplot2::ggplot_build(plot(pred))
  expect_equal(nrow(b$layout$layout), 5 * 2)
  expect_equal(length(unique(b$data[[length(b$data)]]$colour)), 1)            # predictions in black
  expect_equal(length(unique(b$data[[length(b$data) - 1]]$colour)), 2)        # observed: 2 colors
})

test_that("individual participants are shown, averaged over columns not in `by`", {
  iv <- driftd:::individual_values(pred, by = "val", info = attr(pred, "driftd"))
  # one value per participant, val level, response and statistic
  expect_equal(nrow(iv), nrow(unique(iv[, c("ID", "val", "response", "statistic", "quantile")])))
  src <- pred[pred$ID == 1 & pred$val == "pos" & pred$response == "upper" & pred$statistic == "probability", ]
  expect_equal(iv$observed[iv$ID == 1 & iv$val == "pos" & iv$response == "upper" & iv$statistic == "probability"],
               mean(src$observed))
  n_layers <- function(g) length(g$layers)
  expect_equal(n_layers(plot(pred, individuals = FALSE)), n_layers(plot(pred)) - 1)
  expect_equal(n_layers(plot(pred, participant = 1)), n_layers(plot(pred, individuals = FALSE)))
})

test_that("plot options select, color, panels, labels and units work as requested", {
  build <- function(g) suppressMessages(ggplot2::ggplot_build(g))
  expect_no_error(build(plot(pred, type = "rt", by = "val")))
  expect_no_error(build(plot(pred, type = "choice", x_var = "val", color = "prob")))
  expect_no_error(build(plot(pred, type = "choice", response = "lower")))
  expect_no_error(build(plot(pred, type = "qp", panels = "val")))
  expect_no_error(build(plot(pred, type = "scatter", level = "participant", color = "val")))
  expect_no_error(build(plot(pred, type = "scatter", intervals = FALSE)))
  expect_no_error(build(plot(pred, participant = 2)))
  expect_no_error(build(plot(pred, select = list(prob = c(0.4, 0.6)))))

  g <- plot(pred, rt_unit = "s", response_labels = c(upper = "pos choice", lower = "neg choice"))
  expect_equal(g$labels$x, "RT (s)")
  expect_equal(levels(g$data$response_label), c("pos choice", "neg choice"))
  expect_equal(plot(pred)$labels$x, "RT (ms)")
  expect_equal(plot(pred, type = "choice", response_labels = c(upper = "yes"))$labels$y, "P(yes)")

  expect_error(plot(pred, color = "nope"), "condition columns")
  expect_error(plot(pred, by = "val", panels = "prob"), "condition columns")
  expect_error(plot(pred, participant = 99), "No participant")
  expect_error(plot(pred, select = list(prob = 0.99)), "No data left")
})
