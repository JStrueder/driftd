# Data shapes, codings and model settings beyond the example data. Each case should either
# give a sensible fit or fail with a clear message.

d0 <- edge_data()
fit_edge <- function(d, ...) {
  quiet_driftd(d, sub_id = "id", response = "resp", rt = "rt", compute_se = FALSE, ...)
}
ref <- fit_edge(d0, vary_v = "cond")

test_that("fixed values that put every trial at likelihood 0 are rejected", {
  expect_error(fit_edge(d0, fix_z = TRUE, z_fixed = 0), "strictly between 0 and 1")
  expect_error(fit_edge(d0, fix_z = TRUE, z_fixed = 1), "strictly between 0 and 1")
  # sz is on the scale of a: with a = 0.8 and z = 0.5, starting points range over 0.8 * sz
  expect_error(fit_edge(d0, fix_a = TRUE, a_fixed = 0.8, fix_z = TRUE, sz_fixed = 0.9),
               "sz_fixed \\(0.9\\) must be smaller than")
  expect_s3_class(fit_edge(d0, fix_a = TRUE, a_fixed = 1.2, fix_z = TRUE, sz_fixed = 0.3),
                  "driftd_results")
})

test_that("every response coding gives the same estimates", {
  codings <- list(
    numeric = as.numeric(d0$resp),
    logical = d0$resp == 1,
    factor = factor(d0$resp),
    text = ifelse(d0$resp == 1, "upper", "lower"),
    text_factor = factor(ifelse(d0$resp == 1, "upper", "lower"))
  )
  for (coding in names(codings)) {
    d <- d0
    d$resp <- codings[[coding]]
    expect_equal(estimate_values(fit_edge(d, vary_v = "cond")), estimate_values(ref), label = coding)
  }
})

test_that("participant IDs of any type keep their estimates", {
  ids <- list(
    text = c("s3", "s1", "s2")[d0$id],
    gaps = c(107, 9, 2000)[d0$id],
    factor = factor(d0$id, levels = 3:1)
  )
  for (type in names(ids)) {
    d <- d0
    d$id <- ids[[type]]
    fit <- fit_edge(d, vary_v = "cond")
    # each participant's estimates are those of the same participant in the reference fit
    original <- match(as.character(fit$individual_estimates$id), as.character(unique(ids[[type]])))
    expect_equal(estimate_values(fit), estimate_values(ref)[original, ], ignore_attr = TRUE, label = type)
  }
})

test_that("tibbles, grouped tibbles and condition columns of any type give the same estimates", {
  expect_equal(estimate_values(fit_edge(tibble::as_tibble(d0), vary_v = "cond")), estimate_values(ref))
  grouped <- dplyr::group_by(tibble::as_tibble(d0), id)
  expect_equal(estimate_values(fit_edge(grouped, vary_v = "cond")), estimate_values(ref))

  d <- d0
  d$cond <- factor(d$cond, levels = c("hard", "easy"))
  est <- estimate_values(fit_edge(d, vary_v = "cond"))
  # reordered levels reorder the parameters, so the optimizer takes a slightly different path
  expect_equal(est[names(estimate_values(ref))], estimate_values(ref), tolerance = 1e-5)

  d$cond <- ifelse(d0$cond == "easy", 1, 2)
  est <- estimate_values(fit_edge(d, vary_v = "cond"))
  expect_equal(est$v_1, ref$individual_estimates$v_easy, tolerance = 1e-5)
  expect_equal(est$v_2, ref$individual_estimates$v_hard, tolerance = 1e-5)
})

test_that("a parameter of a condition a participant has no trials in is NA for that participant", {
  d <- d0[!(d0$id == 2 & d0$cond == "hard"), ]
  fit <- quiet_driftd(d, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond")
  est <- fit$individual_estimates
  expect_true(is.na(est$v_hard[est$id == 2]))
  expect_true(all(!is.na(est$v_hard[est$id != 2])))
  # the other parameters of that participant are estimated, with standard errors
  expect_equal(est[est$id == 2, c("a", "t0", "z", "v_easy")],
               estimate_values(fit_edge(d[d$id == 2, ]))[, c("a", "t0", "z", "v")],
               ignore_attr = TRUE, tolerance = 1e-4)
  se <- fit$standard_errors
  expect_true(all(is.finite(unlist(se[se$id == 2, c("a", "t0", "z", "v_easy")]))))
  # AIC and BIC count the parameters actually estimated
  expect_equal(est$AIC[est$id == 2], 2 * 4 - 2 * est$LLE[est$id == 2])
  cmp <- compare_models(constant = fit_edge(d), varying = fit)
  expect_equal(cmp$participants$k[cmp$participants$model == "varying"], c(5, 4, 5))
  expect_equal(cmp$lrt$df, 2)
  # predictions and simulations only cover the cells the participant has data for
  pred <- predict(fit, n_sim = 20, seed = 1)
  expect_false(any(pred$id == 2 & pred$cond == "hard"))
  sim <- simulate(fit, seed = 1)
  expect_equal(nrow(sim[sim$id == 2 & sim$cond == "hard", ]), 0)
  expect_downstream_works(fit)
})

test_that("unbalanced designs, a single participant and few trials are fitted", {
  unbalanced <- fit_edge(d0[-(1:120), ], vary_v = "cond")
  expect_equal(unbalanced$individual_estimates$converge, c(0, 0, 0))

  single <- fit_edge(d0[d0$id == 1, ], vary_v = "cond")
  expect_equal(estimate_values(single), estimate_values(ref)[1, ], ignore_attr = TRUE)
  expect_downstream_works(single)

  few <- do.call(rbind, lapply(split(d0, d0$id), utils::head, 15))
  msgs <- driftd_messages(few, sub_id = "id", response = "resp", rt = "rt", compute_se = FALSE)
  expect_true(any(grepl("very few total trials \\(<20\\): 1 \\(n=15\\)", msgs)))
  fit_few <- fit_edge(few)
  expect_true(all(is.finite(unlist(estimate_values(fit_few)))))
})

test_that("a condition without errors is fitted, and only one-sided participants are warned about", {
  d <- d0
  d$resp[d$id == 1 & d$cond == "easy"] <- 1L
  fit <- quiet_driftd(d, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond")
  expect_equal(fit$individual_estimates$converge, c(0, 0, 0))
  expect_true(all(is.finite(unlist(fit$standard_errors[1, -1]))))
  msgs <- driftd_messages(d, sub_id = "id", response = "resp", rt = "rt", compute_se = FALSE)
  expect_false(any(grepl("one boundary", msgs)))

  d$resp[d$id == 1] <- 1L
  msgs <- driftd_messages(d, sub_id = "id", response = "resp", rt = "rt", compute_se = FALSE)
  expect_true(any(grepl("all their responses at one boundary: 1\\.", msgs)))
})

test_that("data in which every response is at one boundary keep that boundary", {
  # rtdists reads a factor response by its level number, so a factor with the single
  # level "upper" was read as "lower" and gave a mirrored fit (a negative drift rate)
  d <- d0[d0$id == 1, ]
  for (r in list(1L, "upper", 0L, "lower")) {
    d$resp <- r
    fit <- fit_edge(d)
    ll <- sum(log(rtdists::ddiffusion(d$rt, response = rep(if (r %in% list(1L, "upper")) "upper" else "lower", nrow(d)),
      a = fit$individual_estimates$a, v = fit$individual_estimates$v, t0 = fit$individual_estimates$t0,
      z = fit$individual_estimates$z * fit$individual_estimates$a)))
    expect_equal(fit$individual_estimates$LLE, ll, tolerance = 1e-6)
    expect_equal(sign(fit$individual_estimates$v), if (r %in% list(1L, "upper")) 1 else -1)
  }
})

test_that("very fast RTs are warned about, because they pull t0 down", {
  d <- d0
  d$rt[1] <- 0.05
  msgs <- driftd_messages(d, sub_id = "id", response = "resp", rt = "rt", compute_se = FALSE)
  expect_true(any(grepl("RTs below 0.15 s: 1 \\(n=1\\)", msgs)))
  fit <- fit_edge(d, vary_v = "cond")
  expect_lt(fit$individual_estimates$t0[1], 0.05)
  expect_equal(fit$individual_estimates$t0[2:3], ref$individual_estimates$t0[2:3])

  expect_false(any(grepl("RTs below", driftd_messages(d0, sub_id = "id", response = "resp", rt = "rt",
                                                      compute_se = FALSE))))
})

test_that("tied (rounded) RTs and a very slow outlier are fitted", {
  d <- d0
  d$rt <- round(d$rt, 2)
  expect_near(estimate_values(fit_edge(d, vary_v = "cond"))[, c("a", "t0", "z", "v_easy", "v_hard")],
              estimate_values(ref)[, c("a", "t0", "z", "v_easy", "v_hard")], 0.05)

  d <- d0
  d$rt[1] <- 30
  fit <- fit_edge(d, vary_v = "cond")
  expect_equal(fit$individual_estimates$converge, c(0, 0, 0))
  expect_true(all(is.finite(unlist(estimate_values(fit)))))
})

test_that("negative, near-zero and large drift rates and wide boundaries are recovered", {
  cases <- list(
    negative = list(a = 1.2, t0 = 0.3, v = -1.5),
    chance = list(a = 1.0, t0 = 0.3, v = 0),
    strong = list(a = 1.0, t0 = 0.3, v = 4),
    wide = list(a = 3.0, t0 = 0.4, v = 1)
  )
  for (case in names(cases)) {
    p <- cases[[case]]
    d <- simulate_ddm(1000, a = p$a, t0 = p$t0, z = 0.5, v = c(c = p$v), seed = 4)
    fit <- quiet_driftd(d, sub_id = "id", response = "resp", rt = "rt")
    expect_equal(fit$individual_estimates$converge, 0, label = case)
    # each estimate is within 3.5 standard errors of the true value
    truth <- c(a = p$a, t0 = p$t0, z = 0.5, v = p$v)
    z_scores <- (unlist(fit$individual_estimates[names(truth)]) - truth) /
      unlist(fit$standard_errors[names(truth)])
    expect_true(all(abs(z_scores) < 3.5), label = paste(case, paste(round(z_scores, 1), collapse = " ")))
  }
})
