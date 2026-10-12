test_that("parameters are recovered from simulated data", {
  truth <- list(a = 1.5, t0 = 0.3, z = 0.45, v = c(easy = 2.0, hard = 0.7))
  d <- simulate_ddm(1500, a = truth$a, t0 = truth$t0, z = truth$z, v = truth$v, seed = 11)
  fit <- quiet_driftd(d, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond", compute_se = FALSE)
  est <- fit$individual_estimates

  expect_equal(est$converge, 0)
  expect_near(est$a, truth$a, 0.1)
  expect_near(est$t0, truth$t0, 0.02)
  expect_near(est$z, truth$z, 0.03)
  expect_near(est$v_easy, truth$v[["easy"]], 0.25)
  expect_near(est$v_hard, truth$v[["hard"]], 0.25)
})

test_that("estimates are a maximum of an independently written likelihood", {
  # The same model written directly with rtdists: drift varies by val, a, t0 and z are
  # shared, no variability parameters.
  d <- subset(example_data, ID == 1)
  fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_v = "val", compute_se = FALSE)
  est <- fit$individual_estimates
  pars <- c("a", "t0", "z", "v_neg", "v_pos")
  driftd_par <- unlist(est[, pars])

  resp <- ifelse(d$R == 1, "upper", "lower")
  is_pos <- d$val == "pos"
  nll <- function(p) {
    v <- ifelse(is_pos, p[["v_pos"]], p[["v_neg"]])
    -sum(log(rtdists::ddiffusion(d$RT, resp, a = p[["a"]], t0 = p[["t0"]], z = p[["z"]] * p[["a"]], v = v)))
  }
  lower <- c(0.01, 0, 0, -Inf, -Inf)
  upper <- c(Inf, min(d$RT), 1, Inf, Inf)

  # driftd's reported log-likelihood is the likelihood of its estimates
  expect_near(est$LLE, -nll(driftd_par), 1e-6)

  # at least as good as a direct fit from a generic start
  direct <- stats::nlminb(c(a = 1, t0 = 0.1, z = 0.5, v_neg = 0, v_pos = 0), nll, lower = lower, upper = upper)
  expect_gte(est$LLE, -direct$objective - 1e-4)

  # a direct fit started at driftd's estimates stays there: they are a maximum
  polish <- stats::nlminb(driftd_par, nll, lower = lower, upper = upper)
  expect_near(polish$par, driftd_par, 1e-3)
  expect_near(-polish$objective, est$LLE, 1e-4)
})

test_that("RTs in milliseconds give the same estimates as seconds", {
  d <- subset(example_data, ID == 1)
  d_ms <- d
  d_ms$RT <- d_ms$RT * 1000
  fit_s <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_v = "val", compute_se = FALSE)
  expect_message(fit_ms <- driftd_no_output(d_ms, sub_id = "ID", response = "R", rt = "RT",
                                            vary_v = "val", compute_se = FALSE),
                 "milliseconds")
  expect_equal(fit_ms$individual_estimates, fit_s$individual_estimates, tolerance = 1e-5)
})

test_that("overlapping level names in different factors are matched correctly", {
  # x has levels "0.5"/"0.6", which are also levels of prob. Relabelling x to "A"/"B"
  # must give exactly the same fit.
  d <- subset(example_data, ID == 1)
  d$x <- ifelse(d$val == "pos", "0.5", "0.6")
  d_relabelled <- d
  d_relabelled$x <- ifelse(d$val == "pos", "A", "B")

  fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_v = c("prob", "x"), compute_se = FALSE)
  fit_rel <- quiet_driftd(d_relabelled, sub_id = "ID", response = "R", rt = "RT", vary_v = c("prob", "x"), compute_se = FALSE)

  expect_equal(fit$individual_estimates$LLE, fit_rel$individual_estimates$LLE, tolerance = 1e-8)
  expect_equal(fit$individual_estimates$v_0.6_0.5, fit_rel$individual_estimates$v_0.6_A, tolerance = 1e-8)
  # drift rates actually moved away from their starting values
  expect_gt(stats::sd(unlist(fit$individual_estimates[, grep("^v_", names(fit$individual_estimates))])), 0.1)
})

test_that("condition levels containing '_' are handled", {
  d <- subset(example_data, ID == 1)
  d_us <- d
  d_us$val <- paste0(d_us$val, "_cue")
  fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_v = "val", compute_se = FALSE)
  fit_us <- quiet_driftd(d_us, sub_id = "ID", response = "R", rt = "RT", vary_v = "val", compute_se = FALSE)
  expect_equal(fit_us$individual_estimates$LLE, fit$individual_estimates$LLE, tolerance = 1e-8)
  expect_true(all(c("v_neg_cue", "v_pos_cue") %in% names(fit_us$individual_estimates)))
})

test_that("unused levels of a factor ID column are not reported as participants", {
  d <- subset(example_data, ID == 1)
  d$ID <- factor(d$ID, levels = 1:33)
  msgs <- character()
  withCallingHandlers(
    quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", compute_se = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("very few total trials", msgs)))
})

test_that("t0 stays below the fastest RT", {
  d <- subset(example_data, ID %in% c(1, 2))
  fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_v = "val", compute_se = FALSE)
  min_rt <- tapply(d$RT, d$ID, min)
  expect_true(all(fit$individual_estimates$t0 < min_rt[as.character(fit$individual_estimates$ID)]))
})

test_that("the reported likelihood is the true likelihood (no trials written off near t0)", {
  # With the former density floor of 1e-10, the fit for this participant put t0 a fraction
  # of a millisecond below the fastest RT, so that trial's density underflowed to 0 and was
  # floored. The true likelihood of those estimates was 0 (log-likelihood -Inf).
  d <- subset(example_data, ID == 6)
  fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT",
                      vary_z = "val", vary_v = c("prob", "val"), compute_se = FALSE)
  est <- fit$individual_estimates
  v <- unlist(est[paste0("v_", d$prob, "_", d$val)])
  z <- unlist(est[paste0("z_", d$val)])
  dens <- rtdists::ddiffusion(d$RT, ifelse(d$R == 1, "upper", "lower"),
                              a = est$a, v = v, t0 = est$t0, z = z * est$a, precision = 5)
  expect_true(all(dens > 0))
  expect_near(est$LLE, sum(log(dens)), 0.01)
  expect_gt(min(d$RT) - est$t0, 0.001)
})
