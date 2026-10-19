# The outlier mixture (p_outlier), as in HDDM: each trial's likelihood is
# (1 - p_outlier) * ddiffusion + p_outlier * w_outlier.

# One participant, two drift conditions, with 2% fast guesses (0.1-0.2 s, random response)
d_clean <- simulate_ddm(250, a = 1.2, t0 = 0.3, z = 0.5, v = c(easy = 2, hard = 1), seed = 3)
d_guess <- d_clean
set.seed(4)
guesses <- sample(nrow(d_guess), round(0.02 * nrow(d_guess)))
d_guess$rt[guesses] <- runif(length(guesses), 0.1, 0.2)
d_guess$resp[guesses] <- sample(0:1, length(guesses), replace = TRUE)

fit_out <- function(d, ...) {
  quiet_driftd(d, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond", ...)
}

test_that("p_outlier and w_outlier are validated", {
  expect_error(fit_out(d_clean, p_outlier = -0.1), "p_outlier must be")
  expect_error(fit_out(d_clean, p_outlier = 1), "p_outlier must be")
  expect_error(fit_out(d_clean, p_outlier = c(0.05, 0.1)), "p_outlier must be")
  expect_error(fit_out(d_clean, p_outlier = 0.05, w_outlier = 0), "w_outlier must be")
})

test_that("the log-likelihood is that of the HDDM mixture", {
  fit <- fit_out(d_guess, p_outlier = 0.05, compute_se = FALSE)
  est <- fit$individual_estimates
  dens <- unlist(lapply(c("easy", "hard"), function(cond) {
    d <- d_guess[d_guess$cond == cond, ]
    rtdists::ddiffusion(d$rt, response = ifelse(d$resp == 1, "upper", "lower"),
                        a = est$a, v = est[[paste0("v_", cond)]], t0 = est$t0, z = est$z * est$a)
  }))
  expect_equal(est$LLE, sum(log(0.95 * dens + 0.05 * 0.1)), tolerance = 1e-8)
  expect_equal(fit$settings$p_outlier, 0.05)
  expect_equal(fit$settings$w_outlier, 0.1)
  expect_equal(fit$settings$model$p_outlier, 0.05)
})

test_that("p_outlier = 0 is the default and leaves the fit unchanged", {
  expect_equal(estimate_values(fit_out(d_clean, compute_se = FALSE)),
               estimate_values(fit_out(d_clean, p_outlier = 0, compute_se = FALSE)))
})

test_that("fast guesses distort the fit without the outlier component, but not with it", {
  pars <- c("a", "t0", "z", "v_easy", "v_hard")
  clean <- fit_out(d_clean, compute_se = FALSE)$individual_estimates[pars]
  plain <- fit_out(d_guess, compute_se = FALSE)$individual_estimates
  robust_fit <- fit_out(d_guess, p_outlier = 0.05)
  robust <- robust_fit$individual_estimates
  # without it, t0 must lie below the fastest guess, and a and z compensate
  expect_lt(plain$t0, 0.2)
  expect_gt(plain$a - clean$a, 0.3)
  # with it, t0 can exceed the fastest RT, and the estimates are close to those from the
  # same data without the guesses
  expect_gt(robust$t0, min(d_guess$rt))
  expect_near(robust[pars], clean, 0.05)
  expect_true(all(is.finite(unlist(robust_fit$standard_errors[-1]))))
})

test_that("on clean data, a small outlier component changes the estimates little", {
  plain <- fit_out(d_clean, compute_se = FALSE)$individual_estimates
  robust <- fit_out(d_clean, p_outlier = 0.05, compute_se = FALSE)$individual_estimates
  expect_near(robust[c("a", "t0", "z", "v_easy", "v_hard")],
              plain[c("a", "t0", "z", "v_easy", "v_hard")], 0.15)
})

test_that("fits with an outlier component work with every downstream function", {
  fit <- fit_out(d_guess, p_outlier = 0.05)
  expect_downstream_works(fit)
  expect_s3_class(predict(fit, n_sim = 20, seed = 1, parameter_uncertainty = TRUE), "driftd_prediction")
})

test_that("comparing models with different p_outlier warns", {
  f0 <- fit_out(d_clean, compute_se = FALSE)
  f1 <- fit_out(d_clean, p_outlier = 0.05, compute_se = FALSE)
  expect_warning(compare_models(f0, f1), "different p_outlier")
})

test_that("parameter recovery refits with the same outlier component", {
  fit <- fit_out(d_clean, p_outlier = 0.05, compute_se = FALSE)
  rec <- suppressMessages(recover_parameters(fit, seed = 1, verbose = FALSE))
  expect_s3_class(rec, "driftd_recovery")
})
