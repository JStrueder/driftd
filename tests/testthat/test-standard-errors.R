d <- subset(example_data, ID %in% c(1, 2))
fit <- quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_v = "val")
estimated <- c("a", "t0", "z", "v_neg", "v_pos")

test_that("standard errors are returned for all estimated parameters", {
  se <- fit$standard_errors
  expect_equal(names(se), c("ID", estimated))
  expect_equal(se$ID, fit$individual_estimates$ID)
  expect_true(all(is.finite(as.matrix(se[, estimated]))))
  expect_true(all(as.matrix(se[, estimated]) > 0))
  expect_equal(fit$fit_info$se_status, c("ok", "ok"))
})

test_that("vcov matrices are symmetric, positive definite and match the SEs", {
  expect_equal(names(fit$vcov), c("1", "2"))
  for (id in names(fit$vcov)) {
    V <- fit$vcov[[id]]
    expect_equal(dimnames(V), list(estimated, estimated))
    expect_equal(V, t(V))
    expect_gt(min(eigen(V, only.values = TRUE)$values), 0)
    se_row <- unlist(fit$standard_errors[fit$standard_errors$ID == as.numeric(id), estimated])
    expect_equal(sqrt(diag(V)), se_row)
  }
})

test_that("standard errors agree with an independent Hessian (stats::optimHess)", {
  # Same model written directly with rtdists, Hessian from a different numerical method.
  d1 <- subset(d, ID == 1)
  resp <- ifelse(d1$R == 1, "upper", "lower")
  is_pos <- d1$val == "pos"
  nll <- function(p) {
    v <- ifelse(is_pos, p[["v_pos"]], p[["v_neg"]])
    -sum(log(rtdists::ddiffusion(d1$RT, resp, a = p[["a"]], t0 = p[["t0"]], z = p[["z"]] * p[["a"]], v = v,
                                 precision = 4)))
  }
  est <- unlist(fit$individual_estimates[fit$individual_estimates$ID == 1, estimated])
  H <- stats::optimHess(est, nll, control = list(ndeps = rep(1e-4, length(est))))
  se_direct <- sqrt(diag(solve(H)))
  se_driftd <- unlist(fit$standard_errors[fit$standard_errors$ID == 1, estimated])
  expect_equal(se_driftd, se_direct, tolerance = 0.02)
})

test_that("compute_se = FALSE skips the standard errors", {
  fit_no <- quiet_driftd(subset(d, ID == 1), sub_id = "ID", response = "R", rt = "RT", compute_se = FALSE)
  expect_true(all(is.na(as.matrix(fit_no$standard_errors[, -1]))))
  expect_null(fit_no$vcov[["1"]])
  expect_equal(fit_no$fit_info$se_status, "not computed")
})

test_that("a parameter at its bound gets no SE, and an sv stuck at 0 triggers a retry", {
  # Data simulated without sv: the sv estimate is 0, at its lower bound.
  d_sv <- simulate_ddm(300, a = 1.5, t0 = 0.3, z = 0.45, v = c(easy = 2, hard = 0.7), seed = 1)
  fit_sv <- quiet_driftd(d_sv, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond",
                         fix_sv = FALSE, max_retries = 1)
  expect_equal(fit_sv$individual_estimates$sv, 0)
  expect_equal(fit_sv$fit_info$optimizer_runs, 2)         # the start with sv moved to 1
  expect_true(is.na(fit_sv$standard_errors$sv))
  expect_match(fit_sv$fit_info$se_status, "no SE at bound: sv")
  expect_true(all(is.finite(unlist(fit_sv$standard_errors[, c("a", "t0", "z", "v_easy", "v_hard")]))))
})

test_that("max_retries = 0 never retries", {
  d_sv <- simulate_ddm(300, a = 1.5, t0 = 0.3, z = 0.45, v = c(easy = 2, hard = 0.7), seed = 1)
  fit_sv <- quiet_driftd(d_sv, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond",
                         fix_sv = FALSE, max_retries = 0, compute_se = FALSE)
  expect_equal(fit_sv$fit_info$optimizer_runs, 1)
})
