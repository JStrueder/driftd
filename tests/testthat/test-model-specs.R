# Every kind of model specification: each parameter varied across conditions, the
# variability parameters estimated, parameters fixed, and designs with several factors.
# Each fit must converge, keep its estimates inside the valid parameter space, and work
# with every function that takes a fit.

d_spec <- edge_data()
set.seed(2)
d_spec$speed <- sample(c("fast", "slow"), nrow(d_spec), replace = TRUE)
d_spec$block <- sample(c("b1", "b2", "b3"), nrow(d_spec), replace = TRUE)

specs <- list(
  vary_all_main = list(vary_a = "cond", vary_t0 = "cond", vary_z = "cond", vary_v = "cond"),
  crossed_a_v = list(vary_a = "speed", vary_v = "cond"),
  three_factors = list(vary_v = c("cond", "speed"), vary_t0 = "block"),
  estimate_sv = list(vary_v = "cond", fix_sv = FALSE),
  estimate_sz = list(vary_v = "cond", fix_sz = FALSE),
  estimate_st0 = list(vary_v = "cond", fix_st0 = FALSE),
  estimate_all_variability = list(vary_v = "cond", fix_sv = FALSE, fix_sz = FALSE, fix_st0 = FALSE),
  vary_sv = list(vary_v = "cond", fix_sv = FALSE, vary_sv = "cond"),
  fix_all_but_v = list(vary_v = "cond", fix_a = TRUE, a_fixed = 1.2, fix_t0 = TRUE, t0_fixed = 0.25,
                       fix_z = TRUE, z_fixed = 0.5),
  fix_v = list(fix_v = TRUE, v_fixed = 1),
  fixed_variability = list(vary_v = "cond", sv_fixed = 0.5, sz_fixed = 0.2, st0_fixed = 0.1)
)

for (spec in names(specs)) {
  test_that(paste("model specification works end to end:", spec), {
    skip_on_cran()
    fit <- do.call(quiet_driftd, c(list(d_spec, sub_id = "id", response = "resp", rt = "rt",
                                        compute_se = FALSE), specs[[spec]]))
    est <- fit$individual_estimates
    expect_equal(est$converge, c(0, 0, 0))
    expect_true(all(is.finite(est$LLE)))

    # estimates are inside the valid parameter space
    info <- fit$settings$parameters
    value <- function(type) as.matrix(est[info$name[info$parameter == type]])
    expect_true(all(value("a") > 0))
    expect_true(all(value("z") > 0 & value("z") < 1))
    expect_true(all(value("t0") >= 0 & value("t0") < min(d_spec$rt)))
    for (type in c("sv", "sz", "st0")) expect_true(all(value(type) >= 0))
    # all starting points lie between the boundaries (sz is on the scale of a)
    cells <- estimates(fit, format = "cells")
    expect_true(all(cells$sz < 2 * pmin(cells$z, 1 - cells$z) * cells$a + 1e-8))

    # the stored specification refits to the same estimates
    expect_equal(estimate_values(driftd:::refit_model(fit, fit$data)), estimate_values(fit),
                 tolerance = 1e-6)
    expect_downstream_works(fit)
  })
}

test_that("standard errors and parameter recovery work with estimated variability parameters", {
  skip_on_cran()
  d <- simulate_ddm(400, a = 1.4, t0 = 0.3, z = 0.5, v = c(easy = 2, hard = 0.8), sv = 1, ids = 1:2, seed = 5)
  fit <- quiet_driftd(d, sub_id = "id", response = "resp", rt = "rt", vary_v = "cond", fix_sv = FALSE)
  expect_true(all(fit$individual_estimates$sv > 0.2))
  expect_true(all(is.finite(unlist(fit$standard_errors[, c("a", "t0", "z", "v_easy", "v_hard", "sv")]))))
  rec <- recover_parameters(fit, seed = 1, verbose = FALSE)
  expect_setequal(rec$summary$name, c("a", "t0", "z", "v_easy", "v_hard", "sv"))
  expect_true(all(rec$data$converge == 0))
})
