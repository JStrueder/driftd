test_that("EZ-diffusion reproduces the worked example of Wagenmakers et al. (2007)", {
  # The paper reports v = 0.0999, a = 0.1400, Ter = 0.300 for Pc = .802, VRT = .112,
  # MRT = .723, with s = 0.1. driftd uses s = 1, which scales v and a by 10.
  ez <- driftd:::ez_diffusion(prop_upper = 0.802, mean_rt = 0.723, var_rt = 0.112, n_trials = 1000)
  expect_equal(ez$v, 0.999, tolerance = 1e-3)
  expect_equal(ez$a, 1.400, tolerance = 1e-3)
  expect_equal(ez$t0, 0.300, tolerance = 1e-3)
})

test_that("EZ-diffusion handles edge cases", {
  # below .5 upper responses gives a negative drift with the same boundary separation
  ez <- driftd:::ez_diffusion(c(0.802, 0.198), c(0.723, 0.723), c(0.112, 0.112), c(1000, 1000))
  expect_equal(ez$v[2], -ez$v[1], tolerance = 1e-6)
  expect_equal(ez$a[2], ez$a[1], tolerance = 1e-6)

  # all responses at one boundary are edge-corrected instead of giving infinite values
  expect_true(all(is.finite(unlist(driftd:::ez_diffusion(1, 0.6, 0.05, 100)))))

  # undefined cells (one trial, zero variance) give NA
  expect_true(all(is.na(unlist(driftd:::ez_diffusion(0.8, 0.6, NA, 1)))))
  expect_true(all(is.na(unlist(driftd:::ez_diffusion(0.8, 0.6, 0, 50)))))
})

test_that("halton() gives the van der Corput sequence", {
  expect_equal(sapply(1:6, driftd:::halton, base = 2), c(0.5, 0.25, 0.75, 0.125, 0.625, 0.375))
  expect_equal(sapply(1:3, driftd:::halton, base = 3), c(1 / 3, 2 / 3, 1 / 9))
})
