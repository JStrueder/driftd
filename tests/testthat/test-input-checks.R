d <- subset(example_data, ID %in% c(1, 2))

test_that("data must be a data frame with the named columns", {
  expect_error(quiet_driftd(as.matrix(d), sub_id = "ID", response = "R", rt = "RT"), "must be a data frame")
  expect_error(quiet_driftd(d, sub_id = "ID", response = "R", rt = "reaction_time"), "not found in data: reaction_time")
})

test_that("missing values and non-positive RTs are rejected", {
  d_na <- d
  d_na$RT[5] <- NA
  expect_error(quiet_driftd(d_na, sub_id = "ID", response = "R", rt = "RT"), "NA values found in RT column")

  d_neg <- d
  d_neg$RT[5] <- -0.1
  expect_error(quiet_driftd(d_neg, sub_id = "ID", response = "R", rt = "RT"), "must be positive")
})

test_that("responses must be coded 0/1 or lower/upper", {
  d_bad <- d
  d_bad$R <- ifelse(d_bad$R == 1, "yes", "no")
  expect_error(quiet_driftd(d_bad, sub_id = "ID", response = "R", rt = "RT"), "must be coded as 0/1")
})

test_that("a parameter cannot be both fixed and varied", {
  expect_error(quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_a = "val", fix_a = TRUE),
               "both fixed and set to vary across conditions: a")
})

test_that("condition columns must exist and have no missing values", {
  expect_error(quiet_driftd(d, sub_id = "ID", response = "R", rt = "RT", vary_v = "nope"),
               "vary_ columns were not found in data: nope")

  d_na <- d
  d_na$val[3] <- NA
  expect_error(quiet_driftd(d_na, sub_id = "ID", response = "R", rt = "RT", vary_v = "val"),
               "NA values found in condition column")
})

test_that("condition labels that collide because of '_' are rejected", {
  d_col <- d[1:300, ]
  d_col$f1 <- rep(c("a_b", "a"), 150)
  d_col$f2 <- rep(c("c", "b_c"), 150)
  expect_error(quiet_driftd(d_col, sub_id = "ID", response = "R", rt = "RT", vary_v = c("f1", "f2")),
               "contain '_'")
})
