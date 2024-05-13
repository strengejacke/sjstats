test_that("wtd", {
  data(efc)
  set.seed(123)
  efc$weight <- abs(rnorm(nrow(efc), 1, 0.3))
  expect_equal(weighted_se(efc$c12hour, weights = efc$weight), 1.704182, tolerance = 1e-5)
  expect_equal(weighted_se(efc$c12hour, weights = NULL), 1.691623, tolerance = 1e-5)
})
