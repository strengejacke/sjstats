if (require("testthat") && require("sjstats")) {
  data(efc)
  set.seed(123)
  efc$weight <- abs(rnorm(nrow(efc), 1, .3))

  test_that("wtd", {
    expect_equal(weighted_se(efc$c12hour, weights = efc$weight), 1.704182, tolerance = 1e-5)
    expect_equal(weighted_se(efc$c12hour, weights = NULL), 1.691623, tolerance = 1e-5)
  })

  test_that("weighted_ttest", {
    weighted_ttest(efc, e17age, weights = weight)
    weighted_ttest(efc, e17age, c160age, weights = weight)
    weighted_ttest(e17age ~ e16sex + weight, efc)
    weighted_ttest(efc, e17age, c160age, weights = weight, ci.lvl = .8)
  })
}
