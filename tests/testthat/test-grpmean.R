if (require("testthat") && require("sjstats") && require("dplyr")) {
  context("sjstats, grpmean")

  data(efc)
  efc$weight <- abs(rnorm(n = nrow(efc), mean = 1, sd = .5))
  efc_grouped <- group_by(efc, c172code)

  test_that("grpmean", {
    grpmean(efc, c12hour, e42dep)
  })

  test_that("grpmean, weighting", {
    w <- "weight"
    grpmean(efc, c12hour, e42dep, weights = weight)
    grpmean(efc, c12hour, e42dep, weights = "weight")
    grpmean(efc, c12hour, e42dep, weights = w)
  })

  test_that("grpmean, grouping", {
    grpmean(efc_grouped, c12hour, e42dep)
  })

  test_that("grpmean, grouped weighting", {
    w <- "weight"
    grpmean(efc_grouped, c12hour, e42dep, weights = weight)
    grpmean(efc_grouped, c12hour, e42dep, weights = "weight")
    grpmean(efc_grouped, c12hour, e42dep, weights = w)
  })

}
