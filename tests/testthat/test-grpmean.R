if (require("testthat") && require("sjstats") && require("dplyr")) {
  data(efc)
  efc$weight <- abs(rnorm(n = nrow(efc), mean = 1, sd = .5))
  efc_grouped <- group_by(efc, c172code)

  test_that("means_by_group", {
    means_by_group(efc, c12hour, e42dep)
  })

  test_that("means_by_group, weighting", {
    w <- "weight"
    means_by_group(efc, c12hour, e42dep, weights = weight)
    means_by_group(efc, c12hour, e42dep, weights = "weight")
    means_by_group(efc, c12hour, e42dep, weights = w)
  })

  test_that("means_by_group, grouping", {
    means_by_group(efc_grouped, c12hour, e42dep)
  })

  test_that("means_by_group, grouped weighting", {
    w <- "weight"
    means_by_group(efc_grouped, c12hour, e42dep, weights = weight)
    means_by_group(efc_grouped, c12hour, e42dep, weights = "weight")
    means_by_group(efc_grouped, c12hour, e42dep, weights = w)
  })

}
