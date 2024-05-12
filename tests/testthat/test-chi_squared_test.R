skip_if_not_installed("effectsize")
test_that("chi_squared_test", {
  data(efc)
  set.seed(123)
  efc$weight <- abs(rnorm(nrow(efc), 1, 0.3))
  chi_squared_test(efc, "c161sex", by = "e16sex")
  chi_squared_test(efc, "c161sex", by = "e16sex", weights = "weight")
  chi_squared_test(efc, "c161sex", probabilities = c(0.3, 0.7))

})
