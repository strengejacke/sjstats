skip_if_not_installed("datawizard")
skip_if_not_installed("effectsize")

test_that("t_test", {
  data(efc)
  set.seed(123)
  efc$weight <- abs(rnorm(nrow(efc), 1, 0.3))
  expect_snapshot(t_test(efc, "e17age"))
  expect_snapshot(t_test(efc, "e17age", "e16sex"))
  expect_snapshot(t_test(efc, c("e17age", "c160age")))
  expect_snapshot(t_test(efc, c("e17age", "c160age"), paired = TRUE))

  expect_snapshot(t_test(efc, "e17age", weights = "weight"))
  expect_snapshot(t_test(efc, "e17age", "e16sex", weights = "weight"))
  expect_snapshot(t_test(efc, c("e17age", "c160age"), weights = "weight"))
  expect_snapshot(t_test(efc, c("e17age", "c160age"), weights = "weight", paired = TRUE))
})
