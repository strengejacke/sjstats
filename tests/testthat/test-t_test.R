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

  out1 <- t_test(efc, "e17age")
  out2 <- t.test(efc$e17age ~ 1)
  expect_equal(out1$statistic, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$effect_size, 9.774916, tolerance = 1e-4, ignore_attr = TRUE)

  out1 <- t_test(efc, "e17age", "e16sex")
  out2 <- t.test(efc$e17age ~ efc$e16sex)
  expect_equal(out1$statistic, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$effect_size, -0.5641989, tolerance = 1e-4, ignore_attr = TRUE)
})
