skip_if_not_installed("survey")
skip_if_not_installed("datawizard")

test_that("kruskal_wallis_test", {
  data(efc)
  set.seed(123)
  efc$weight <- abs(rnorm(nrow(efc), 1, 0.3))
  out1 <- kruskal_wallis_test(efc, "e17age", by = "c172code")
  out2 <- kruskal.test(e17age ~ c172code, data = efc)
  expect_equal(out1$Chi2, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_snapshot(print(out1))

  set.seed(123)
  wide_data <- data.frame(
    scale1 = runif(20),
    scale2 = runif(20),
    scale3 = runif(20)
  )
  long_data <- data.frame(
    scales = c(wide_data$scale1, wide_data$scale2, wide_data$scale3),
    groups = as.factor(rep(c("A", "B", "C"), each = 20)),
    stringsAsFactors = FALSE
  )
  out1 <- kruskal_wallis_test(wide_data, select = c("scale1", "scale2", "scale3"))
  out2 <- kruskal_wallis_test(long_data, select = "scales", by = "groups")
  out3 <- kruskal.test(scales ~ groups, data = long_data)
  expect_equal(out1$Chi2, out2$Chi2, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$Chi2, out3$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out3$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_snapshot(print(out1))

  out1 <- kruskal_wallis_test(efc, "e17age", by = "c172code", weights = "weight")
})
