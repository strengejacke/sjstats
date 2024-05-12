skip_if_not_installed("survey")
skip_if_not_installed("datawizard")

test_that("kruskal_wallis_test", {
#' data(efc)
#' # Kruskal-Wallis test for elder's age by education
#' kruskal_wallis_test(efc, "e17age", by = "c172code")
#'
#' # when data is in wide-format, specify all relevant continuous
#' # variables in `select` and omit `by`
#' set.seed(123)
#' wide_data <- data.frame(
#'   scale1 = runif(20),
#'   scale2 = runif(20),
#'   scale3 = runif(20)
#' )
#' kruskal_wallis_test(wide_data, select = c("scale1", "scale2", "scale3"))
#'
#' # same as if we had data in long format, with grouping variable
#' long_data <- data.frame(
#'   scales = c(wide_data$scale1, wide_data$scale2, wide_data$scale3),
#'   groups = rep(c("A", "B", "C"), each = 20)
#' )
#' kruskal_wallis_test(long_data, select = "scales", by = "groups")
#' # base R equivalent
#' kruskal.test(scales ~ groups, data = long_data)
  data(efc)
  set.seed(123)
  efc$weight <- abs(rnorm(nrow(efc), 1, 0.3))
  out1 <- mann_whitney_test(efc, "e17age", by = "e16sex")
  out2 <- wilcox.test(e17age ~ e16sex, data = efc)
  expect_equal(out1$w, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$estimate, -1561, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$r, 0.2571254, tolerance = 1e-4, ignore_attr = TRUE)

  set.seed(123)
  wide_data <- data.frame(scale1 = runif(20), scale2 = runif(20))
  out1 <- mann_whitney_test(wide_data, select = c("scale1", "scale2"))
  out2 <- wilcox.test(wide_data$scale1, wide_data$scale2)
  expect_equal(out1$w, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$r, 0.05132394, tolerance = 1e-4, ignore_attr = TRUE)

  out <- mann_whitney_test(efc, "e17age", by = "e16sex", weights = "weight")
  expect_equal(out$p, 1.976729e-14, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out$estimate, 0.1594972, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out$r, 0.2599877, tolerance = 1e-4, ignore_attr = TRUE)
})
