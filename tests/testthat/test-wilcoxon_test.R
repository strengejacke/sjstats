skip_if_not_installed("survey")
skip_if_not_installed("datawizard")
skip_if_not_installed("coin")

test_that("wilcoxon_test", {
  data(mtcars)
  out1 <- wilcoxon_test(mtcars, "mpg")
  out2 <- suppressWarnings(wilcox.test(mtcars$mpg ~ 1))
  expect_equal(out1$v, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_snapshot(print(out1))

  out1 <- wilcoxon_test(mtcars, c("mpg", "hp"))
  out2 <- suppressWarnings(wilcox.test(mtcars$mpg, mtcars$hp, paired = TRUE))
  expect_equal(out1$v, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_snapshot(print(out1))

  data(iris)
  d <- iris[iris$Species != "setosa", ]
  out <- wilcoxon_test(d, "Sepal.Width", by = "Species")
  expect_snapshot(print(out))
})
