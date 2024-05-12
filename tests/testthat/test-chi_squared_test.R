skip_if_not_installed("effectsize")
skip_if_not_installed("datawizard")

test_that("chi_squared_test", {
  data(efc)
  set.seed(123)
  efc$weight <- abs(rnorm(nrow(efc), 1, 0.3))
  out1 <- chi_squared_test(efc, "c161sex", by = "e16sex")
  out2 <- chisq.test(efc$c161sex, efc$e16sex)
  expect_equal(out1$statistic, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)

  out <- chi_squared_test(efc, "c161sex", by = "e16sex", weights = "weight")
  expect_equal(out$statistic, 2.415755, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out$effect_size, 0.05448519, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out$p, 0.1201201, tolerance = 1e-4, ignore_attr = TRUE)

  out1 <- chi_squared_test(efc, "c161sex", probabilities = c(0.3, 0.7))
  out2 <- chisq.test(table(efc$c161sex), p = c(0.3, 0.7))
  expect_equal(out1$statistic, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)

  out <- chi_squared_test(efc, "c161sex", probabilities = c(0.3, 0.7), weights = "weight")
  expect_equal(out$statistic, 20.07379, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out$effect_size, 0.0974456, tolerance = 1e-4, ignore_attr = TRUE)

  set.seed(1234)
  d <- data.frame(
    survey_1 = sample(c("Approve", "Disapprove"), size = 1000, replace = TRUE, prob = c(0.45, 0.55)),
    survey_2 = sample(c("Approve", "Disapprove"), size = 1000, replace = TRUE, prob = c(0.42, 0.58))
  )
  out1 <- chi_squared_test(d, "survey_1", "survey_2", paired = TRUE)
  out2 <- mcnemar.test(table(d))
  expect_equal(out1$statistic, out2$statistic, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$p, out2$p.value, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(out1$effect_size, 0.03170437, tolerance = 1e-4, ignore_attr = TRUE)
})
