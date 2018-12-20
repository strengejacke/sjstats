if (require("testthat") && require("sjstats") && require("splines")) {
  context("sjstats, model_frame")

  data(efc)

  m1 <- lm(neg_c_7 ~ e42dep + ns(c160age), data = efc)
  m2 <- lm(neg_c_7 ~ e42dep + ns(c160age, knots = 2), data = efc)
  m3 <- lm(neg_c_7 ~ e42dep + bs(c160age, degree = 3), data = efc)
  m4 <- lm(neg_c_7 ~ e42dep + bs(c160age, degree = 1), data = efc)

  m5 <- lm(neg_c_7 ~ e42dep + c160age, data = efc)


  test_that("model_frame", {
    mf1 <- model_frame(m1)
    mf2 <- model_frame(m2)
    mf3 <- model_frame(m3)
    mf4 <- model_frame(m4)
    mf5 <- model.frame(m5)

    expect_equal(as.vector(mf1$c160age), as.vector(mf5$c160age))
    expect_equal(as.vector(mf2$c160age), as.vector(mf5$c160age))
    expect_equal(as.vector(mf3$c160age), as.vector(mf5$c160age))
    expect_equal(as.vector(mf4$c160age), as.vector(mf5$c160age))
  })
}

