.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {

  if (require("testthat") && require("sjstats")) {
    # fit linear model
    data(efc)
    m <- aov(
      c12hour ~ as.factor(e42dep) + as.factor(c172code) + c160age,
      data = efc
    )

    test_that("eta_sq", {
      expect_equal(eta_sq(m, partial = FALSE)$etasq, c(0.266, 0.005, 0.048), tolerance = 1e-2)
      expect_equal(eta_sq(m, partial = TRUE)$partial.etasq, c(0.281, 0.008, 0.066), tolerance = 1e-2)
      eta_sq(m, partial = FALSE, ci.lvl = .5)
      eta_sq(m, partial = TRUE, ci.lvl = .6)
    })

    test_that("omega_sq", {
      omega_sq(m, partial = FALSE)
      omega_sq(m, partial = TRUE)
      omega_sq(m, partial = FALSE, ci.lvl = .5)
      omega_sq(m, partial = TRUE, ci.lvl = .6)
    })

    test_that("cohens_f", {
      cohens_f(m)
    })

    test_that("anova_stats", {
      anova_stats(m, digits = 3)
      anova_stats(m, digits = 5)
      anova_stats(car::Anova(m, type = 2))
      anova_stats(car::Anova(m, type = 3))
    })


    set.seed(123)
    fit <- aov(
      c12hour ~ as.factor(e42dep) + as.factor(c172code) + c160age,
      data = efc
    )

    test_that("omega_sq", {
      omega_sq(fit, partial = TRUE, ci.lvl = 0.95)
    })

    set.seed(123)
    data(mtcars)
    m <-
      stats::aov(
        formula = mpg ~ wt + qsec + Error(disp / am),
        data = mtcars
      )

    test_that("anova_stats", {
      anova_stats(m, digits = 3)
      anova_stats(m, digits = 5)
      eta_sq(m, partial = TRUE, ci.lvl = 0.95)
      eta_sq(m, partial = FALSE, ci.lvl = 0.95)
      omega_sq(m, partial = TRUE, ci.lvl = 0.95)
      omega_sq(m, partial = FALSE, ci.lvl = 0.95)
    })
  }

}
