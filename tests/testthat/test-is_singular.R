.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {
  if (require("testthat") && require("sjstats") && require("glmmTMB") && require("lme4")) {
    context("sjstats, is_singular")

    # fit linear model
    data(fish)

    m <- glmmTMB(
      count ~ child + camper + (1 | persons),
      ziformula = ~ child + camper + (1 | persons),
      data = fish,
      family = truncated_poisson()
    )

    m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

    test_that("is_singular", {
      expect_false(is_singular(m2))
      expect_false(is_singular(m))
    })

    test_that("converge_ok", {
      expect_true(converge_ok(m2))
    })
  }
}
