.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {
  if (require("testthat") && require("sjstats") && require("glmmTMB")) {
    context("sjstats, icc")

    # fit linear model
    data(fish)

    m <- glmmTMB(
      count ~ child + camper + (1 | persons),
      ziformula = ~ child + camper + (1 | persons),
      data = fish,
      family = truncated_poisson()
    )

    test_that("icc", {
      expect_warning(r2(m))
      expect_warning(icc(m))
      expect_warning(icc(m, adjusted = TRUE))
    })
  }
}
