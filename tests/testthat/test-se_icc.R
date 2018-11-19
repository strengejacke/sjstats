.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {
  if (require("testthat") && require("sjstats") && require("lme4")) {
    context("sjstats, se_icc")

    # load sample data
    data(sleepstudy)

    # fit linear mixed model
    m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)

    set.seed(2018)
    sleepstudy$mygrp <- sample(1:45, size = 180, replace = TRUE)
    m2 <- lmer(Reaction ~ Days + (1 | mygrp) + (1 | Subject), sleepstudy)

    test_that("se_icc", {
      sjstats::se(sjstats::icc(m), nsim = 50)
      sjstats::se(sjstats::icc(m2), nsim = 50)
      sjstats::se(sjstats::icc(m, adjusted = TRUE), nsim = 50)
      sjstats::se(sjstats::icc(m2, adjusted = TRUE), nsim = 50)
    })
  }
}
