.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {
  if (require("testthat") && require("sjstats") && require("lme4")) {
    context("sjstats, get_re_var")

    # load sample data
    data(sleepstudy)

    # fit linear mixed model
    m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)

    set.seed(2018)
    sleepstudy$mygrp <- sample(1:45, size = 180, replace = TRUE)
    m2 <- lmer(Reaction ~ Days + (1 | mygrp) + (1 | Subject), sleepstudy)

    test_that("get_re_var", {
      get_re_var(m, "tau.00")
      get_re_var(m2, "tau.00")
      get_re_var(icc(m), "tau.00")
      get_re_var(icc(m2), "tau.00")
      get_re_var(m, "rho.01")
      get_re_var(m2, "rho.01")
      get_re_var(icc(m), "rho.01")
      get_re_var(icc(m2), "rho.01")
      get_re_var(m, "sigma_2")
      get_re_var(m2, "sigma_2")
      get_re_var(icc(m), "sigma_2")
      get_re_var(icc(m2), "sigma_2")
    })
  }
}
