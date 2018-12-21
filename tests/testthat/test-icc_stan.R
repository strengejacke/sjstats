.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {

  if (suppressWarnings(
    require("testthat") &&
    require("sjstats") &&
    require("rstanarm") &&
    require("lme4") &&
    require("sjmisc") &&
    require("dplyr")
  )) {
    context("sjstats, icc-stan")


    # fit linear model
    data(sleepstudy)
    data(efc)
    efc$e15relat <- as.factor(efc$e15relat)

    set.seed(123)
    sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
    sleepstudy <- sleepstudy %>%
      dplyr::group_by(mygrp) %>%
      dplyr::mutate(mysubgrp = sample(1:30, size = n(), replace = TRUE))

    m1 <- rstanarm::stan_glmer(
      Reaction ~ Days + (1 + Days | Subject),
      data = sleepstudy,
      iter = 500,
      chains = 1
    )

    m2 <- rstanarm::stan_glmer(
      Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
      data = sleepstudy,
      iter = 500,
      chains = 1
    )

    m3 <- rstanarm::stan_glmer(
      Reaction ~ Days + (1 | mysubgrp) + (1 | Subject),
      data = sleepstudy,
      iter = 500,
      chains = 1
    )

    m4 <- rstanarm::stan_glmer(
      Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
      data = sleepstudy,
      iter = 500,
      chains = 1
    )

    m5 <- rstanarm::stan_glmer(
      tot_sc_e ~ e42dep + c160age + (1 | e15relat),
      data = efc,
      family = poisson(),
      iter = 500,
      chains = 1
    )

    test_that("icc-stan", {
      icc(m1, adjusted = TRUE)
      icc(m1, adjusted = FALSE)
      icc(m1, ppd = TRUE, adjusted = FALSE)
      icc(m2, adjusted = TRUE)
      icc(m2, adjusted = FALSE)
      icc(m2, ppd = TRUE, adjusted = FALSE)
      icc(m3, adjusted = TRUE)
      icc(m3, adjusted = FALSE)
      icc(m3, ppd = TRUE, adjusted = FALSE)
      icc(m4, adjusted = TRUE)
      icc(m4, adjusted = FALSE)
      icc(m4, ppd = TRUE, adjusted = FALSE)
    })
  }
}
