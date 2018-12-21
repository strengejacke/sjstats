.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {

  if (suppressWarnings(
        require("testthat") &&
        require("sjstats") &&
        require("rstanarm") &&
        require("lme4") &&
        require("sjmisc") &&
        require("brms") &&
        require("dplyr")
    )) {
    context("sjstats, icc-bayes")


    # fit linear model
    data(sleepstudy)

    set.seed(123)
    sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
    sleepstudy <- sleepstudy %>%
      group_by(mygrp) %>%
      mutate(mysubgrp = sample(1:30, size = n(), replace = TRUE))

    m1 <- brms::brm(
      Reaction ~ Days + (1 + Days | Subject),
      data = sleepstudy,
      iter = 500,
      chains = 1
    )

    m2 <- brms::brm(
      Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
      data = sleepstudy,
      iter = 500,
      chains = 1
    )

    m3 <- brms::brm(
      Reaction ~ Days + (1 | mysubgrp) + (1 | Subject),
      data = sleepstudy,
      iter = 500,
      chains = 1
    )

    m4 <- brms::brm(
      Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
      data = sleepstudy,
      iter = 500,
      chains = 1
    )

    test_that("icc-bayes", {
      icc(m1, adjusted = TRUE)
      icc(m1, adjusted = FALSE)
      icc(m2, adjusted = TRUE)
      icc(m2, adjusted = FALSE)
      icc(m3, adjusted = TRUE)
      icc(m3, adjusted = FALSE)
      icc(m4, adjusted = TRUE)
      icc(m4, adjusted = FALSE)
    })



    data(sleepstudy)
    sleepstudy$age <- round(runif(nrow(sleepstudy), min = 20, max = 60))
    sleepstudy$Rdicho <- dicho(sleepstudy$Reaction, as.num = TRUE)

    m1 <- stan_glmer(
      Rdicho ~ Days + age + (1 | Subject),
      data = sleepstudy,
      family = binomial,
      chains = 2, iter = 500
    )

    m2 <- brm(
      Rdicho ~ Days + age + (1 | Subject),
      data = sleepstudy,
      family = binomial,
      chains = 2, iter = 500
    )

    test_that("icc-bayes", {
      icc(m1, adjusted = TRUE)
      icc(m1, adjusted = FALSE)
      icc(m2, adjusted = TRUE)
      icc(m2, adjusted = FALSE)
    })
  }
}
