.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {
  if (require("testthat") &&
      require("sjstats") &&
      require("lme4") &&
      require("sjmisc") &&
      require("dplyr") &&
      require("glmmTMB")) {

    context("sjstats, icc")

    # fit linear model
    data(fish)
    data(sleepstudy)
    data(efc)
    efc$e15relat <- as.factor(efc$e15relat)

    set.seed(123)
    sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
    sleepstudy <- sleepstudy %>%
      dplyr::group_by(mygrp) %>%
      dplyr::mutate(mysubgrp = sample(1:30, size = n(), replace = TRUE))

    m1 <- glmmTMB(
      count ~ child + camper + (1 | persons),
      ziformula = ~ child + camper + (1 | persons),
      data = fish,
      family = truncated_poisson()
    )

    m2 <- lmer(
      Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
      data = sleepstudy
    )

    m3 <- lmer(
      Reaction ~ Days + (1 | mysubgrp) + (1 | Subject),
      data = sleepstudy
    )

    m4 <- glmer(
      tot_sc_e ~ e42dep + c160age + (1 | e15relat),
      data = data,
      family = poisson()
    )

    test_that("icc", {
      expect_warning(r2(m1))
      expect_warning(icc(m1))
      expect_warning(icc(m1, adjusted = TRUE))
    })

    test_that("icc", {
      icc(m2, adjusted = TRUE)
      icc(m2, adjusted = FALSE)
      icc(m3, adjusted = TRUE)
      icc(m3, adjusted = FALSE)
      icc(m4, adjusted = TRUE)
      icc(m4, adjusted = FALSE)
    })
  }
}
