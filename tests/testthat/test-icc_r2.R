.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {
  if (require("testthat") &&
      require("sjstats") &&
      require("lme4") &&
      require("sjmisc") &&
      require("dplyr") &&
      require("MASS") &&
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

    data("Salamanders")

    m <- glmmTMB(
      count ~ mined + (1 + sample | site),
      zi =  ~ mined,
      family = poisson,
      data = Salamanders
    )

    test_that("r2", {
      expect_warning(r2(m))
    })

    set.seed(2)
    n <- 10
    x1 <- mvrnorm(mu = rep(0, n), Sigma = .7 ^ as.matrix(dist(1:n)))
    y1 <- x1 + rnorm(n)
    x2 <- mvrnorm(mu = rep(0, n), Sigma = .7 ^ as.matrix(dist(1:n)))
    y2 <- x2 + rnorm(n)
    x3 <- mvrnorm(mu = rep(0, n), Sigma = .7 ^ as.matrix(dist(1:n)))
    y3 <- x3 + rnorm(n)
    x4 <- mvrnorm(mu = rep(0, n), Sigma = .7 ^ as.matrix(dist(1:n)))
    y4 <- x4 + rnorm(n)
    x5 <- mvrnorm(mu = rep(0, n), Sigma = .7 ^ as.matrix(dist(1:n)))
    y5 <- x5 + rnorm(n)

    x <- c(x1, x2, x3, x4, x5)
    y <- c(y1, y2, y3, y4, y5)
    times <- factor(c(1:n, 1:n, 1:n, 1:n, 1:n))
    group <- factor(c(rep(1, n), rep(2, n), rep(3, n), rep(4, n), rep(5, n)))

    dat0 <- data.frame(y, times, group)

    j11 <- glmmTMB(y ~ times + (times | group), data = dat0)
    j12 <- glmmTMB(y ~ times + ar1(times + 0 | group) , data = dat0)
    j13 <- glmmTMB(y ~ times + group +  ar1(times + 0 | group) , data = dat0)

    test_that("r2", {
      r2(j11)
      expect_warning(r2(j12))
      expect_warning(r2(j13))
    })
  }
}
