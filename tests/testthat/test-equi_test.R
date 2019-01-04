.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {

  if (suppressWarnings(
        require("testthat") &&
        require("sjstats") &&
        require("rstanarm") &&
        require("lme4") &&
        require("brms") &&
        require("sjmisc")
    )) {
    context("sjstats, stan")


    # fit linear model
    data(sleepstudy)
    sleepstudy$age <- round(runif(nrow(sleepstudy), min = 20, max = 60))
    sleepstudy$Rdicho <- dicho(sleepstudy$Reaction)

    m <- stan_glmer(
      Reaction ~ Days + age + (1 | Subject),
      data = sleepstudy, QR = TRUE,
      # this next line is only to keep the example small in size!
      chains = 2, cores = 1, seed = 12345, iter = 500
    )

    m2 <- stan_glmer(
      Rdicho ~ Days + age + (1 | Subject),
      data = sleepstudy, QR = TRUE,
      family = binomial,
      chains = 2, iter = 500
    )


    m3 <- brm(
      Reaction ~ Days + age + (1 | Subject),
      data = sleepstudy, chains = 2, iter = 500
    )

    m4 <- brm(
      Rdicho ~ Days + age + (1 | Subject),
      data = sleepstudy, chains = 2, iter = 500,
      family = bernoulli()
    )


    expect_string <- function(object) {
      # 2. Call expect()
      expect_is(
        object$term,
        "character"
      )
    }


    test_that("equi_test", {
      expect_string(equi_test(m))
      expect_string(equi_test(m, rope = c(-1, 1)))
      expect_string(equi_test(m, rope = c(-1, 1), eff_size = .5))
      expect_string(equi_test(m, eff_size = .5))

      expect_string(equi_test(m2))
      expect_string(equi_test(m2, rope = c(-1, 1)))
      expect_string(equi_test(m2, rope = c(-1, 1), eff_size = .5))
      expect_string(equi_test(m2, eff_size = .5))

      p <- equi_test(m, out = "plot")
      p <- equi_test(m, rope = c(-1, 1), out = "plot")
      p <- equi_test(m, rope = c(-1, 1), eff_size = .5, out = "plot")
      p <- equi_test(m, eff_size = .5, out = "plot")

      p <- equi_test(m2, out = "plot")
      p <- equi_test(m2, rope = c(-1, 1), out = "plot")
      p <- equi_test(m2, rope = c(-1, 1), eff_size = .5, out = "plot")
      p <- equi_test(m2, eff_size = .5, out = "plot")

      p <- equi_test(m, out = "viewer")
      p <- equi_test(m, rope = c(-1, 1), out = "viewer")
      p <- equi_test(m, rope = c(-1, 1), eff_size = .5, out = "viewer")
      p <- equi_test(m, eff_size = .5, out = "viewer")

      p <- equi_test(m2, out = "viewer")
      p <- equi_test(m2, rope = c(-1, 1), out = "viewer")
      p <- equi_test(m2, rope = c(-1, 1), eff_size = .5, out = "viewer")
      p <- equi_test(m2, eff_size = .5, out = "viewer")
    })


    test_that("equi_test", {
      expect_string(equi_test(m3))
      expect_string(equi_test(m3, rope = c(-1, 1)))
      expect_string(equi_test(m3, rope = c(-1, 1), eff_size = .5))
      expect_string(equi_test(m3, eff_size = .5))

      expect_string(equi_test(m4))
      expect_string(equi_test(m4, rope = c(-1, 1)))
      expect_string(equi_test(m4, rope = c(-1, 1), eff_size = .5))
      expect_string(equi_test(m4, eff_size = .5))
    })

  }
}
