.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest && Sys.getenv("USER") != "travis") {

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


    test_that("tidy_stan", {
      expect_string(tidy_stan(m))
      expect_string(tidy_stan(m, prob = c(.5, .8)))
      expect_string(tidy_stan(m, type = "random"))
      expect_string(tidy_stan(m, type = "all"))
      expect_string(tidy_stan(m, prob = c(.5, .8), type = "random"))
      expect_string(tidy_stan(m, prob = c(.5, .8), type = "all"))
      expect_string(tidy_stan(m, prob = c(.5, .8), typical = "mean"))
      expect_string(tidy_stan(m, prob = c(.5, .8), type = "random", typical = "mean"))
      expect_string(tidy_stan(m, prob = c(.5, .8), type = "all", typical = "mean"))

      expect_string(tidy_stan(m2))
      expect_string(tidy_stan(m2, prob = c(.5, .8)))
      expect_string(tidy_stan(m2, type = "random"))
      expect_string(tidy_stan(m2, type = "all"))

      expect_string(tidy_stan(m2, trans = "exp"))
      expect_string(tidy_stan(m2, prob = c(.5, .8), trans = "exp"))
      expect_string(tidy_stan(m2, type = "random", trans = "exp"))
      expect_string(tidy_stan(m2, type = "all", trans = "exp"))
    })

    test_that("tidy_stan", {
      expect_string(tidy_stan(m3))
      expect_string(tidy_stan(m3, prob = c(.5, .8)))
      expect_string(tidy_stan(m3, type = "random"))
      expect_string(tidy_stan(m3, type = "all"))
      expect_string(tidy_stan(m3, prob = c(.5, .8), type = "random"))
      expect_string(tidy_stan(m3, prob = c(.5, .8), type = "all"))
      expect_string(tidy_stan(m3, prob = c(.5, .8), typical = "mean"))
      expect_string(tidy_stan(m3, prob = c(.5, .8), type = "random", typical = "mean"))
      expect_string(tidy_stan(m3, prob = c(.5, .8), type = "all", typical = "mean"))

      expect_string(tidy_stan(m4))
      expect_string(tidy_stan(m4, prob = c(.5, .8)))
      expect_string(tidy_stan(m4, type = "random"))
      expect_string(tidy_stan(m4, type = "all"))

      expect_string(tidy_stan(m4, trans = "exp"))
      expect_string(tidy_stan(m4, prob = c(.5, .8), trans = "exp"))
      expect_string(tidy_stan(m4, type = "random", trans = "exp"))
      expect_string(tidy_stan(m4, type = "all", trans = "exp"))
    })
  }
}
