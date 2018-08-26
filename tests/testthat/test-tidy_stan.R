.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {

  stopifnot(require("testthat"),
            require("sjstats"),
            require("rstanarm"),
            require("lme4"),
            require("sjmisc"))

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


  test_that("hdi", {
    expect_string(hdi(m))
    expect_string(hdi(m2))

    expect_string(hdi(m, prob = c(.3, .7)))
    expect_string(hdi(m2, prob = c(.3, .7)))

    expect_string(hdi(m, prob = c(.3, .7), type = "all"))
    expect_string(hdi(m2, prob = c(.3, .7), type = "all"))

    expect_string(hdi(m2, prob = c(.3, .7), type = "all", trans = "exp"))
  })


  test_that("mcse", {
    expect_string(mcse(m))
    expect_string(mcse(m2))
    expect_string(mcse(m, type = "all"))
    expect_string(mcse(m2, type = "all"))
  })


  test_that("n_eff", {
    expect_string(n_eff(m))
    expect_string(n_eff(m2))
    expect_string(n_eff(m, type = "all"))
    expect_string(n_eff(m2, type = "all"))
  })

}
