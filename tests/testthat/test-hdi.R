.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {

  if (suppressWarnings(
        require("testthat") &&
        require("sjstats") &&
        require("rstanarm") &&
        require("lme4") &&
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


    expect_string <- function(object) {
      # 2. Call expect()
      expect_is(
        object$term,
        "character"
      )
    }


    test_that("hdi", {
      expect_string(hdi(m))
      expect_string(hdi(m2))

      expect_string(hdi(m, prob = c(.3, .7)))
      expect_string(hdi(m2, prob = c(.3, .7)))

      expect_string(hdi(m, prob = c(.3, .7), type = "all"))
      expect_string(hdi(m2, prob = c(.3, .7), type = "all"))

      expect_string(hdi(m2, prob = c(.3, .7), type = "all", trans = "exp"))
    })
  }
}
