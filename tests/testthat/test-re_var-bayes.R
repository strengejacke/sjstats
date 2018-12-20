.runThisTest <- Sys.getenv("RunAllsjstatsTests") == "yes"

if (.runThisTest) {
  if (suppressWarnings(
    require("testthat") &&
    require("sjstats") &&
    require("rstanarm") &&
    require("lme4") &&
    require("sjmisc")
  )) {
    context("sjstats, get_re_var")

    # load sample data
    data(sleepstudy)
    sleepstudy$age <- round(runif(nrow(sleepstudy), min = 20, max = 60))
    sleepstudy$Rdicho <- dicho(sleepstudy$Reaction)

    m1 <- stan_glmer(
      Reaction ~ Days + age + (1 | Subject),
      data = sleepstudy, QR = TRUE,
      chains = 2, iter = 500
    )

    m2 <- stan_glmer(
      Rdicho ~ Days + age + (1 | Subject),
      data = sleepstudy, QR = TRUE,
      family = binomial,
      chains = 2, iter = 500
    )

    m3 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
    m4 <- glmer(RDicho ~ Days + (Days | Subject), data = sleepstudy, family = binomial)

    test_that("re_var", {
      re_var(m1)
      re_var(m2)
      re_var(m3)
      re_var(m4)
      re_var(m1, adjusted = TRUE)
      re_var(m2, adjusted = TRUE)
      re_var(m3, adjusted = TRUE)
      re_var(m4, adjusted = TRUE)
    })
  }
}
