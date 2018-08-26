library(testthat)
library(sjstats)

if (length(strsplit(packageDescription("sjstats")$Version, "\\.")[[1]]) > 3) {
  Sys.setenv("RunAllsjstatsTests" = "yes")
}

test_check("sjstats")
