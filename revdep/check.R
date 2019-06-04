library("devtools")

revdep_check()
revdep_check_save_summary()
revdep_check_print_problems()

# tools::check_packages_in_dir(
#   dir = "C:/Users/mail/Documents/R/sj",
#   check_args = c("--as-cran", ""),
#   reverse = list(repos = getOption("repos")["CRAN"])
# )
