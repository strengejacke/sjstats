# Help-functions

data_frame <- function(...) {
  x <- data.frame(..., stringsAsFactors = FALSE)
  rownames(x) <- NULL
  x
}


is_merMod <- function(fit) {
  inherits(fit, c("lmerMod", "glmerMod", "nlmerMod", "merModLmerTest"))
}


.compact_character <- function(x) {
  x[!sapply(x, function(i) is.null(i) || !nzchar(i, keepNA = TRUE) || is.na(i) || any(i == "NULL", na.rm = TRUE))]
}


.format_symbols <- function(x) {
  if (.unicode_symbols()) {
    x <- gsub("Delta", "\u0394", x, ignore.case = TRUE)
    x <- gsub("Phi", "\u03D5", x, ignore.case = TRUE)
    x <- gsub("Eta", "\u03B7", x, ignore.case = TRUE)
    x <- gsub("Epsilon", "\u03b5", x, ignore.case = TRUE)
    x <- gsub("Omega", "\u03b5", x, ignore.case = TRUE)
    x <- gsub("R2", "R\u00b2", x, ignore.case = TRUE)
    x <- gsub("Chi2", "\u03C7\u00b2", x, ignore.case = TRUE)
    x <- gsub("Chi-squared", "\u03C7\u00b2", x, ignore.case = TRUE)
    x <- gsub("Chi", "\u03C7", x, ignore.case = TRUE)
    x <- gsub("Sigma", "\u03C3", x, ignore.case = TRUE)
    x <- gsub("Rho", "\u03C1", x, ignore.case = TRUE)
    x <- gsub("Mu", "\u03BC", x, ignore.case = TRUE)
    x <- gsub("Theta", "\u03B8", x, ignore.case = TRUE)
    x <- gsub("Fei", "\u05E4\u200E", x, ignore.case = TRUE)
  }
  x
}


.unicode_symbols <- function() {
  win_os <- tryCatch(
    {
      si <- Sys.info()
      if (is.null(si["sysname"])) {
        FALSE
      } else {
        si["sysname"] == "Windows" || startsWith(R.version$os, "mingw")
      }
    },
    error = function(e) {
      TRUE
    }
  )
  l10n_info()[["UTF-8"]] && ((win_os && getRversion() >= "4.2") || (!win_os && getRversion() >= "4.0"))
}


.is_pseudo_numeric <- function(x) {
  (is.character(x) && !anyNA(suppressWarnings(as.numeric(stats::na.omit(x[nzchar(x, keepNA = TRUE)]))))) || (is.factor(x) && !anyNA(suppressWarnings(as.numeric(levels(x))))) # nolint
}


.misspelled_string <- function(source, searchterm, default_message = NULL) {
  if (is.null(searchterm) || length(searchterm) < 1) {
    return(default_message)
  }
  # used for many matches
  more_found <- ""
  # init default
  msg <- ""
  # remove matching strings
  same <- intersect(source, searchterm)
  searchterm <- setdiff(searchterm, same)
  source <- setdiff(source, same)
  # guess the misspelled string
  possible_strings <- unlist(lapply(searchterm, function(s) {
    source[.fuzzy_grep(source, s)] # nolint
  }), use.names = FALSE)
  if (length(possible_strings)) {
    msg <- "Did you mean "
    if (length(possible_strings) > 1) {
      # make sure we don't print dozens of alternatives for larger data frames
      if (length(possible_strings) > 5) {
        more_found <- sprintf(
          " We even found %i more possible matches, not shown here.",
          length(possible_strings) - 5
        )
        possible_strings <- possible_strings[1:5]
      }
      msg <- paste0(msg, "one of ", toString(paste0("\"", possible_strings, "\"")))
    } else {
      msg <- paste0(msg, "\"", possible_strings, "\"")
    }
    msg <- paste0(msg, "?", more_found)
  } else {
    msg <- default_message
  }
  # no double white space
  insight::trim_ws(msg)
}


.fuzzy_grep <- function(x, pattern, precision = NULL) {
  if (is.null(precision)) {
    precision <- round(nchar(pattern) / 3)
  }
  if (precision > nchar(pattern)) {
    return(NULL)
  }
  p <- sprintf("(%s){~%i}", pattern, precision)
  grep(pattern = p, x = x, ignore.case = FALSE)
}
