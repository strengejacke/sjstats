# Help-functions


is_merMod <- function(fit) {
  return(any(class(fit) %in% c("lmerMod", "glmerMod", "nlmerMod", "merModLmerTest")))
}

