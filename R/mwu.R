#' @title Mann-Whitney-U-Test
#' @name mwu
#' @description This function performs a Mann-Whitney-U-Test (or \code{Wilcoxon rank sum test},
#'                see \code{\link{wilcox.test}} and \code{\link[coin]{wilcox_test}})
#'                for \code{x}, for each group indicated by \code{grp}. If \code{grp}
#'                has more than two categories, a comparison between each combination of
#'                two groups is performed. \cr \cr
#'                The function reports U, p and Z-values as well as effect size r
#'                and group-rank-means.
#'
#' @param x Numeric vector / variable, where the Mann-Whitney-U-Test should be applied to
#' @param grp Grouping variable indicating the groups that should be used for comparison
#' @param distribution Indicates how the null distribution of the test statistic should be computed.
#'          May be one of \code{"exact"}, \code{"approximate"} or \code{"asymptotic"}
#'          (default). See \code{\link[coin]{wilcox_test}} for details.
#' @param weights Integer valued weights for the observations. By default,
#'          this is \code{NULL}.
#' @return (Invisibly) returns a data frame with U, p and Z-values for each group-comparison
#'         as well as effect-size r; additionally, group-labels and groups' n's are
#'         also included.
#'
#' @note This function calls the \code{\link[coin]{wilcox_test}} with formula. If \code{grp}
#'         has more than two groups, additionally a Kruskal-Wallis-Test (see \code{\link{kruskal.test}})
#'         is performed. \cr \cr
#'         Interpretation of effect sizes, as a rule-of-thumb:
#'         \itemize{
#'          \item small effect >= 0.1
#'          \item medium effect >= 0.3
#'          \item large effect >= 0.5
#'        }
#'
#' @examples
#' data(efc)
#' # Mann-Whitney-U-Tests for elder's age by elder's dependency.
#' mwu(efc$e17age, efc$e42dep)
#'
#' @importFrom stats na.omit wilcox.test kruskal.test
#' @importFrom coin wilcox_test pvalue statistic
#' @export
mwu <- function(x, grp, distribution = "asymptotic", weights = NULL) {
  # coerce factor and character to numeric
  if (is.factor(grp) || is.character(grp)) grp <- sjmisc::to_value(grp)
  # group "counter" (index) should start with 1, not 0
  if (min(grp, na.rm = TRUE) == 0) grp <- grp + 1
  # retrieve unique group values. need to iterate all values
  grp_values <- sort(unique(stats::na.omit(grp)))
  # length of value range
  cnt <- length(grp_values)
  labels <- sjmisc::get_labels(grp, attr.only = F, include.values = NULL,
                               include.non.labelled = T)
  message("Performing Mann-Whitney-U-Test...")
  message("---------------------------------")
  message("showing statistics between groups (x|y)")
  df <- data.frame()
  for (i in 1:cnt) {
    for (j in i:cnt) {
      if (i != j) {
        # retrieve cases (rows) of subgroups
        xsub <- x[which(grp == grp_values[i] | grp == grp_values[j])]
        ysub <- grp[which(grp == grp_values[i] | grp == grp_values[j])]
        # only use rows with non-missings
        ysub <- ysub[which(!is.na(xsub))]
        # adjust weights, pick rows from subgroups (see above)
        if (!is.null(weights)) {
          wsub <- as.integer(stats::na.omit(weights[which(!is.na(xsub))]))
        }
        # remove missings
        xsub <- as.numeric(stats::na.omit(xsub))
        ysub.n <- stats::na.omit(ysub)
        # grouping variable is a factor
        ysub <- as.factor(ysub.n)
        if (is.null(weights)) {
          wt <- coin::wilcox_test(xsub ~ ysub, distribution = distribution)
        } else {
          wt <- coin::wilcox_test(xsub ~ ysub,
                                  distribution = distribution,
                                  weights = as.formula("~wsub"))
        }
        # compute statistics
        u <- as.numeric(coin::statistic(wt, type = "linear"))
        z <- as.numeric(coin::statistic(wt, type = "standardized"))
        p <- coin::pvalue(wt)
        r <- abs(z / sqrt(length(x)))
        w <- stats::wilcox.test(xsub, ysub.n, paired = TRUE)$statistic
        rkm.i <- mean(rank(xsub)[which(ysub.n == grp_values[i])], na.rm = TRUE)
        rkm.j <- mean(rank(xsub)[which(ysub.n == grp_values[j])], na.rm = TRUE)
        # compute n for each group
        n_grp1 <- length(xsub[which(ysub.n == grp_values[i])])
        n_grp2 <- length(xsub[which(ysub.n == grp_values[j])])
        # print to console
        if (is.null(labels)) {
          cat(sprintf("Groups (%i|%i), n = %i/%i:\n", grp_values[i],
                      grp_values[j], n_grp1, n_grp2))
        } else {
          cat(sprintf("Groups %i = %s (n = %i) | %i = %s (n = %i):\n",
                      grp_values[i], labels[i], n_grp1, grp_values[j],
                      labels[j], n_grp2))
        }
        if (p < 0.001) {
          p <- 0.001
          p.string <- "<"
        } else {
          p.string <- "="
        }
        cat(sprintf("  U = %.3f, W = %.3f, p %s %.3f, Z = %.3f\n  effect-size r = %.3f\n  rank-mean(%i) = %.2f\n  rank-mean(%i) = %.2f\n\n", u, w, p.string, p, z, r, i, rkm.i, j, rkm.j))
        df <- rbind(df, cbind(grp1 = grp_values[i], grp1.label = labels[i],
                              grp1.n = n_grp1, grp2 = grp_values[j], grp2.label = labels[j],
                              grp2.n = n_grp2, u = u, w = w, p = p, z = z, r = r,
                              rank.mean.grp1 = rkm.i, rank.mean.grp2 = rkm.j))
      }
    }
  }
  # if we have more than 2 groups, also perfom kruskal-wallis-test
  if (cnt > 2) {
    message("Performing Kruskal-Wallis-Test...")
    message("---------------------------------")
    kw <- stats::kruskal.test(x, grp)
    cat(sprintf("chi-squared = %.3f\n", kw$statistic))
    cat(sprintf("df = %i\n", kw$parameter))
    if (kw$p.value < 0.001) {
      p  <- 0.001
      p.string <- "<"
    } else {
      p <- kw$p.value
      p.string <- "="
    }
    cat(sprintf("p %s %.3f\n", p.string, p))
  }
  # prepare a data frame that can be used for 'sjt.df'.
  tab.df <- data.frame(Groups = sprintf("%s<br>%s", df$grp1.label, df$grp2.label),
                       N = sprintf("%s<br>%s", df$grp1.n, df$grp2.n),
                       'Mean Rank' = sprintf("%.2f<br>%.2f",
                                             as.numeric(as.character(df$rank.mean.grp1)),
                                             as.numeric(as.character(df$rank.mean.grp2))),
                       'Mann-Whitney-U' = df$u,
                       'Wilcoxon-W' = df$w,
                       Z = sprintf("%.3f", as.numeric(as.character(df$z))),
                       'Effect Size' = sprintf("%.3f", as.numeric(as.character(df$r))),
                       p = sprintf("%.3f", as.numeric(as.character(df$p))))
  # replace 0.001 with <0.001
  levels(tab.df$p)[which(levels(tab.df$p) == "0.001")] <- "<0.001"
  # return both data frames
  invisible(structure(class = "mwu",list(df = df, tab.df = tab.df)))
}
