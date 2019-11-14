#' @title Mann-Whitney-U-Test
#' @name mwu
#' @description This function performs a Mann-Whitney-U-Test (or Wilcoxon rank sum test,
#'                see \code{\link[stats]{wilcox.test}} and \code{\link[coin]{wilcox_test}})
#'                for \code{x}, for each group indicated by \code{grp}. If \code{grp}
#'                has more than two categories, a comparison between each combination of
#'                two groups is performed. \cr \cr
#'                The function reports U, p and Z-values as well as effect size r
#'                and group-rank-means.
#'
#' @param x Bare (unquoted) variable name, or a character vector with the variable name.
#' @param distribution Indicates how the null distribution of the test statistic should be computed.
#'          May be one of \code{"exact"}, \code{"approximate"} or \code{"asymptotic"}
#'          (default). See \code{\link[coin]{wilcox_test}} for details.
#'
#' @inheritParams wtd_se
#' @inheritParams grpmean
#'
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
#' mwu(efc, e17age, e42dep)
#'
#' @importFrom stats na.omit wilcox.test kruskal.test
#' @importFrom sjmisc recode_to is_empty
#' @importFrom sjlabelled get_labels as_numeric
#' @importFrom rlang quo_name enquo
#' @export
mwu <- function(data,
                x,
                grp,
                distribution = "asymptotic",
                out = c("txt", "viewer", "browser"),
                encoding = "UTF-8",
                file = NULL) {

  out <- match.arg(out)

  if (out != "txt" && !requireNamespace("sjPlot", quietly = TRUE)) {
    message("Package `sjPlot` needs to be loaded to print HTML tables.")
    out <- "txt"
  }

  if (!requireNamespace("coin", quietly = TRUE)) {
    stop("Package `coin` needs to be installed to compute the Mann-Whitney-U test.", call. = FALSE)
  }


  # create quosures
  grp.name <- rlang::quo_name(rlang::enquo(grp))
  dv.name <- rlang::quo_name(rlang::enquo(x))

  # create string with variable names
  vars <- c(grp.name, dv.name)

  # get data
  data <- suppressMessages(dplyr::select(data, !! vars))

  grp <- data[[grp.name]]
  dv <- data[[dv.name]]

  # coerce factor and character to numeric
  if (is.factor(grp) || is.character(grp)) grp <- sjlabelled::as_numeric(grp)

  # group "counter" (index) should start with 1, not 0
  if (min(grp, na.rm = TRUE) < 1) grp <- sjmisc::recode_to(grp, lowest = 1, append = FALSE)

  # retrieve unique group values. need to iterate all values
  grp_values <- sort(unique(stats::na.omit(grp)))

  # length of value range
  cnt <- length(grp_values)
  labels <- sjlabelled::get_labels(
    grp, attr.only = F, values = NULL, non.labelled = T
  )

  df <- data.frame()
  for (i in seq_len(cnt)) {
    for (j in i:cnt) {
      if (i != j) {
        # retrieve cases (rows) of subgroups
        xsub <- dv[which(grp == grp_values[i] | grp == grp_values[j])]
        ysub <- grp[which(grp == grp_values[i] | grp == grp_values[j])]

        # this is for unpaired wilcox.test()
        xsub_2 <- stats::na.omit(dv[which(grp == grp_values[i])])
        ysub_2 <- stats::na.omit(dv[which(grp == grp_values[j])])

        # only use rows with non-missings
        ysub <- ysub[which(!is.na(xsub))]

        # remove missings
        xsub <- as.numeric(stats::na.omit(xsub))
        ysub.n <- stats::na.omit(ysub)

        # grouping variable is a factor
        ysub <- as.factor(ysub.n)

        wcdat <- data.frame(
          x = xsub,
          y = ysub
        )

        # perfom wilcox test
        wt <- coin::wilcox_test(x ~ y, data = wcdat, distribution = distribution)

        # compute statistics
        u <- as.numeric(coin::statistic(wt, type = "linear"))
        z <- as.numeric(coin::statistic(wt, type = "standardized"))
        p <- coin::pvalue(wt)
        r <- abs(z / sqrt(length(ysub)))
        w <- stats::wilcox.test(xsub_2, ysub_2, paired = FALSE)$statistic
        rkm.i <- mean(rank(xsub)[which(ysub.n == grp_values[i])], na.rm = TRUE)
        rkm.j <- mean(rank(xsub)[which(ysub.n == grp_values[j])], na.rm = TRUE)

        # compute n for each group
        n_grp1 <- length(xsub[which(ysub.n == grp_values[i])])
        n_grp2 <- length(xsub[which(ysub.n == grp_values[j])])

        # generate result data frame
        df <-
          rbind(
            df,
            cbind(
              grp1 = grp_values[i],
              grp1.label = labels[i],
              grp1.n = n_grp1,
              grp2 = grp_values[j],
              grp2.label = labels[j],
              grp2.n = n_grp2,
              u = u,
              w = w,
              p = p,
              z = z,
              r = r,
              rank.mean.grp1 = rkm.i,
              rank.mean.grp2 = rkm.j
            )
          )
      }
    }
  }

  # convert variables
  df[["grp1"]] <- as.numeric(as.character(df[["grp1"]]))
  df[["grp2"]] <- as.numeric(as.character(df[["grp2"]]))
  df[["grp1.n"]] <- as.numeric(as.character(df[["grp1.n"]]))
  df[["grp2.n"]] <- as.numeric(as.character(df[["grp2.n"]]))
  df[["grp1.label"]] <- as.character(df[["grp1.label"]])
  df[["grp2.label"]] <- as.character(df[["grp2.label"]])
  df[["u"]] <- as.numeric(as.character(df[["u"]]))
  df[["w"]] <- as.numeric(as.character(df[["w"]]))
  df[["p"]] <- as.numeric(as.character(df[["p"]]))
  df[["z"]] <- as.numeric(as.character(df[["z"]]))
  df[["r"]] <- as.numeric(as.character(df[["r"]]))
  df[["rank.mean.grp1"]] <- as.numeric(as.character(df[["rank.mean.grp1"]]))
  df[["rank.mean.grp2"]] <- as.numeric(as.character(df[["rank.mean.grp2"]]))

  # prepare a data frame that can be used for 'sjt.df'.
  tab.df <-
    data_frame(
      Groups = sprintf("%s<br>%s", df$grp1.label, df$grp2.label),
      N = sprintf("%s<br>%s", df$grp1.n, df$grp2.n),
      'Mean Rank' = sprintf("%.2f<br>%.2f", df$rank.mean.grp1, df$rank.mean.grp2),
      'Mann-Whitney-U' = as.character(df$u),
      'Wilcoxon-W' = as.character(df$w),
      Z = sprintf("%.3f", df$z),
      'Effect Size' = sprintf("%.3f", df$r),
      p = sprintf("%.3f", df$p)
    )

  # replace 0.001 with <0.001
  tab.df$p[which(tab.df$p == "0.001")] <- "<0.001"

  ret.df <- list(df = df, tab.df = tab.df, data = data.frame(dv, grp))

  # save how to print output
  attr(ret.df, "print") <- out
  attr(ret.df, "encoding") <- encoding
  attr(ret.df, "file") <- file

  if (out %in% c("viewer", "browser"))
    class(ret.df) <- c("mwu", "sjt_mwu")
  else
    class(ret.df) <- c("mwu", "sj_mwu")

  ret.df
}
