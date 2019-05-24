#' @title Rescale design weights for multilevel analysis
#' @name scale_weights
#'
#' @description Most functions to fit multilevel and mixed effects models only
#'    allow to specify frequency weights, but not design (i.e. sampling or probability)
#'    weights, which should be used when analyzing complex samples and survey data.
#'    \code{scale_weights()} implements an algorithm proposed by Aaparouhov (2006)
#'    and Carle (2009) to rescale design weights in survey data to account for
#'    the grouping structure of multilevel models, which then can be used for
#'    multilevel modelling.
#'
#' @param x A data frame.
#' @param cluster.id Variable indicating the grouping structure (strata) of
#'    the survey data (level-2-cluster variable).
#' @param pweight Variable indicating the probability (design or sampling)
#'    weights of the survey data (level-1-weight).
#'
#' @return \code{x}, with two new variables: \code{svywght_a} and \code{svywght_b},
#'    which represent the rescaled design weights to use in multilevel models
#'    (use these variables for the \code{weights} argument).
#'
#' @details Rescaling is based on two methods: For \code{svywght_a}, the sample
#'    weights \code{pweight} are adjusted by a factor that represents the proportion
#'    of cluster size divided by the sum of sampling weights within each cluster.
#'    The adjustment factor for \code{svywght_b} is the sum of sample weights
#'    within each cluster devided by the sum of squared sample weights within
#'    each cluster (see Carle (2009), Appendix B).
#'    \cr \cr
#'    Regarding the choice between scaling methods A and B, Carle suggests
#'    that "analysts who wish to discuss point estimates should report results
#'    based on weighting method A. For analysts more interested in residual
#'    between-cluster variance, method B may generally provide the least biased
#'    estimates". In general, it is recommended to fit a non-weighted model
#'    and weighted models with both scaling methods and when comparing the
#'    models, see whether the "inferential decisions converge", to gain
#'    confidence in the results.
#'    \cr \cr
#'    Though the bias of scaled weights decreases with increasing cluster size,
#'    method A is preferred when insufficient or low cluster size is a concern.
#'    \cr \cr
#'    The cluster ID and probably PSU may be used as random effects (e.g.
#'    nested design, or cluster and PSU as varying intercepts), depending
#'    on the survey design that should be mimicked.
#'
#' @references Carle AC. \emph{Fitting multilevel models in complex survey data with design weights: Recommendations} BMC Medical Research Methodology 2009, 9(49): 1-13
#'    \cr \cr
#'    Asparouhov T. \emph{General Multi-Level Modeling with Sampling Weights} Communications in Statistics - Theory and Methods 2006, 35: 439-460
#'
#' @examples
#' data(nhanes_sample)
#' scale_weights(nhanes_sample, SDMVSTRA, WTINT2YR)
#'
#' library(lme4)
#' nhanes_sample <- scale_weights(nhanes_sample, SDMVSTRA, WTINT2YR)
#' glmer(
#'   total ~ factor(RIAGENDR) * (log(age) + factor(RIDRETH1)) + (1 | SDMVPSU),
#'   family = poisson(),
#'   data = nhanes_sample,
#'   weights = svywght_a
#' )
#'
#' @importFrom dplyr group_by summarise n right_join enquo filter quo_name slice n
#' @importFrom rlang .data
#' @importFrom sjmisc is_empty
#' @export
scale_weights <- function(x, cluster.id, pweight) {

  # quote cluster.id and get name as string

  quo.id <- dplyr::enquo(cluster.id)
  id.name <- dplyr::quo_name(quo.id)


  # quote sample weights and get name as string

  quo.weights <- dplyr::enquo(pweight)
  pw.name <- dplyr::quo_name(quo.weights)


  # check if weight has missings. we need to remove them first,
  # and add back weights to correct cases later

  weight_missings <- which(is.na(x[[pw.name]]))
  weight_non_na <- which(!is.na(x[[pw.name]]))

  if (!sjmisc::is_empty(weight_missings)) {
    dummy_x <- dplyr::slice(x, !! weight_non_na)
  } else {
    dummy_x <- x
  }


  # copy data set, so we only append the two new weights

  tmp <- dummy_x
  tmp$s_q_w <- dummy_x[[pw.name]]^2


  # compute sum of weights per cluster

  tmp <- tmp %>%
    dplyr::group_by(!! quo.id) %>%
    dplyr::summarise(
      sum_wij = sum(!! quo.weights),
      sum_sqwij = sum(.data$s_q_w),
      nj = dplyr::n()
    ) %>%
    dplyr::right_join(dummy_x, by = id.name)


  # multiply the original weight by the fraction of the
  # sampling unit total population based on Carle 2009

  w_a <- tmp[[pw.name]] * tmp$nj / tmp$sum_wij
  w_b <- tmp[[pw.name]] * tmp$sum_wij / tmp$sum_sqwij

  x$svywght_a <- NA
  x$svywght_b <- NA

  x$svywght_a[weight_non_na] <- w_a
  x$svywght_b[weight_non_na] <- w_b

  x
}

