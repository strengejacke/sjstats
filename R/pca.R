#' @title Tidy summary of Principal Component Analysis
#' @name pca
#' @description ...
#'
#' @param x A data frame or a \code{\link[stats]{prcomp}} object.
#' @param nf Number of components to extract. If \code{rotation = "varmiax"}
#'    and \code{nf = NULL}, number of components is based on the Kaiser-criteria.
#' @param rotation Rotation of the factor loadings. May be \code{"varimax"} for
#'    orthogonal rotation or \code{"oblimin"} for oblique transformation.
#'
#' @return A tidy data frame with either all loadings of principal components
#'    (for \code{pca()}) or a rotated loadings matrix (for \code{pca_rotate()}).
#'
#' @details The \code{print()}-method for \code{pca_rotate()} has a
#'    \code{cutoff}-argument, which is a scalar between 0 and 1, indicating
#'    which (absolute) values from the loadings should be blank in the
#'    output. By default, all loadings below .1 (or -.1) are not shown.
#'
#' @examples
#' data(efc)
#' # recveive first item of COPE-index scale
#' start <- which(colnames(efc) == "c82cop1")
#' # recveive last item of COPE-index scale
#' end <- which(colnames(efc) == "c90cop9")
#'
#' # extract principal components
#' pca(efc[, start:end])
#'
#' # extract principal components, varimax-rotation.
#' # number of components based on Kaiser-criteria
#' pca_rotate(efc[, start:end])
#'
#' @importFrom stats prcomp na.omit
#' @importFrom tibble tibble
#' @importFrom sjmisc rotate_df
#' @importFrom rlang .data
#' @export
pca <- function(x) {

  if (!inherits(x, c("prcomp", "data.frame")))
    stop("`x` must be of class `prcomp` or a data frame.", call. = F)


  # if x is a df, run prcomp

  if (inherits(x, "data.frame"))
    x <- stats::prcomp(stats::na.omit(x), retx = TRUE, center = TRUE, scale. = TRUE)


  # get tidy summary of prcomp object

  tmp <- tibble::tibble(
    comp = sprintf("PC%i", seq_len(length(x$sdev))),
    std.dev = x$sdev,
    eigen = .data$std.dev ^ 2,
    prop.var = .data$eigen / sum(.data$eigen),
    cum.var = cumsum(.data$prop.var)
  )

  # add information on Kaiser criteria and loadings
  attr(tmp, "kaiser") <- which(tmp$eigen < 1)[1] - 1
  attr(tmp, "loadings") <- x$rotation %*% diag(x$sdev)


  # rotate df for proper output
  tmp <- sjmisc::rotate_df(tmp, cn = T)

  # add class-attribute for printing
  class(tmp) <- c("sj_pca", class(tmp))

  tmp
}


#' @rdname pca
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export
pca_rotate <- function(x, nf = NULL, rotation = c("varimax", "oblimin")) {

  rotation <- match.arg(rotation)

  if (!inherits(x, c("prcomp", "data.frame")))
    stop("`x` must be of class `prcomp` or a data frame.", call. = F)

  if (!inherits(x, "data.frame") && rotation == "oblimin")
    stop("`x` must be a data frame for oblimin-rotation.", call. = F)

  if (rotation == "oblimin" && !requireNamespace("psych", quietly = TRUE))
    stop("Package `psych` required for `oblimin`-rotation.", call. = F)


  # rotate loadings

  if (rotation == "oblimin")
    tmp <- psych::principal(r = x, nfactors = nf, rotate = "oblimin")
  else {
    if (!inherits(x, "sj_pca")) x <- pca(x)

    loadings <- attr(x, "loadings", exact = TRUE)
    if (is.null(nf)) nf <- attr(x, "kaiser", exact = TRUE)

    tmp <- stats::varimax(loadings[, seq_len(nf)])
  }


  # tweak column names and class attributes

  tmp <- as.data.frame(unclass(tmp$loadings))
  colnames(tmp) <- sprintf("PC%i", 1:ncol(tmp))
  class(tmp) <- c("sj_pca_rotate", "data.frame")


  # add explained proportions and proportional and cumulative variance

  attr(tmp, "variance") <- tibble::tibble(
    prop.var = colSums(tmp ^ 2) / nrow(tmp),
    cum.var = cumsum(.data$prop.var),
    prop.exp = .data$prop.var / sum(.data$prop.var),
    cum.exp = cumsum(.data$prop.exp)
  )

  tmp
}

