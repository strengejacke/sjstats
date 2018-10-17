#' @rdname hdi
#' @export
equi_test <- function(x, ...) {
  UseMethod("equi_test")
}


#' @export
equi_test.default <- function(x, rope, eff_size, out = c("txt", "viewer", "browser", "plot"), ...) {
  out <- match.arg(out)
  equi_test_worker(x = x, rope = rope, eff_size = eff_size, out = out, fm = NULL, ...)
}


#' @rdname hdi
#' @export
equi_test.stanreg <- function(x, rope, eff_size, out = c("txt", "viewer", "browser", "plot"), ...) {
  out <- match.arg(out)
  equi_test_worker(x = x, rope = rope, eff_size = eff_size, out = out, fm = model_family(x), ...)
}


#' @rdname hdi
#' @export
equi_test.brmsfit <- function(x, rope, eff_size, out = c("txt", "viewer", "browser", "plot"), ...) {
  # check for pkg availability, else function might fail
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Please install and load package `brms` first.")

  out <- match.arg(out)
  equi_test_worker(x = x, rope = rope, eff_size = eff_size, out = out, fm = model_family(x), ...)
}


#' @export
equi_test.stanfit <- function(x, rope, eff_size, out = c("txt", "viewer", "browser", "plot"), ...) {
  out <- match.arg(out)
  equi_test_worker(x = x, rope = rope, eff_size = eff_size, out = out, fm = model_family(x), ...)
}


#' @export
equi_test.data.frame <- function(x, rope, eff_size, out = c("txt", "viewer", "browser", "plot"), ...) {
  out <- match.arg(out)
  equi_test_worker(x = x, rope = rope, eff_size = eff_size, out = out, fm = NULL, ...)
}


#' @importFrom purrr map_df
#' @importFrom sjmisc add_columns var_rename is_empty add_variables
#' @importFrom dplyr case_when select pull
#' @importFrom stats sd
#' @importFrom bayesplot neff_ratio
equi_test_worker <- function(x, rope, eff_size, out, fm, ...) {

  if (fm$is_multivariate)
    stop("Multivariate response models not supported yet.", call. = F)


  if (out != "txt" && !requireNamespace("sjPlot", quietly = TRUE)) {
    message("Package `sjPlot` needs to be loaded to print HTML tables.")
    out <- "txt"
  }


  dat <- as.data.frame(x, stringsAsFactors = FALSE)
  if (missing(eff_size)) eff_size <- .1

  if (!is.null(fm)) {
    if (missing(rope)) {
      if (fm$is_linear)
        eff_range <- stats::sd(resp_val(x)) * eff_size
      else if (fm$is_bin)
        eff_range <- stats::sd(dat[[1]]) * eff_size / 4
      else
        eff_range <- stats::sd(dat[[1]]) * eff_size

      rope <- c(-1, 1) * eff_range
    }
  }

  .hdi <- hdi(x = x, prob = .95, trans = NULL, type = "fixed")
  .rope <- rope(x = x, rope = rope, trans = NULL, type = "fixed")
  .neff <- nrow(dat)

  result <- dplyr::case_when(
    .hdi$hdi.low > rope[2] ~ "reject",
    .hdi$hdi.high < rope[1] ~ "reject",
    .hdi$hdi.low >= rope[1] & .hdi$hdi.high <= rope[2] ~ "accept",
    TRUE ~ "undecided"
  )

  # for convenience reasons, also add proportion of values outside rope
  dat <- .hdi %>%
    dplyr::select(-1) %>%
    sjmisc::add_columns(.rope) %>%
    dplyr::select(-3) %>%
    sjmisc::add_variables(decision = result, .after = 1) %>%
    sjmisc::var_rename(rope = "inside.rope")

  # indicate parameters with critical number of effective samples

  critical <- NULL
  if (inherits(x, c("stanfit", "stanreg", "brmsfit"))) {
    nratio <- bayesplot::neff_ratio(x)
    nratio <- nratio[names(nratio) %in% dat$term]
    critical <- which(nratio < .7)
    if (!sjmisc::is_empty(critical))
      dat$term[critical] <- sprintf("%s (*)", dat$term[critical])
  }

  attr(dat, "rope") <- rope
  attr(dat, "eff_size") <- eff_size
  attr(dat, "nsamples") <- .neff
  attr(dat, "critical") <- !sjmisc::is_empty(critical)

  if (is_stan_model(x)) {
    dat <- remove_effects_from_stan(dat, type = "fixed", is.brms = inherits(x, "brmsfit"))
  }

  # save how to print output
  attr(dat, "print") <- out

  if (out == "plot") {
    plot_sj_equi_test(dat, x, ...)
  } else if (out == "txt") {
    class(dat) <- c("sj_equi_test", class(dat))
    dat
  } else {
    class(dat) <- c("sjt_equi_test", class(dat))
    dat
  }
}


#' @importFrom dplyr slice select case_when
#' @importFrom purrr map2_df
#' @importFrom tidyr gather
#' @importFrom sjlabelled as_factor get_dv_labels
plot_sj_equi_test <- function(x, model, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE) && !requireNamespace("ggridges", quietly = TRUE)) {
    warning("Packages 'ggplot2' and 'ggridges' required to plot test for practical equivalence.", call. = FALSE)
    return(x)
  }

  remove <- c(1, string_contains("sigma", x$term))

  # if we have intercept-only models, keep at least the intercept
  if (length(remove) == nrow(x)) remove <- remove[-1]

  x <- dplyr::slice(x, -!! remove)

  # remove indicator for insufficient sample chains
  x$term <- gsub(" (*)", "", x = x$term, fixed = TRUE)

  tmp <- model %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::select(!! x$term) %>%
    purrr::map2_df(
      x$hdi.low,
      function(.x, .y) {
        .x[.x < .y] <- NA
        .x
      }) %>%
    purrr::map2_df(
      x$hdi.high,
      function(.x, .y) {
        .x[.x > .y] <- NA
        .x
      }) %>%
    tidyr::gather(
      key = "predictor",
      value = "estimate"
    ) %>%
    sjlabelled::as_factor(.data$predictor)


  x$decision <- dplyr::case_when(
    x$decision == "accept" ~ "reject",
    x$decision == "reject" ~ "accept",
    TRUE ~ "undecided"
  )

  tmp$grp <- NA
  for (i in 1:nrow(x))
    tmp$grp[tmp$predictor == x$term[i]] <- x$decision[i]

  rope <- attr(x, "rope")

  tmp$predictor <- factor(tmp$predictor)
  tmp$predictor <- factor(tmp$predictor, levels = rev(levels(tmp$predictor)))


  # check for user defined arguments

  fill.color <- c("#00b159", "#d11141", "#cccccc")
  rope.color <- "#004D80"
  rope.alpha <- 0.15
  x.title <- "95% Highest Density Region of Posterior Samples"
  legend.title <- "Decision on Parameters"
  labels <- levels(tmp$predictor)
  names(labels) <- labels

  fill.color <- fill.color[sort(unique(match(x$decision, c("accept", "reject", "undecided"))))]

  add.args <- lapply(match.call(expand.dots = F)$`...`, function(x) x)
  if ("colors" %in% names(add.args)) fill.color <- eval(add.args[["colors"]])
  if ("x.title" %in% names(add.args)) x.title <- eval(add.args[["x.title"]])
  if ("rope.color" %in% names(add.args)) rope.color <- eval(add.args[["rope.color"]])
  if ("rope.alpha" %in% names(add.args)) rope.alpha <- eval(add.args[["rope.alpha"]])
  if ("legend.title" %in% names(add.args)) legend.title <- eval(add.args[["legend.title"]])
  if ("labels" %in% names(add.args)) labels <- eval(add.args[["labels"]])

  rope.line.alpha <- 1.25 * rope.alpha
  if (rope.line.alpha > 1) rope.line.alpha <- 1


  ggplot2::ggplot(tmp, ggplot2::aes_string(x = "estimate", y = "predictor", fill = "grp")) +
    ggplot2::annotate("rect", xmin = rope[1], xmax = rope[2], ymin = 0, ymax = Inf, fill = rope.color, alpha = rope.alpha) +
    ggplot2::geom_vline(xintercept = 0, colour = rope.color, size = .8, alpha = rope.line.alpha) +
    ggridges::geom_density_ridges2(rel_min_height = 0.01, scale = 2, alpha = .5) +
    ggplot2::scale_fill_manual(values = fill.color) +
    ggplot2::labs(x = x.title, y = NULL, fill = legend.title) +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::theme(legend.position = "bottom")
}
