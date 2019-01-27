#' @rdname pred_vars
#' @importFrom stats model.frame formula getCall na.omit
#' @importFrom purrr map_lgl map reduce
#' @importFrom dplyr select bind_cols full_join
#' @importFrom sjmisc add_columns is_empty
#' @export
model_frame <- function(x, fe.only = TRUE) {
  # we may store model weights here later
  mw <- NULL

  tryCatch(
    {
      if (inherits(x, "stanmvreg"))
        fitfram <- suppressMessages(
          purrr::reduce(stats::model.frame(x), ~ dplyr::full_join(.x, .y))
        )
      else if (inherits(x, "clm2"))
        fitfram <- x$location
      else if (inherits(x, c("merMod", "lmerMod", "glmerMod", "nlmerMod", "merModLmerTest")))
        fitfram <- stats::model.frame(x, fixed.only = fe.only)
      else if (inherits(x, "lme"))
        fitfram <- x$data
      else if (inherits(x, "vgam"))
        fitfram <- get(x@misc$dataname, envir = parent.frame())
      else if (inherits(x, c("gee", "gls")))
        fitfram <- eval(x$call$data, envir = parent.frame())
      else if (inherits(x, "Zelig-relogit"))
        fitfram <- get_zelig_relogit_frame(x)
      else if (inherits(x, "gmnl"))
        fitfram <- x$mf
      else if (inherits(x, "vglm")) {
        if (!length(x@model)) {
          env <- environment(x@terms$terms)
          if (is.null(env)) env <- parent.frame()
          fcall <- x@call
          fcall$method <- "model.frame"
          fcall$smart <- FALSE
          fitfram <- eval(fcall, env, parent.frame())
        } else {
          fitfram <- x@model
        }
      } else if (inherits(x, "MixMod")) {
        fitfram <- x$model_frames$mfX
        if (!sjmisc::is_empty(x$model_frames$mfZ))
          fitfram <- sjmisc::add_columns(x$model_frames$mfZ, fitfram, replace = TRUE)
        if (!sjmisc::is_empty(x$model_frames$mfX_zi))
          fitfram <- sjmisc::add_columns(x$model_frames$mfX_zi, fitfram, replace = TRUE)
        if (!sjmisc::is_empty(x$model_frames$mfZ_zi))
          fitfram <- sjmisc::add_columns(x$model_frames$mfZ_zi, fitfram, replace = TRUE)

        fitfram <- sjmisc::add_variables(fitfram, grp__id = x$id)
        colnames(fitfram)[ncol(fitfram)] <- x$id_name[1]
      } else if (inherits(x, "MCMCglmm")) {
        env_dataframes <- names(which(unlist(eapply(.GlobalEnv, is.data.frame))))
        env_dataframes <- names(which(unlist(eapply(.GlobalEnv, is.data.frame))))
        pv <- pred_vars(x, fe.only = FALSE)
        matchframe <- unlist(lapply(env_dataframes, function(.x) {
          dat <- get(.x)
          all(pv %in% colnames(dat))
        }))
        mf <- env_dataframes[matchframe][1]
        if (!is.na(mf))
          fitfram <- get(mf)
        else
          fitfram <- NULL
      } else
        fitfram <- stats::model.frame(x)
    },
    error = function(x) { fitfram <- NULL }
  )


  if (is.null(fitfram)) {
    warning("Could not get model frame.", call. = F)
    return(NULL)
  }


  # do we have an offset, not specified in the formula?

  if ("(offset)" %in% colnames(fitfram)) {
    if (obj_has_name(x, "call")) {
      if (obj_has_name(x$call, "offset")) {
        offcol <- which(colnames(fitfram) == "(offset)")
        colnames(fitfram)[offcol] <- var_names(deparse(x$call$offset, width.cutoff = 500L))
      }
    }
  }


  # clean 1-dimensional matrices

  fitfram <- purrr::modify_if(fitfram, is.matrix, function(x) {
    if (dim(x)[2] == 1 && !inherits(x, c("ns", "bs")))
      as.vector(x)
    else
      x
  })


  # check if we have any matrix columns, e.g. from splines

  mc <- purrr::map_lgl(fitfram, is.matrix)


  # don't change response value, if it's a matrix
  # bound with cbind()
  rn <- resp_var(x, combine = TRUE)
  trials.data <- NULL

  if (mc[1] && rn == colnames(fitfram)[1]) {
    mc[1] <- FALSE
    tryCatch(
      {
        trials.data <- as.data.frame(fitfram[[1]])
        colnames(trials.data) <- resp_var(x, combine = FALSE)
      },
      error = function(x) { NULL }
    )
  }


  # if we have any matrix columns, we remove them from original
  # model frame and convert them to regular data frames, give
  # proper column names and bind them back to the original model frame

  if (any(mc)) {
    # try to get model data from environment
    md <- tryCatch(
      {
        eval(stats::getCall(x)$data, environment(stats::formula(x)))
      },
      error = function(x) { NULL }
    )

    # if data not found in environment, reduce matrix variables into regular vectors
    if (is.null(md)) {
      # first, we select the non-matrix variables. calling "as_tibble" would
      # remove their column name, so we us as_tibble to convert matrix
      # to vectors only for the matrix-columns
      fitfram_matrix <- dplyr::select(fitfram, which(mc))
      fitfram_nonmatrix <- dplyr::select(fitfram, -which(mc))
      fitfram_matrix <- dplyr::bind_cols(purrr::map(fitfram_matrix, ~ as.data.frame(.x, stringsAsFactors = FALSE)))
      fitfram <- dplyr::bind_cols(fitfram_nonmatrix, fitfram_matrix)
    } else {

      # fix NA in column names

      if (any(is.na(colnames(md)))) colnames(md) <- make.names(colnames(md))

      # get "matrix" terms and "normal" predictors, but exclude
      # response variable(s)

      fitfram_matrix <- dplyr::select(fitfram, -which(mc))
      spline.term <- get_vn_helper(names(which(mc)))
      other.terms <- get_vn_helper(colnames(fitfram_matrix))[-1]

      # now we have all variable names that we need from the original
      # data set

      needed.vars <- c(other.terms, spline.term)

      # if response is a matrix vector (e.g. multivariate response),
      # we need to include all response names as well, because else
      # rows may not match due to additional missings in the response variables

      if (is.matrix(fitfram[[1]])) {
        needed.vars <- c(dimnames(fitfram[[1]])[[2]], needed.vars)
      } else {
        needed.vars <- c(colnames(fitfram)[1], needed.vars)
      }

      # check model weights

      if ("(weights)" %in% needed.vars && !obj_has_name(md, "(weights)")) {
        needed.vars <- needed.vars[-which(needed.vars == "(weights)")]
        mw <- fitfram[["(weights)"]]
      }


      if (inherits(x, "coxph")) {
        fitfram <- md
      } else {
        fitfram <- stats::na.omit(dplyr::select(md, !! needed.vars))
      }

      # add back model weights, if any
      if (!is.null(mw)) fitfram$`(weights)` <- mw
    }

    # check if we really have all formula terms in our model frame now
    pv <- tryCatch(
      {
        pred_vars(x, fe.only = fe.only)
      },
      error = function(x) { NULL }
    )

    if (!is.null(pv) && !all(pv %in% colnames(fitfram))) {
      warning("Some model terms could not be found in model frame. You probably need to load the data into the environment.", call. = FALSE)
    }

  }

  # check if we have monotonic variables, included in formula
  # with "mo()"? If yes, remove from model frame
  mos_eisly <- grepl(pattern = "^mo\\(([^,)]*).*", x = colnames(fitfram))
  if (any(mos_eisly)) fitfram <- fitfram[!mos_eisly]

  # clean variable names
  cvn <- get_vn_helper(colnames(fitfram))

  # do we have duplicated names?
  dupes <- which(duplicated(cvn))
  if (!sjmisc::is_empty(dupes)) cvn[dupes] <- sprintf("%s.%s", cvn[dupes], 1:length(dupes))

  colnames(fitfram) <- cvn

  # add back possible trials-data
  if (!is.null(trials.data)) {
    new.cols <- setdiff(colnames(trials.data), colnames(fitfram))
    if (!sjmisc::is_empty(new.cols)) fitfram <- cbind(fitfram, trials.data[, new.cols, drop = FALSE])
  }


  # for glmmtmb, check dispersion and zi-formula
  # and add variables to model frame

  if (inherits(x, "glmmTMB")) {
    disp <- tryCatch(
      {all.vars(x$modelInfo$allForm$dispformula[[2L]])},
      error = function(x) { NULL}
    )

    if (!is.null(disp)) {
      fitfram <- tryCatch(
        {
          eval(x$call$data, envir = parent.frame()) %>%
            dplyr::select(!! disp) %>%
            sjmisc::add_columns(fitfram, replace = TRUE)
        },
        error = function(x) { fitfram }
      )
    }


    zi <- tryCatch(
      {all.vars(x$modelInfo$allForm$ziformula[[2L]])},
      error = function(x) { NULL}
    )

    if (!is.null(zi)) {
      fitfram <- tryCatch(
        {
          eval(x$call$data, envir = parent.frame()) %>%
            dplyr::select(!! zi) %>%
            sjmisc::add_columns(fitfram, replace = TRUE)
        },
        error = function(x) { fitfram }
      )
    }
  }


  fitfram
}

#' @importFrom dplyr select
get_zelig_relogit_frame <- function(x) {
  vars <- c(resp_var(x), pred_vars(x))
  dplyr::select(x$data, !! vars)
}
