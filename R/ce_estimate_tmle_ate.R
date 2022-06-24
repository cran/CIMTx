#' Causal inference with multiple treatments using TMLE for ATE effects
#'
#' The function \code{ce_estimate_tmle_ate} implements
#' TMLE to estimate ATE effect with
#' multiple treatments using observational data.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param w A numeric vector representing the treatment groups.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param sl_library A character vector of prediction algorithms.
#' A list of functions included in the SuperLearner package
#' can be found with \code{\link[SuperLearner:listWrappers]{listWrappers}}.
#' @param ... Other parameters that can be passed through to functions.
#'
#' @return A summary of the effect estimates can be obtained
#' with \code{summary} function.
#' @importFrom dplyr as_tibble mutate case_when select bind_cols
#' @import SuperLearner SuperLearner
#' @import tmle tmle
#' @references
#' Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021).
#' \emph{dplyr: A Grammar of Data Manipulation}.
#' R package version 1.0.7.
#' URL: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Eric Polley, Erin LeDell, Chris Kennedy and Mark van der Laan (2021).
#' \emph{SuperLearner: Super Learner Prediction}.
#' R package version 2.0-28.
#' URL:\url{https://CRAN.R-project.org/package=SuperLearner}
#'
#' Susan Gruber, Mark J. van der Laan (2012).
#' tmle: An R Package for Targeted Maximum Likelihood Estimation.
#'  \emph{Journal of Statistical Software}, \strong{51}(13), 1-35.
ce_estimate_tmle_ate <- function(y, w, x, sl_library, ...) {
  if (any((sl_library %in%
           getNamespaceExports("SuperLearner")[
             grepl(pattern = "^[S]L",
                   getNamespaceExports("SuperLearner"))]) == F))
    stop(
      "sl_library argument unrecgonized;
      please use listWrappers() in SuperLearner to
      find the list of supported values",
      call. = FALSE
    )
  n_trt <- length(unique(w))
  x_mat <- cbind(w, x)
  for (i in 1:n_trt) {
    x_mat <- x_mat %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(dplyr::case_when(w == i ~ 1,
                                     TRUE ~ 0)) %>%
      as.data.frame()
    names(x_mat)[length(x_mat)] <- paste0("w", i)
  }

  k <- n_trt

  n <- dim(x_mat)[1]
  t_mat <- NULL
  for (i in 1:n_trt) {
    t_mat_once <- x_mat %>%
      dplyr::select(paste0("w", i))
    t_mat <- dplyr::bind_cols(t_mat, t_mat_once)
  }

  w <- x
  #---------------------------------------------#
  ### Create Counterfactual Treatment Scenarios###
  #---------------------------------------------#
  for (i in 1:n_trt) {
    assign(paste0("w", i, "_countfactual"), w)
  }

  for (j in 1:n_trt) {
    for (i in 1:n_trt) {
      if (i == j) {
        names_w_countfactual <-
          names(eval(parse(text = (
            paste0("w", i, "_countfactual")
          ))))
        assign(paste0("w", i, "_countfactual"),
               eval(parse(text = paste0(
                 "w", i, "_countfactual"
               ))) %>% dplyr::mutate(1))
        assign(paste0("w", i, "_countfactual"),
               stats::setNames(eval(parse(
                 text = (paste0("w", i, "_countfactual"))
               )), c(
                 names_w_countfactual, paste0("w", j)
               )))
      } else {
        names_w_countfactual <-
          names(eval(parse(text = (
            paste0("w", i, "_countfactual")
          ))))
        assign(paste0("w", i, "_countfactual"),
               eval(parse(text = paste0(
                 "w", i, "_countfactual"
               ))) %>% dplyr::mutate(0))
        assign(paste0("w", i, "_countfactual"),
               stats::setNames(eval(parse(
                 text = (paste0("w", i, "_countfactual"))
               )), c(
                 names_w_countfactual, paste0("w", j)
               )))
      }
    }
  }

  w_countfactual_combined <- NULL
  for (i in 1:n_trt) {
    w_countfactual_combined <-
      as.data.frame(rbind(w_countfactual_combined, eval(parse(
        text = paste0("w", i, "_countfactual")
      ))))
  }



  # Run Super Learner Once, Obtain Initial Predicted Values
  # for All Counterfactual Settings###

  # Step 1: Estimating the outcome regression using super learner

  sl_fit <-
    SuperLearner::SuperLearner(
      Y = y,
      X = x_mat[, -1],
      newX = w_countfactual_combined,
      SL.library = sl_library,
      family = stats::binomial(),
      verbose = FALSE,
      ...
    )

  q_0 <- rep(0, n * k)
  q_tvector <- cbind(q_0, sl_fit$SL.predict)
  q_tmat <-
    matrix(unlist(split(
      as.data.frame(q_tvector), rep(1:k, each = n)
    )), ncol = 2 * k)


  # Run TMLE to Calculate Point Estimates of each T=t

  # Steps 2-5 are performed in this code chunk#

  w_results <- matrix(NA, nrow = k, ncol = 4)
  w_results_row_names <- NULL
  for (i in 1:n_trt) {
    w_results_row_names_once <- paste0("w", i)
    w_results_row_names <-
      c(w_results_row_names, w_results_row_names_once)
  }
  rownames(w_results) <- w_results_row_names

  colnames(w_results) <- c("EYt", "SE", "CI1", "CI2")
  start <- 1
  end <- 2
  for (t in 1:k) {
    # Step 2: Super learner fit for P(T_k=t|W) specified with
    # g.sl_library=sl_library in tmle call#
    # Steps 3-4: Performed in tmle call, target  EYt parameter
    # using A=NULL and Delta=Tmat[,t]#
    fit <- tmle::tmle(
      Y = y,
      A = NULL,
      Delta = t_mat[, t],
      W = w,
      Q = q_tmat[, c(start, end)],
      g.SL.library = sl_library,
      family = "binomial",
      verbose = FALSE,
      ...
    )
    # Step 5: The parameter estimates are stored in fit$estimates$EY1$psi#
    w_results[t, 1] <- fit$estimates$EY1$psi
    w_results[t, 2] <- fit$estimates$EY1$var.psi
    w_results[t, 3] <- fit$estimates$EY1$CI[1]
    w_results[t, 4] <- fit$estimates$EY1$CI[2]
    start <- start + 2
    end <- end + 2
  }

  w_result_one_repetition <- w_results %>%
    dplyr::as_tibble(rownames = "treatment")

  result <-
    list(
      method = "TMLE",
      result_TMLE = w_result_one_repetition,
      estimand = "ATE",
      n_trt = n_trt
    )
  class(result) <- "CIMTx_nonIPTW_once"
  return(result)
}
