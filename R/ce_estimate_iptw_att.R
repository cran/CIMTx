#' Causal inference with multiple treatments using IPTW for ATT effects
#'
#' The function \code{ce_estimate_iptw_att} implements
#' IPTW to estimate ATT effect with
#' multiple treatments using observational data.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param w A numeric vector representing the treatment groups.
#' @param method A character string. Users can selected from the
#' following methods including \code{"IPTW-Multinomial"},
#' \code{"IPTW-GBM"}, \code{"IPTW-SL"}.
#' @param reference_trt A numeric value indicating reference treatment group
#' for ATT effect.
#' @param ... Other parameters that can be passed through to functions.
#'
#' @return A summary of the effect estimates can be obtained
#' with \code{summary} function. The weight distributions can be
#' visualized using \code{plot} function.
#' @importFrom nnet multinom
#' @importFrom WeightIt weightit
#' @references
#'
#' Venables, W. N. & Ripley, B. D. (2002)
#' \emph{Modern Applied Statistics with S}.
#' Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#'
#' Matthew Cefalu, Greg Ridgeway, Dan McCaffrey,
#' Andrew Morral, Beth Ann Griffin and Lane Burgette (2021).
#' \emph{twang: Toolkit for Weighting and Analysis of Nonequivalent Groups}.
#' R package version 2.5. URL:\url{https://CRAN.R-project.org/package=twang}
#'
#' Noah Greifer (2021).
#' \emph{WeightIt: Weighting for Covariate Balance in Observational Studies}.
#' R package version 0.12.0.
#' URL:\url{https://CRAN.R-project.org/package=WeightIt}
ce_estimate_iptw_att <-
  function(y, x, w, method, reference_trt, ...) {
    xwdata <- as.data.frame(cbind(x, w = w))
    n_trt <- length(unique(w))
    trt_indicator <- 1:n_trt
    trt_indicator_no_reference <-
      trt_indicator[trt_indicator != reference_trt]
    trim_perc <- parent.frame()$trim_perc
    if (method == "IPTW-SL" && is.null(trim_perc)) {
      sl_library <- parent.frame()$sl_library
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
      # Fit a SL model with treatment indicator as the outcome
      weightit_superlearner <-
        WeightIt::weightit(
          w ~ .,
          data = xwdata %>% mutate(w = as.factor(w)),
          focal = reference_trt,
          method = "super",
          estimand = "ATT",
          SL.library = sl_library,
          ...
        )
      # Extract the weights
      weight_superlearner <- weightit_superlearner$weights
      assign(paste0("mu_", reference_trt, "_hat_iptw_superlearner"),
             mean(y[w == reference_trt]))
      # Calculate the weighted means
      for (i in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0(
            "mu_",
            trt_indicator_no_reference[i],
            "_hat_iptw_superlearner"
          ),
          sum(y[w == trt_indicator_no_reference[i]] *
                weight_superlearner[w == trt_indicator_no_reference[i]]) /
            sum(weight_superlearner[w == trt_indicator_no_reference[i]])
        )
      }
      # Obtain the causal effects based on RD, OR and RR
      result_list_superlearner <- NULL
      for (j in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("RD", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_superlearner")
          ))) - eval(parse(text = (
            paste0(
              "mu_",
              trt_indicator_no_reference[j],
              "_hat_iptw_superlearner"
            )
          )))
        )
        assign(
          paste0("RR", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_superlearner")
          ))) / eval(parse(text = (
            paste0(
              "mu_",
              trt_indicator_no_reference[j],
              "_hat_iptw_superlearner"
            )
          )))
        )
        assign(
          paste0("OR", reference_trt, trt_indicator_no_reference[j]),
          (eval(parse(
            text = (paste0(
              "mu_", reference_trt, "_hat_iptw_superlearner"
            ))
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", reference_trt, "_hat_iptw_superlearner")
            ))
          ))) / (eval(parse(
            text = (
              paste0(
                "mu_",
                trt_indicator_no_reference[j],
                "_hat_iptw_superlearner"
              )
            )
          )) / (1 - eval(
            parse(text = (
              paste0(
                "mu_",
                trt_indicator_no_reference[j],
                "_hat_iptw_superlearner"
              )
            ))
          )))
        )
        result_once <-
          rbind(eval(parse(
            text = paste0("RD", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("RR", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("OR", reference_trt, trt_indicator_no_reference[j])
          )))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <-
          paste0("ATT", reference_trt, trt_indicator_no_reference[j])
        result_list_superlearner <-
          c(result_list_superlearner, result_once_list)
      }
      result_list_superlearner <-
        c(
          result_list_superlearner,
          list(weight = weight_superlearner),
          list(method = method),
          list(estimand = "ATT")
        )
      class(result_list_superlearner) <- "CIMTx_IPTW"
      return(result_list_superlearner)
    }
    if (method == "IPTW-Multinomial" && is.null(trim_perc)) {
      # Fit a multinomial logistic regression model with
      # treatment indicator as the outcome
      psmod2 <- nnet::multinom(w ~ ., data = xwdata, trace = FALSE)
      pred_ps <- stats::fitted(psmod2)
      # Get the weights
      for (j in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("att_wt_", reference_trt, trt_indicator_no_reference[j]),
          pred_ps[, reference_trt] / pred_ps[, trt_indicator_no_reference[j]]
        )
      }
      # Record the weights
      weight_glm <- NULL
      for (i in seq_len(length(trt_indicator_no_reference))) {
        weight_glm <-
          c(weight_glm, eval(parse(
            text = paste0("att_wt_", reference_trt,
                          trt_indicator_no_reference[i])
          ))[w == trt_indicator_no_reference[i]])
      }
      assign(paste0("mu_", reference_trt, "_hat_iptw"),
             mean(y[w == reference_trt]))
      # Calculate the weighted means
      for (i in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("mu_", trt_indicator_no_reference[i], "_hat_iptw"),
          sum(y[w == trt_indicator_no_reference[i]] * eval(parse(
            text = paste0("att_wt_", reference_trt,
                          trt_indicator_no_reference[i])
          ))[w == trt_indicator_no_reference[i]]) / sum(eval(parse(
            text = paste0("att_wt_", reference_trt,
                          trt_indicator_no_reference[i])
          ))[w == trt_indicator_no_reference[i]])
        )
      }
      # Obtain the causal effects based on RD, OR and RR
      result_list_multinomial <- NULL
      for (j in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("RD", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw")
          ))) - eval(parse(text = (
            paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw")
          )))
        )
        assign(
          paste0("RR", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw")
          ))) / eval(parse(text = (
            paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw")
          )))
        )
        assign(
          paste0("OR", reference_trt, trt_indicator_no_reference[j]),
          (eval(parse(
            text = (paste0("mu_", reference_trt, "_hat_iptw"))
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", reference_trt, "_hat_iptw")
            ))
          ))) / (eval(parse(
            text = (
              paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw")
            )
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw")
            ))
          )))
        )
        result_once <-
          rbind(eval(parse(
            text = paste0("RD", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("RR", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("OR", reference_trt, trt_indicator_no_reference[j])
          )))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <-
          paste0("ATT", reference_trt, trt_indicator_no_reference[j])
        result_list_multinomial <-
          c(result_list_multinomial, result_once_list)
      }
      result_list_multinomial <-
        c(
          result_list_multinomial,
          list(weight = weight_glm),
          list(method = method),
          list(estimand = "ATT")
        )
      class(result_list_multinomial) <- "CIMTx_IPTW"
      return(result_list_multinomial)
    }
    if (method == "IPTW-GBM" && is.null(trim_perc)) {
      # Fit a GBM model with treatment indicator as the outcome
      temp <- noquote(names(x))
      str_formula <-
        sprintf("w~%s", paste(temp, sep = "", collapse = "+"))
      psmod <- twang::mnps(
        stats::as.formula(str_formula),
        data = xwdata %>% mutate(w = as.factor(w)),
        estimand = "ATT",
        treatATT = reference_trt,
        ...
      )
      # Get the weights
      wt_hat <- twang::get.weights(psmod, estimand = "ATT")
      assign(paste0("mu_", reference_trt, "_hat_iptw_gbm"),
             mean(y[w == reference_trt]))
      # Calculate the weighted means
      for (i in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("mu_", trt_indicator_no_reference[i], "_hat_iptw_gbm"),
          sum(y[w == trt_indicator_no_reference[i]] *
                wt_hat[w == trt_indicator_no_reference[i]]) /
            sum(wt_hat[w == trt_indicator_no_reference[i]])
        )
      }
      # Obtain the causal effects based on RD, OR and RR
      result_list_gbm <- NULL
      for (j in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("RD", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_gbm")
          ))) - eval(parse(text = (
            paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw_gbm")
          )))
        )
        assign(
          paste0("RR", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_gbm")
          ))) / eval(parse(text = (
            paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw_gbm")
          )))
        )
        assign(
          paste0("OR", reference_trt, trt_indicator_no_reference[j]),
          (eval(parse(
            text = (paste0(
              "mu_", reference_trt, "_hat_iptw_gbm"
            ))
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", reference_trt, "_hat_iptw_gbm")
            ))
          ))) / (eval(parse(
            text = (
              paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw_gbm")
            )
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw_gbm")
            ))
          )))
        )
        result_once <-
          rbind(eval(parse(
            text = paste0("RD", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("RR", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("OR", reference_trt, trt_indicator_no_reference[j])
          )))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <-
          paste0("ATT", reference_trt, trt_indicator_no_reference[j])
        result_list_gbm <- c(result_list_gbm, result_once_list)
      }
      result_list_gbm <-
        c(
          result_list_gbm,
          list(weight = wt_hat),
          list(method = method),
          list(estimand = "ATT")
        )
      class(result_list_gbm) <- "CIMTx_IPTW"
      return(result_list_gbm)
    }
    if (method == "IPTW-Multinomial" && !is.null(trim_perc)) {
      # Fit a multinomial logistic regression model with
      # treatment indicator as the outcome
      psmod2 <- nnet::multinom(w ~ ., data = xwdata, trace = FALSE)
      pred_ps <- stats::fitted(psmod2)
      for (j in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("att_wt_", reference_trt, trt_indicator_no_reference[j]),
          pred_ps[, reference_trt] / pred_ps[, trt_indicator_no_reference[j]]
        )
      }
      # Get the trimmed weights
      for (i in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0(
            "att_wt_",
            reference_trt,
            trt_indicator_no_reference[i],
            "_trunc"
          ),
          trunc_fun(eval(parse(
            text = paste0("att_wt_", reference_trt,
                          trt_indicator_no_reference[i])
          )), trim_perc)
        )
      }
      weight_glm <- NULL
      for (i in seq_len(length(trt_indicator_no_reference))) {
        weight_glm <-
          c(weight_glm, eval(parse(
            text = paste0(
              "att_wt_",
              reference_trt,
              trt_indicator_no_reference[i],
              "_trunc"
            )
          ))[w == trt_indicator_no_reference[i]])
      }
      assign(paste0("mu_", reference_trt, "_hat_iptw_trim"),
             mean(y[w == reference_trt]))
      # Calculate the weighted means
      for (i in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("mu_", trt_indicator_no_reference[i], "_hat_iptw_trim"),
          sum(y[w == trt_indicator_no_reference[i]] * eval(parse(
            text = paste0(
              "att_wt_",
              reference_trt,
              trt_indicator_no_reference[i],
              "_trunc"
            )
          ))[w == trt_indicator_no_reference[i]]) / sum(eval(parse(
            text = paste0(
              "att_wt_",
              reference_trt,
              trt_indicator_no_reference[i],
              "_trunc"
            )
          ))[w == trt_indicator_no_reference[i]])
        )
      }
      # Obtain the causal effects based on RD, OR and RR
      result_list_multinomial_trim <- NULL
      for (j in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("RD", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_trim")
          ))) - eval(parse(text = (
            paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw_trim")
          )))
        )
        assign(
          paste0("RR", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_trim")
          ))) / eval(parse(text = (
            paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw_trim")
          )))
        )
        assign(
          paste0("OR", reference_trt, trt_indicator_no_reference[j]),
          (eval(parse(
            text = (paste0(
              "mu_", reference_trt, "_hat_iptw_trim"
            ))
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", reference_trt, "_hat_iptw_trim")
            ))
          ))) / (eval(parse(
            text = (
              paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw_trim")
            )
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", trt_indicator_no_reference[j], "_hat_iptw_trim")
            ))
          )))
        )
        result_once <-
          rbind(eval(parse(
            text = paste0("RD", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("RR", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("OR", reference_trt, trt_indicator_no_reference[j])
          )))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <-
          paste0("ATT", reference_trt, trt_indicator_no_reference[j])
        result_list_multinomial_trim <-
          c(result_list_multinomial_trim, result_once_list)
      }
      result_list_multinomial_trim <-
        c(
          result_list_multinomial_trim,
          list(weight = weight_glm),
          list(method = paste0(method, "-Trim")),
          list(estimand = "ATT")
        )
      class(result_list_multinomial_trim) <- "CIMTx_IPTW"
      return(result_list_multinomial_trim)
    }
    if (method == "IPTW-GBM" && !is.null(trim_perc)) {
      # Fit a GBM model with treatment indicator as the outcome
      temp <- noquote(names(x))
      str_formula <-
        sprintf("w~%s", paste(temp, sep = "", collapse = "+"))
      psmod <- twang::mnps(
        stats::as.formula(str_formula),
        data = xwdata %>% mutate(w = as.factor(w)),
        estimand = "ATT",
        treatATT = reference_trt,
        ...
      )
      wt_hat <- twang::get.weights(psmod, estimand = "ATT")
      # Get the trimmed weights
      for (i in seq_len(length(trt_indicator_no_reference))) {
        wt_hat[w == i] <- trunc_fun(wt_hat[w == i], trim_perc)
      }
      assign(paste0("mu_", reference_trt, "_hat_iptw_gbm_trim"),
             mean(y[w == reference_trt]))
      # Calculate the weighted means
      for (i in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0(
            "mu_",
            trt_indicator_no_reference[i],
            "_hat_iptw_gbm_trim"
          ),
          sum(y[w == trt_indicator_no_reference[i]] *
                wt_hat[w == trt_indicator_no_reference[i]]) /
            sum(wt_hat[w == trt_indicator_no_reference[i]])
        )
      }
      # Obtain the causal effects based on RD, OR and RR
      result_list_gbm_trim <- NULL
      for (j in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("RD", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_gbm_trim")
          ))) - eval(parse(text = (
            paste0(
              "mu_",
              trt_indicator_no_reference[j],
              "_hat_iptw_gbm_trim"
            )
          )))
        )
        assign(
          paste0("RR", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_gbm_trim")
          ))) / eval(parse(text = (
            paste0(
              "mu_",
              trt_indicator_no_reference[j],
              "_hat_iptw_gbm_trim"
            )
          )))
        )
        assign(
          paste0("OR", reference_trt, trt_indicator_no_reference[j]),
          (eval(parse(
            text = (paste0(
              "mu_", reference_trt, "_hat_iptw_gbm_trim"
            ))
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", reference_trt, "_hat_iptw_gbm_trim")
            ))
          ))) / (eval(parse(
            text = (
              paste0(
                "mu_",
                trt_indicator_no_reference[j],
                "_hat_iptw_gbm_trim"
              )
            )
          )) / (1 - eval(
            parse(text = (
              paste0(
                "mu_",
                trt_indicator_no_reference[j],
                "_hat_iptw_gbm_trim"
              )
            ))
          )))
        )
        result_once <-
          rbind(eval(parse(
            text = paste0("RD", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("RR", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("OR", reference_trt, trt_indicator_no_reference[j])
          )))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <-
          paste0("ATT", reference_trt, trt_indicator_no_reference[j])
        result_list_gbm_trim <-
          c(result_list_gbm_trim, result_once_list)
      }
      result_list_gbm_trim <-
        c(
          result_list_gbm_trim,
          list(weight = wt_hat),
          list(method = paste0(method, "-Trim")),
          list(estimand = "ATT")
        )
      class(result_list_gbm_trim) <- "CIMTx_IPTW"
      return(result_list_gbm_trim)
    }


    if (method == "IPTW-SL" && !is.null(trim_perc)) {
      # Fit a SL model with treatment indicator as the outcome
      sl_library <- parent.frame()$sl_library
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
      # Extract the weights
      weightit_superlearner <-
        WeightIt::weightit(
          w ~ .,
          data = xwdata %>% mutate(w = as.factor(w)),
          focal = reference_trt,
          method = "super",
          estimand = "ATT",
          SL.library = sl_library,
          ...
        )
      weight_superlearner <- weightit_superlearner$weights
      # Trim the weights
      weight_superlearner_trunc <-
        trunc_fun(weight_superlearner, trim_perc)
      for (i in seq_len(length(trt_indicator_no_reference))) {
        weight_superlearner_trunc[w == i] <-
          trunc_fun(weight_superlearner_trunc[w == i], trim_perc)
      }
      assign(paste0("mu_", reference_trt, "_hat_iptw_superlearner_trim"),
             mean(y[w == reference_trt]))
      # Calculate the weighted means
      for (i in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0(
            "mu_",
            trt_indicator_no_reference[i],
            "_hat_iptw_superlearner_trim"
          ),
          sum(y[w == trt_indicator_no_reference[i]] *
                weight_superlearner_trunc[w == trt_indicator_no_reference[i]]) /
            sum(weight_superlearner_trunc[w == trt_indicator_no_reference[i]])
        )
      }
      # Obtain the causal effects based on RD, OR and RR
      result_list_superlearner_trim <- NULL
      for (j in seq_len(length(trt_indicator_no_reference))) {
        assign(
          paste0("RD", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_superlearner_trim")
          ))) - eval(parse(text = (
            paste0(
              "mu_",
              trt_indicator_no_reference[j],
              "_hat_iptw_superlearner_trim"
            )
          )))
        )
        assign(
          paste0("RR", reference_trt, trt_indicator_no_reference[j]),
          eval(parse(text = (
            paste0("mu_", reference_trt, "_hat_iptw_superlearner_trim")
          ))) / eval(parse(text = (
            paste0(
              "mu_",
              trt_indicator_no_reference[j],
              "_hat_iptw_superlearner_trim"
            )
          )))
        )
        assign(
          paste0("OR", reference_trt, trt_indicator_no_reference[j]),
          (eval(parse(
            text = (
              paste0("mu_", reference_trt, "_hat_iptw_superlearner_trim")
            )
          )) / (1 - eval(
            parse(text = (
              paste0("mu_", reference_trt, "_hat_iptw_superlearner_trim")
            ))
          ))) / (eval(parse(
            text = (
              paste0(
                "mu_",
                trt_indicator_no_reference[j],
                "_hat_iptw_superlearner_trim"
              )
            )
          )) / (1 - eval(
            parse(text = (
              paste0(
                "mu_",
                trt_indicator_no_reference[j],
                "_hat_iptw_superlearner_trim"
              )
            ))
          )))
        )
        result_once <-
          rbind(eval(parse(
            text = paste0("RD", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("RR", reference_trt, trt_indicator_no_reference[j])
          )), eval(parse(
            text = paste0("OR", reference_trt, trt_indicator_no_reference[j])
          )))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <-
          paste0("ATT", reference_trt, trt_indicator_no_reference[j])
        result_list_superlearner_trim <-
          c(result_list_superlearner_trim, result_once_list)
      }
      result_list_superlearner_trim <-
        c(
          result_list_superlearner_trim,
          list(weight = weight_superlearner_trunc),
          list(method = paste0(method, "-Trim")),
          list(estimand = "ATT")
        )
      class(result_list_superlearner_trim) <- "CIMTx_nonIPTW_once"
      return(result_list_superlearner_trim)
    }
  }
