#' Causal inference with multiple treatments using RAMS for ATT effects
#'
#' The function \code{ce_estimate_rams_att} implements
#' RAMS to estimate ATT effect with
#' multiple treatments using observational data.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param w A numeric vector representing the treatment groups.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param method A character string. Users can selected from the
#' following methods including \code{"RAMS-Multinomial"},
#' \code{"RAMS-GBM"}, \code{"RAMS-SL"}.
#' @param reference_trt A numeric value indicating
#' reference treatment group for ATT effect.
#' @param ... Other parameters that can be passed through to functions.
#'
#' @return A summary of the effect estimates can be obtained
#' with \code{summary} function.
#' @importFrom mgcv gam
#' @importFrom nnet multinom
#' @importFrom WeightIt weightit
#' @references
#' Matthew Cefalu, Greg Ridgeway, Dan McCaffrey, Andrew Morral,
#'  Beth Ann Griffin and Lane Burgette (2021).
#'  \emph{twang: Toolkit for Weighting and Analysis of Nonequivalent Groups}.
#'  R package version 2.5.
#'  URL:\url{https://CRAN.R-project.org/package=twang}
#'
#' Venables, W. N. & Ripley, B. D. (2002)
#' \emph{Modern Applied Statistics with S}.
#' Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#'
#' Noah Greifer (2021).
#' \emph{WeightIt: Weighting for Covariate Balance in Observational Studies}.
#' R package version 0.12.0.
#' URL:\url{https://CRAN.R-project.org/package=WeightIt}
#'
#' Wood, S.N. (2011)
#' Fast stable restricted maximum likelihood and
#' marginal likelihood estimation of semiparametric generalized linear models.
#' \emph{Journal of the Royal Statistical Society (B)} \strong{73}(1):3-36
ce_estimate_rams_att <-
  function(y, w, x, method, reference_trt, ...) {
    n.trees <- parent.frame()$n.trees
    interaction.depth <- parent.frame()$interaction.depth
    n_trt <- length(unique(w))
    xwydata <- as.data.frame(cbind(y = y, x, w = w))
    xwdata <- as.data.frame(cbind(x, w = w))
    trt_indicator <- 1:n_trt
    trt_indicator_no_reference <-
      trt_indicator[trt_indicator != reference_trt]
    n_trt <- length(unique(w))
    trim_perc <- parent.frame()$trim_perc
    for (i in 1:n_trt) {
      assign(paste0("n", i), sum(w == i))
    }

    if (method == "RAMS-Multinomial" && is.null(trim_perc)) {
      # Fit a multinomial logistic regression model with
      # treatment indicator as the outcome
      psmod2 <- nnet::multinom(w ~ ., data = xwdata, trace = FALSE)
      pred_ps <- stats::fitted(psmod2)
      for (i in 1:n_trt) {
        assign(paste0("ps", i), pred_ps[, i])
      }
    } else if (method == "RAMS-Multinomial-Trim" &&
               !is.null(trim_perc)) {
      # Fit a multinomial logistic regression model with
      # treatment indicator as the outcome
      psmod2 <- nnet::multinom(w ~ ., data = xwdata, trace = FALSE)
      pred_ps <- stats::fitted(psmod2)
      # Trim the PS
      for (i in 1:n_trt) {
        assign(paste0("ps", i), trunc_fun(pred_ps[, i]))
      }
    } else if (method == "RAMS-GBM" && is.null(trim_perc)) {
      es_max_ate <- NULL
      # Fit a GBM model with treatment indicator as the outcome
      temp <- noquote(names(x))
      str_formula <-
        sprintf("w~%s", paste(temp, sep = "", collapse = "+"))
      psmod <- twang::mnps(stats::as.formula(str_formula),
                           data = xwdata %>% mutate(w = as.factor(w)),
                           ...)
      for (i in 1:n_trt) {
        assign(paste0("ps", i), psmod$psList[[i]]$ps %>% pull(es_max_ate))
      }
    } else if (method == "RAMS-GBM-Trim" && !is.null(trim_perc)) {
      temp <- noquote(names(x))
      str_formula <-
        sprintf("w~%s", paste(temp, sep = "", collapse = "+"))
      psmod <- twang::mnps(stats::as.formula(str_formula),
                           data = xwdata %>% mutate(w = as.factor(w)),
                           ...)
      # Trim the PS
      for (i in 1:n_trt) {
        assign(paste0("ps", i), trunc_fun(psmod$psList[[i]]$ps %>%
                                            pull(es_max_ate)))
      }
    } else if (method == "RAMS-SL" && is.null(trim_perc)) {
      # Fit a SL model with treatment indicator as the outcome
      sl_library <- parent.frame()$sl_library
      weightit_superlearner <- WeightIt::weightit(as.factor(w) ~ .,
                                                  data = xwdata,
                                                  method = "super",
                                                  SL.library = sl_library,
                                                  ...)
      for (i in 1:n_trt) {
        assign(paste0("ps", i), 1 / weightit_superlearner$weights)
      }
    } else if (method == "RAMS-SL-Trim" && !is.null(trim_perc)) {
      # Fit a SL model with treatment indicator as the outcome
      sl_library <- parent.frame()$sl_library
      weightit_superlearner <- WeightIt::weightit(as.factor(w) ~ .,
                                                  data = xwdata,
                                                  method = "super",
                                                  SL.library = sl_library,
                                                  ...)
      # Trim the PS
      for (i in 1:n_trt) {
        assign(paste0("ps", i),
               trunc_fun(1 / weightit_superlearner$weights))
      }
    }
    # logit of propensity scores
    logit_ps1 <- NULL
    logit_ps2 <- NULL
    for (i in 1:n_trt) {
      assign(paste0("logit_ps", i), stats::qlogis(eval(parse(text = paste0(
        "ps", i
      )))))
    }

    mod_splinedat <-
      as.data.frame(cbind(
        w = xwydata$w,
        logit_ps1 = logit_ps1,
        logit_ps32 = logit_ps2
      ))
    # Fit a generalized additive model using the treatment indicator and
    # multivariate spline function of the logit of GPS as the predictors
    mod_spline <-
      mgcv::gam(
        y ~ w + te(logit_ps1, logit_ps2),
        family = stats::binomial(link = "logit"),
        data = mod_splinedat
      )
    # Predict the potential outcomes using the fitted GAM model
    for (i in 1:n_trt) {
      assign(
        paste0("newdata", i),
        data.frame(
          w = rep(i, sum(w == reference_trt)),
          logit_ps1 = logit_ps1[w == reference_trt],
          logit_ps2 = logit_ps2[w == reference_trt]
        )
      )
      assign(paste0("spline.pred", i),
             stats::plogis(stats::predict(mod_spline,
                                          newdata = eval(
        parse(text = paste0("newdata", i))
      ))))
      assign(paste0("y", i, ".hat"), mean(eval(parse(
        text = paste0("spline.pred", i)
      ))))
    }
    # Estimate the causal effects in terms of OR, RR and RD
    result_list_rams_att <- NULL
    for (j in seq_len(length(trt_indicator_no_reference))) {
      assign(
        paste0("RD", reference_trt, trt_indicator_no_reference[j]),
        eval(parse(text = (
          paste0("y", reference_trt, ".hat")
        ))) - eval(parse(text = (
          paste0("y", trt_indicator_no_reference[j], ".hat")
        )))
      )
      assign(
        paste0("RR", reference_trt, trt_indicator_no_reference[j]),
        eval(parse(text = (
          paste0("y", reference_trt, ".hat")
        ))) / eval(parse(text = (
          paste0("y", trt_indicator_no_reference[j], ".hat")
        )))
      )
      assign(paste0("OR", reference_trt, trt_indicator_no_reference[j]),
             (eval(parse(text = (
               paste0("y", reference_trt, ".hat")
             ))) / (1 - eval(parse(
               text = (paste0("y", reference_trt, ".hat"))
             )))) / (eval(parse(text = (
               paste0("y", trt_indicator_no_reference[j], ".hat")
             ))) / (1 - eval(parse(
               text = (paste0(
                 "y", trt_indicator_no_reference[j], ".hat"
               ))
             )))))
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
      result_list_rams_att <-
        c(result_list_rams_att, result_once_list)
    }
    if (!is.null(trim_perc)) {
      method <- paste0(method, "-Trim")
    }
    result_list_rams_att <-
      c(result_list_rams_att, list(estimand = "ATT"), method = method)
    class(result_list_rams_att) <- "CIMTx_nonIPTW_once"
    return(result_list_rams_att)
  }
