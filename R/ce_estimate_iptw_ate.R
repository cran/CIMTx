#' Inverse probability of treatment weighting (IPTW) for ATE estimation
#'
#'This function implements the IPTW method when estimand is ATE. Please use our main function ce_estimate.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x a dataframe, including all the covariates but not treatments
#' @param w numeric vector for the treatment indicator
#' @param method a character string. Users can selected from the following methods including "IPTW-Multinomial", "IPTW-GBM", "IPTW-SL"
#' @param ... Other parameters that can be passed through the twang::GBM() function
#'
#' @return a list with w*(w-1)/2 elements for ATE effect. Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR
#' @export
#' @examples
#'lp_w_all <-
#'  c(".4*x1 + .1*x2  - .1*x4 + .1*x5",    # w = 1
#'    ".2 * x1 + .2 * x2  - .2 * x4 - .3 * x5")  # w = 2
#'nlp_w_all <-
#'  c("-.5*x1*x4  - .1*x2*x5", # w = 1
#'    "-.3*x1*x4 + .2*x2*x5")# w = 2
#'lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - .1*x4 - .2*x5", 3)
#'nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
#'X_all <- c(
#'  "rnorm(300, 0, 0.5)",# x1
#'  "rbeta(300, 2, .4)",   # x2
#'  "runif(300, 0, 0.5)",# x3
#'  "rweibull(300,1,2)",  # x4
#'  "rbinom(300, 1, .4)"# x5
#')

#'set.seed(111111)
#'data <- data_sim(
#'  sample_size = 300,
#'  n_trt = 3,
#'  X = X_all,
#'  lp_y = lp_y_all,
#'  nlp_y  = nlp_y_all,
#'  align = FALSE,
#'  lp_w = lp_w_all,
#'  nlp_w = nlp_w_all,
#'  tau = c(-1.5,0,1.5),
#'  delta = c(0.5,0.5),
#'  psi = 1
#')
#'ce_estimate_iptw_ate(y = data$y, x = data$covariates,
#'w = data$w, ndpost=100, method = "IPTW-Multinomial")
ce_estimate_iptw_ate <- function (y, w, x, method,...){

  xwdata = as.data.frame(cbind(x, w = w))
  n_trt <- length(unique(w))

  if (method == "IPTW-Multinomial") {
    psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
    pred_ps <- stats::fitted(psmod2)
    for (i in 1:n_trt){
      assign(paste0("ate_wt_",i), 1/pred_ps[,i])
    }
    weight_glm <- NULL
    for (i in 1:n_trt){
      weight_glm <- c(weight_glm, eval(parse(text = paste0("ate_wt_",i)))[w == i])
    }
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hat_iptw"), sum(y[w == i] * weight_glm[w == i]) / sum(weight_glm[w == i]))
    }
    result_list_multinomial <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat_iptw"))) - eval(parse(text = paste0("mu_",j, "_hat_iptw"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat_iptw"))) / eval(parse(text = paste0("mu_",j, "_hat_iptw"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat_iptw"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat_iptw"))))) / (eval(parse(text = paste0("mu_",j, "_hat_iptw"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat_iptw"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_multinomial <- c(result_list_multinomial, result_once_list)
      }
    }
    result_list_multinomial <- c(result_list_multinomial, list(weight = weight_glm), list(method = method))
    return(result_list_multinomial)
  } else if (method == "IPTW-Multinomial-Trim") {
    trim_perc <- parent.frame()$trim_perc
    psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
    pred_ps <- stats::fitted(psmod2)
    for (i in 1:n_trt){
      assign(paste0("ate_wt_",i), 1/pred_ps[,i])
    }
    weight_glm <- NULL
    for (i in 1:n_trt){
      weight_glm <- c(weight_glm, eval(parse(text = paste0("ate_wt_",i)))[w == i])
    }
    weight_glm_trim <- trunc_fun(weight_glm, trim_perc)
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hat_iptw_trim"), sum(y[w == i] * weight_glm_trim[w == i]) / sum(weight_glm_trim[w == i]))
    }
    result_list_multinomial_trim <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat_iptw_trim"))) - eval(parse(text = paste0("mu_",j, "_hat_iptw_trim"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat_iptw_trim"))) / eval(parse(text = paste0("mu_",j, "_hat_iptw_trim"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat_iptw_trim"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat_iptw_trim"))))) / (eval(parse(text = paste0("mu_",j, "_hat_iptw_trim"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat_iptw_trim"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_multinomial_trim <- c(result_list_multinomial_trim, result_once_list)
      }
    }
    result_list_multinomial_trim <- c(result_list_multinomial_trim, list(weight = weight_glm_trim), list(method = method))
    return(result_list_multinomial_trim)

  } else if (method == "IPTW-GBM") {
    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod <- twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)), estimand = "ATE",)
    for (i in 1:n_trt){
      assign(paste0("ps", i), psmod$psList[[i]]$ps %>% pull(es.max.ATE))
    }
    wt_hat<- twang::get.weights(psmod, estimand = "ATE")

    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hatgbm"), sum(y[w == i] * wt_hat[w == i]) / sum(wt_hat[w == i]))
    }
    result_list_gbm <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hatgbm"))) - eval(parse(text = paste0("mu_",j, "_hatgbm"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hatgbm"))) / eval(parse(text = paste0("mu_",j, "_hatgbm"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hatgbm"))) /(1 - eval(parse(text = paste0("mu_",i, "_hatgbm"))))) / (eval(parse(text = paste0("mu_",j, "_hatgbm"))) /(1 - eval(parse(text = paste0("mu_",j, "_hatgbm"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_gbm <- c(result_list_gbm, result_once_list)
      }
    }
    result_list_gbm <- c(result_list_gbm, list(weight = wt_hat), list(method = method))
    return(result_list_gbm)
  } else if (method == "IPTW-GBM-Trim") {
    # n.trees <- parent.frame()$n.trees
    # interaction.depth <- parent.frame()$interaction.depth
    es.max.ATE <- NULL
    trim_perc <- parent.frame()$trim_perc
    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod<-twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)), estimand = "ATE",...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), psmod$psList[[i]]$ps %>% pull(es.max.ATE))
    }
    wt_hat<- twang::get.weights(psmod,estimand = "ATE")

    wt_hat_trunc <- trunc_fun(wt_hat, trim_perc)
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hatgbm_trim"), sum(y[w == i] * wt_hat_trunc[w == i]) / sum(wt_hat_trunc[w == i]))
    }
    result_list_gbm_trim <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hatgbm_trim"))) - eval(parse(text = paste0("mu_",j, "_hatgbm_trim"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hatgbm_trim"))) / eval(parse(text = paste0("mu_",j, "_hatgbm_trim"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hatgbm_trim"))) /(1 - eval(parse(text = paste0("mu_",i, "_hatgbm_trim"))))) / (eval(parse(text = paste0("mu_",j, "_hatgbm_trim"))) /(1 - eval(parse(text = paste0("mu_",j, "_hatgbm_trim"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_gbm_trim <- c(result_list_gbm_trim, result_once_list)
      }
    }
    result_list_gbm_trim <- c(result_list_gbm_trim, list(weight = wt_hat_trunc), list(method = method))
    return(result_list_gbm_trim)
  } else if (method == "IPTW-SL") {
    SL.library <- parent.frame()$SL.library
      weightit_superlearner <- WeightIt::weightit(w~., data = xwdata,
                                                  method = "super", estimand = "ATE",SL.library = SL.library,...)
    weight_superlearner <- weightit_superlearner$weights

    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hat_superlearner"), sum(y[w == i] * weight_superlearner[w == i]) / sum(weight_superlearner[w == i]))
    }
    result_list_superlearner <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat_superlearner"))) - eval(parse(text = paste0("mu_",j, "_hat_superlearner"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat_superlearner"))) / eval(parse(text = paste0("mu_",j, "_hat_superlearner"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat_superlearner"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat_superlearner"))))) / (eval(parse(text = paste0("mu_",j, "_hat_superlearner"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat_superlearner"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_superlearner <- c(result_list_superlearner, result_once_list)
      }
    }
    result_list_superlearner <- c(result_list_superlearner, list(weight = weight_superlearner), list(method = method))
    return(result_list_superlearner)
  } else if (method == "IPTW-SL-Trim") {
    SL.library <- parent.frame()$SL.library
    trim_perc <- parent.frame()$trim_perc
      weightit_superlearner <- WeightIt::weightit(w~., data = xwdata,
                                                  method = "super", estimand = "ATE",
                                                  SL.library = SL.library,...)

      weight_superlearner_trim <- trunc_fun(weightit_superlearner$weights,trim_perc)
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hat_superlearner_trim"), sum(y[w == i] * weight_superlearner_trim[w == i]) / sum(weight_superlearner_trim[w == i]))
    }
    result_list_superlearner_trim <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat_superlearner_trim"))) - eval(parse(text = paste0("mu_",j, "_hat_superlearner_trim"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat_superlearner_trim"))) / eval(parse(text = paste0("mu_",j, "_hat_superlearner_trim"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat_superlearner_trim"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat_superlearner_trim"))))) / (eval(parse(text = paste0("mu_",j, "_hat_superlearner_trim"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat_superlearner_trim"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_superlearner_trim <- c(result_list_superlearner_trim, result_once_list)
      }
    }
    result_list_superlearner_trim <- c(result_list_superlearner_trim, list(weight = weight_superlearner_trim), list(method = method))
    return(result_list_superlearner_trim)

  }
}
