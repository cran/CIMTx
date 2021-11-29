
ce_estimate_iptw_ate <- function (y, w, x, method,...){

  xwdata = as.data.frame(cbind(x, w = w))
  n_trt <- length(unique(w))
  trim_perc <- parent.frame()$trim_perc
  if (method == "IPTW-Multinomial" && is.null(trim_perc)) {
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
    class(result_list_multinomial) <- "CIMTx_IPTW"
    return(result_list_multinomial)
  } else if (method == "IPTW-Multinomial" && !is.null(trim_perc)) {

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
    result_list_multinomial_trim <- c(result_list_multinomial_trim, list(weight = weight_glm_trim), list(method = paste0(method, "-Trim")))
    class(result_list_multinomial_trim) <- "CIMTx_IPTW"
    return(result_list_multinomial_trim)

  } else if (method == "IPTW-GBM" && is.null(trim_perc)) {
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
    class(result_list_gbm) <- "CIMTx_IPTW"
    return(result_list_gbm)
  } else if (method == "IPTW-GBM"&& !is.null(trim_perc)) {
    # n.trees <- parent.frame()$n.trees
    # interaction.depth <- parent.frame()$interaction.depth
    es.max.ATE <- NULL
    # trim_perc <- parent.frame()$trim_perc
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
    result_list_gbm_trim <- c(result_list_gbm_trim, list(weight = wt_hat_trunc), list(method = paste0(method, "-Trim")))
    class(result_list_gbm_trim) <- "CIMTx_IPTW"
    return(result_list_gbm_trim)
  } else if (method == "IPTW-SL"&& is.null(trim_perc)) {
    SL.library <- parent.frame()$SL.library
    if (any((SL.library %in% getNamespaceExports("SuperLearner")[grepl(pattern = "^[S]L", getNamespaceExports("SuperLearner"))]) == F)) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values", call. = FALSE)
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
    class(result_list_superlearner) <- "CIMTx_IPTW"
    return(result_list_superlearner)
  } else if (method == "IPTW-SL"&& !is.null(trim_perc)) {
    SL.library <- parent.frame()$SL.library
    if (any((SL.library %in% getNamespaceExports("SuperLearner")[grepl(pattern = "^[S]L", getNamespaceExports("SuperLearner"))]) == F)) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values", call. = FALSE)
    # trim_perc <- parent.frame()$trim_perc
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
    result_list_superlearner_trim <- c(result_list_superlearner_trim, list(weight = weight_superlearner_trim), list(method = paste0(method, "-Trim")))
    class(result_list_superlearner_trim) <- "CIMTx_IPTW"
    return(result_list_superlearner_trim)

  }
}
