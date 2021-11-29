

ce_estimate_iptw_att <- function (y, x, w, method, reference_trt,...){

  xwdata = as.data.frame(cbind(x, w = w))
  n_trt <- length(unique(w))
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference_trt]
  trim_perc <- parent.frame()$trim_perc
  if (method == "IPTW-SL" && is.null(trim_perc)) {
    SL.library <- parent.frame()$SL.library
    if (any((SL.library %in% getNamespaceExports("SuperLearner")[grepl(pattern = "^[S]L", getNamespaceExports("SuperLearner"))]) == F)) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values", call. = FALSE)
    weightit_superlearner <- WeightIt::weightit(w~., data = xwdata %>% mutate(w = as.factor(w)),focal = reference_trt, method = "super", estimand = "ATT",SL.library = SL.library,...)
    weight_superlearner <- weightit_superlearner$weights
    assign(paste0("mu_",reference_trt,"_hat_iptw_superlearner"), mean(y[w == reference_trt]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_superlearner"), sum(y[w == trt_indicator_no_reference[i]] * weight_superlearner[w == trt_indicator_no_reference[i]]) / sum(weight_superlearner[w == trt_indicator_no_reference[i]]))
    }
    result_list_superlearner <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_superlearner")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner")))))
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_superlearner")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner")))))
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_superlearner")))) / (1 - eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_superlearner")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference_trt,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference_trt,trt_indicator_no_reference[j])
      result_list_superlearner <- c(result_list_superlearner, result_once_list)
    }
    return(c(result_list_superlearner, list(weight = weight_superlearner), list(method = method)))
  }

  #1a_ to compute ATTs using LR estimated weights
  if (method == "IPTW-Multinomial" && is.null(trim_perc)) {
    # use logistic Multinomial model with main effects only to estimate ps
    psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
    pred_ps <- stats::fitted(psmod2)
    for (j in 1:length(trt_indicator_no_reference)){
      assign(paste0("att_wt_",reference_trt, trt_indicator_no_reference[j]), pred_ps[,reference_trt]/pred_ps[,trt_indicator_no_reference[j]])
    }
    weight_glm <- NULL
    for (i in 1:length(trt_indicator_no_reference)){
      weight_glm <- c(weight_glm, eval(parse(text = paste0("att_wt_",reference_trt, trt_indicator_no_reference[i])))[w == trt_indicator_no_reference[i]])
    }
    # att_wt_12 <- pred_ps[,1]/pred_ps[,2]
    # att_wt_13 <- pred_ps[,1]/pred_ps[,3]
    assign(paste0("mu_",reference_trt,"_hat_iptw"), mean(y[w == reference_trt]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw"), sum(y[w == trt_indicator_no_reference[i]] * eval(parse(text = paste0("att_wt_",reference_trt, trt_indicator_no_reference[i])))[w == trt_indicator_no_reference[i]]) / sum(eval(parse(text = paste0("att_wt_",reference_trt,trt_indicator_no_reference[i])))[w == trt_indicator_no_reference[i]]))
    }
    result_list_multinomial <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw")))))
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw")))))
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw")))) / (1 - eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference_trt,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference_trt,trt_indicator_no_reference[j])
      result_list_multinomial <- c(result_list_multinomial, result_once_list)
    }
    return(c(result_list_multinomial, list(weight = weight_glm), list(method = method)))
  }


  #1b_ to compute ATTs using GBM estimated weights
  if (method == "IPTW-GBM" && is.null(trim_perc)) {
    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod<-twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)), estimand = "ATT",
                       treatATT = reference_trt,...)
    wt_hat<- twang::get.weights(psmod,estimand = "ATT")
    assign(paste0("mu_",reference_trt,"_hat_iptw_gbm"), mean(y[w == reference_trt]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_gbm"), sum(y[w == trt_indicator_no_reference[i]] * wt_hat[w == trt_indicator_no_reference[i]]) / sum(wt_hat[w == trt_indicator_no_reference[i]]))
    }
    result_list_gbm <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_gbm")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm")))))
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_gbm")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm")))))
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_gbm")))) / (1 - eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_gbm")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference_trt,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference_trt,trt_indicator_no_reference[j])
      result_list_gbm <- c(result_list_gbm, result_once_list)
    }
    return(c(result_list_gbm, list(weight = wt_hat), list(method = method)))
  }

  if (method == "IPTW-Multinomial" && !is.null(trim_perc)) {
    # trim_perc <- parent.frame()$trim_perc
    psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
    pred_ps <- stats::fitted(psmod2)
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("att_wt_",reference_trt, trt_indicator_no_reference[j]), pred_ps[,reference_trt]/pred_ps[,trt_indicator_no_reference[j]])
    }

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("att_wt_",reference_trt,trt_indicator_no_reference[i],"_trunc"), trunc_fun(eval(parse(text = paste0("att_wt_",reference_trt, trt_indicator_no_reference[i]))), trim_perc))
    }

    weight_glm <- NULL
    for (i in 1:length(trt_indicator_no_reference)){
      weight_glm <- c(weight_glm, eval(parse(text = paste0("att_wt_",reference_trt, trt_indicator_no_reference[i],"_trunc")))[w == trt_indicator_no_reference[i]])
    }
    assign(paste0("mu_",reference_trt,"_hat_iptw_trim"), mean(y[w == reference_trt]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_trim"), sum(y[w == trt_indicator_no_reference[i]] * eval(parse(text = paste0("att_wt_",reference_trt, trt_indicator_no_reference[i],"_trunc")))[w == trt_indicator_no_reference[i]]) / sum(eval(parse(text = paste0("att_wt_",reference_trt,trt_indicator_no_reference[i],"_trunc")))[w == trt_indicator_no_reference[i]]))
    }
    result_list_multinomial_trim <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_trim")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_trim")))))
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_trim")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_trim")))))
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_trim")))) / (1 - eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_trim")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_trim")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_trim")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference_trt,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference_trt,trt_indicator_no_reference[j])
      result_list_multinomial_trim <- c(result_list_multinomial_trim, result_once_list)
    }
    return(c(result_list_multinomial_trim, list(weight = weight_glm), list(method = paste0(method, "-Trim"))))
  }


  if (method == "IPTW-GBM"&& !is.null(trim_perc)) {
    # trim_perc <- parent.frame()$trim_perc
    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod<-twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)), estimand = "ATT", treatATT = reference_trt,...)
    wt_hat<- twang::get.weights(psmod, estimand = "ATT")

    for (i in 1:length((trt_indicator_no_reference))){
      wt_hat[w==i] <- trunc_fun(wt_hat[w==i], trim_perc)
    }
    assign(paste0("mu_",reference_trt,"_hat_iptw_gbm_trim"), mean(y[w == reference_trt]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_gbm_trim"), sum(y[w == trt_indicator_no_reference[i]] * wt_hat[w == trt_indicator_no_reference[i]]) / sum(wt_hat[w == trt_indicator_no_reference[i]]))
    }
    result_list_gbm_trim <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_gbm_trim")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm_trim")))))
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_gbm_trim")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm_trim")))))
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_gbm_trim")))) / (1 - eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_gbm_trim")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm_trim")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm_trim")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference_trt,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference_trt,trt_indicator_no_reference[j])
      result_list_gbm_trim <- c(result_list_gbm_trim, result_once_list)
    }
    return(c(result_list_gbm_trim, list(weight = wt_hat), list(method = paste0(method, "-Trim"))))
  }


  if (method == "IPTW-SL"&& !is.null(trim_perc)) {
    # trim_perc <- parent.frame()$trim_perc
    SL.library <- parent.frame()$SL.library
    if (any((SL.library %in% getNamespaceExports("SuperLearner")[grepl(pattern = "^[S]L", getNamespaceExports("SuperLearner"))]) == F)) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values", call. = FALSE)
    weightit_superlearner <- WeightIt::weightit(w~., data = xwdata %>% mutate(w = as.factor(w)),focal = reference_trt, method = "super", estimand = "ATT",SL.library = SL.library,...)
    weight_superlearner <- weightit_superlearner$weights
    weight_superlearner_trunc <- trunc_fun(weight_superlearner, trim_perc)
    for (i in 1:length((trt_indicator_no_reference))){
      weight_superlearner_trunc[w==i] <- trunc_fun(weight_superlearner_trunc[w==i], trim_perc)
    }
    assign(paste0("mu_",reference_trt,"_hat_iptw_superlearner_trim"), mean(y[w == reference_trt]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_superlearner_trim"), sum(y[w == trt_indicator_no_reference[i]] * weight_superlearner_trunc[w == trt_indicator_no_reference[i]]) / sum(weight_superlearner_trunc[w == trt_indicator_no_reference[i]]))
    }
    result_list_superlearner_trim <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_superlearner_trim")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner_trim")))))
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_superlearner_trim")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner_trim")))))
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_superlearner_trim")))) / (1 - eval(parse(text =(paste0("mu_",reference_trt, "_hat_iptw_superlearner_trim")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner_trim")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner_trim")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference_trt,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference_trt,trt_indicator_no_reference[j])
      result_list_superlearner_trim <- c(result_list_superlearner_trim, result_once_list)
    }
    return(c(result_list_superlearner_trim, list(weight = weight_superlearner_trunc), list(method = paste0(method, "-Trim"))))
  }
}

