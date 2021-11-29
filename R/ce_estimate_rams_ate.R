
ce_estimate_rams_ate = function(y, w, x, method,...) {

  n.trees <- parent.frame()$n.trees
  interaction.depth <- parent.frame()$interaction.depth

  n_trt <- length(unique(w))
  xwydata = as.data.frame(cbind(y = y, x, w = w))
  xwdata = as.data.frame(cbind( x, w = w))
  n <- nrow(xwydata)
  trim_perc <- parent.frame()$trim_perc
  if (method == "RAMS-Multinomial" && is.null(trim_perc)) {
    psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
    pred_ps <- stats::fitted(psmod2)
    for (i in 1:n_trt){
      assign(paste0("ps",i), pred_ps[,i])
    }
  } else if (method == "RAMS-Multinomial" && !is.null(trim_perc)) {

      psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
      pred.ps <- stats::fitted(psmod2)
      for (i in 1:n_trt){
        assign(paste0("ps",i), trunc_fun(pred.ps[,i]))
      }
  } else if (method == "RAMS-GBM" && is.null(trim_perc)) {
    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod<-twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)), estimand = "ATE", treatATT = NULL,...)
    for (i in 1:n_trt){
      es.max.ATE <- NULL
      assign(paste0("ps", i), psmod$psList[[i]]$ps %>% pull(es.max.ATE))
    }
  } else if (method == "RAMS-GBM"&& !is.null(trim_perc)) {
    # trim_perc <- parent.frame()$trim_perc
    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod<-twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)), estimand = "ATE",
                       treatATT = NULL,...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), trunc_fun(psmod$psList[[i]]$ps %>% pull(es.max.ATE)))
    }

  } else if (method == "RAMS-SL" && is.null(trim_perc)) {
    SL.library <- parent.frame()$SL.library
    if (any((SL.library %in% getNamespaceExports("SuperLearner")[grepl(pattern = "^[S]L", getNamespaceExports("SuperLearner"))]) == F)) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values", call. = FALSE)
    weightit_superlearner <- WeightIt::weightit(w~., data = xwdata,
                                                method = "super", estimand = "ATE",SL.library = SL.library,...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), 1/weightit_superlearner$weights)
    }

  } else if (method == "RAMS-SL"&& !is.null(trim_perc)) {
    SL.library <- parent.frame()$SL.library
    # trim_perc <- parent.frame()$trim_perc
    if (any((SL.library %in% getNamespaceExports("SuperLearner")[grepl(pattern = "^[S]L", getNamespaceExports("SuperLearner"))]) == F)) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values", call. = FALSE)
    weightit_superlearner <- WeightIt::weightit(w~., data = xwdata,
                                                method = "super", estimand = "ATE",SL.library = SL.library,...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), trunc_fun(1/weightit_superlearner$weights, trim_perc = trim_perc))
    }
  }
  # logit of propensity scores
  logit_ps1 <- NULL
  logit_ps2 <- NULL

  for (i in 1:n_trt){
    assign(paste0("logit_ps", i), stats::qlogis(eval(parse(text = paste0("ps", i)))))
  }

  mod.splinedat = as.data.frame(cbind(w = xwydata$w, logit_ps1, logit_ps2))
  mod.spline = mgcv::gam(y ~ w + te(logit_ps1,logit_ps2), family = stats::binomial(link="logit"), data = mod.splinedat)

  for (i in 1:n_trt){
    assign(paste0("newdata",i), data.frame(w = rep(i,n), logit_ps1, logit_ps2))
    assign(paste0("spline.pred",i), stats::plogis(stats::predict(mod.spline, newdata = eval(parse(text = paste0("newdata",i))))))
    assign(paste0("y",i,".hat"), mean(eval(parse(text = paste0("spline.pred",i)))))
  }

  result_list <- NULL
  for (i in 1:(n_trt-1)){
    result_once <- NULL
    for (j in (i + 1):n_trt){
      assign(paste0("RD",i,j),eval(parse(text = paste0("y",i, ".hat"))) - eval(parse(text = paste0("y",j, ".hat"))))
      assign(paste0("RR",i,j),eval(parse(text = paste0("y",i, ".hat"))) / eval(parse(text = paste0("y",j, ".hat"))))
      assign(paste0("OR",i,j), (eval(parse(text = paste0("y",i, ".hat"))) /(1 - eval(parse(text = paste0("y",i, ".hat"))))) / (eval(parse(text = paste0("y",j, ".hat"))) /(1 - eval(parse(text = paste0("y",j, ".hat"))))))
      result_once <- round(rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j)))),2)
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATE",i,j)
      result_list <- c(result_list, result_once_list)
    }
  }
  return(result_list)
}
