
ce_estimate_rams_att <- function(y, w, x, method,reference_trt,...) {
  n_trt <- length(unique(w))
  xwydata = as.data.frame(cbind(y = y, x, w = w))
  xwdata = as.data.frame(cbind(x, w = w))
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference_trt]
  n_trt <- length(unique(w))
  trim_perc <- parent.frame()$trim_perc
  for (i in 1:n_trt){
    assign(paste0("n",i), sum(w==i))
  }
  n <- nrow(xwydata)

  if (method == "RAMS-Multinomial" && is.null(trim_perc)) {
    psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
    pred_ps <- stats::fitted(psmod2)
    for (i in 1:n_trt){
      assign(paste0("ps",i), pred_ps[,i])
    }
  } else if (method == "RAMS-Multinomial-Trim" && !is.null(trim_perc)) {
    # trim_alpha <- parent.frame()$trim_alpha
    psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
    pred.ps <- stats::fitted(psmod2)
    for (i in 1:n_trt){
      assign(paste0("ps",i), trunc_fun(pred.ps[,i]))
    }
  } else if (method == "RAMS-GBM" && is.null(trim_perc)) {
    es.max.ATE <- NULL
    n.trees <- parent.frame()$n.trees
    interaction.depth <- parent.frame()$interaction.depth

    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod <- twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)),...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), psmod$psList[[i]]$ps %>% pull(es.max.ATE))
    }
  } else if (method == "RAMS-GBM-Trim" && !is.null(trim_perc)) {
    # trim_alpha <- parent.frame()$trim_alpha
    n.trees <- parent.frame()$n.trees
    interaction.depth <- parent.frame()$interaction.depth

    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod<-twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)),...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), trunc_fun(psmod$psList[[i]]$ps %>% pull(es.max.ATE)))
    }

  } else if (method == "RAMS-SL" && is.null(trim_perc)) {
    SL.library <- parent.frame()$SL.library
    weightit_superlearner <- WeightIt::weightit(w~., data = xwdata,
                                                method = "super",SL.library = SL.library,...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), 1/weightit_superlearner$weights)
    }

  } else if (method == "RAMS-SL-Trim" && !is.null(trim_perc)) {
    SL.library <- parent.frame()$SL.library
    weightit_superlearner <- WeightIt::weightit(w~., data = xwdata,
                                                method = "super",SL.library = SL.library,...
    )
    for (i in 1:n_trt){
      assign(paste0("ps", i), trunc_fun(1/weightit_superlearner$weights))
    }
  }
  # logit of propensity scores
  # logit_ps1 = stats::qlogis(ps1)
  #
  # logit_ps2 = stats::qlogis(ps2)
  logit_ps1 <- NULL
  logit_ps2 <- NULL

  for (i in 1:n_trt){
    assign(paste0("logit_ps", i), stats::qlogis(eval(parse(text = paste0("ps", i)))))
  }

  mod.splinedat = as.data.frame(cbind(w = xwydata$w, logit_ps1 = logit_ps1, logit_ps32= logit_ps2))
  mod.spline = mgcv::gam(y ~ w + te(logit_ps1,logit_ps2), family = stats::binomial(link="logit"), data = mod.splinedat)

  # # predict potential outcomes Y(1)
  for (i in 1:n_trt){
    assign(paste0("newdata",i), data.frame(w = rep(i,sum(w == reference_trt)), logit_ps1=logit_ps1[w == reference_trt], logit_ps2=logit_ps2[w == reference_trt]))
    assign(paste0("spline.pred",i), stats::plogis(stats::predict(mod.spline, newdata = eval(parse(text = paste0("newdata",i))))))
    assign(paste0("y",i,".hat"), mean(eval(parse(text = paste0("spline.pred",i)))))
  }

  result_list_rams_att <- NULL
  for (j in 1:length((trt_indicator_no_reference))){
    assign(paste0("RD",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("y",reference_trt, ".hat")))) - eval(parse(text =(paste0("y",trt_indicator_no_reference[j], ".hat")))))
    assign(paste0("RR",reference_trt,trt_indicator_no_reference[j]), eval(parse(text =(paste0("y",reference_trt, ".hat")))) / eval(parse(text =(paste0("y",trt_indicator_no_reference[j], ".hat")))))
    assign(paste0("OR",reference_trt,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("y",reference_trt, ".hat")))) / (1 - eval(parse(text =(paste0("y",reference_trt, ".hat")))))) / (eval(parse(text =(paste0("y",trt_indicator_no_reference[j], ".hat")))) / (1 - eval(parse(text =(paste0("y",trt_indicator_no_reference[j], ".hat")))))))
    result_once <- rbind(eval(parse(text = paste0("RD",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference_trt,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference_trt,trt_indicator_no_reference[j]))))
    colnames(result_once) <- "EST"
    rownames(result_once) <- c("RD", "RR", "OR")
    result_once_list <- list(result_once)
    names(result_once_list) <- paste0("ATT",reference_trt,trt_indicator_no_reference[j])
    result_list_rams_att <- c(result_list_rams_att, result_once_list)
  }
  return(result_list_rams_att)
}
