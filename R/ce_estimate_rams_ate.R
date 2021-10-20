#' Regression adjustment with multivariate spline of GPS (RAMS) for ATE estimation
#'
#'This function implements the RAMS method when estimand is ATE. Please use our main function ce_estimate.R.
#'
#' @param y numeric vector for the binary outcome
#' @param w numeric vector for the treatment indicator
#' @param x dataframe including the treatment indicator and the covariates
#' @param method a character string. Users can selected from the following methods including "RAMS-Multinomial", "RAMS-GBM", "RAMS-SL"
#' @param ... Other paramters to be passed to twang::mnps()
#'
#' @return a list with w*(w-1)/2 elements for ATE effect. Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR
#'
#' @export
#'
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
#'ce_estimate_rams_ate(y = data$y, x = data$covariates ,
#'w = data$w, method = "RAMS-Multinomial")
ce_estimate_rams_ate = function(y, w, x, method,...) {

  n.trees <- parent.frame()$n.trees
  interaction.depth <- parent.frame()$interaction.depth

  n_trt <- length(unique(w))
  xwydata = as.data.frame(cbind(y = y, x, w = w))
  xwdata = as.data.frame(cbind( x, w = w))
  n <- nrow(xwydata)

  if (method == "RAMS-Multinomial") {
    psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
    pred_ps <- stats::fitted(psmod2)
    for (i in 1:n_trt){
      assign(paste0("ps",i), pred_ps[,i])
    }
  } else if (method == "RAMS-Multinomial-Trim") {
    trim_perc <- parent.frame()$trim_perc
      psmod2 <-  nnet::multinom(w~., data = xwdata,trace = FALSE)
      pred.ps <- stats::fitted(psmod2)
      for (i in 1:n_trt){
        assign(paste0("ps",i), trunc_fun(pred.ps[,i]))
      }
  } else if (method == "RAMS-GBM") {
    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod<-twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)), estimand = "ATE", treatATT = NULL,...)
    for (i in 1:n_trt){
      es.max.ATE <- NULL
      assign(paste0("ps", i), psmod$psList[[i]]$ps %>% pull(es.max.ATE))
    }
  } else if (method == "RAMS-GBM-Trim") {
    trim_perc <- parent.frame()$trim_perc
    temp<- noquote(names(x))
    strFormula  = sprintf("w~%s", paste(temp, sep = "",collapse="+"))
    psmod<-twang::mnps(stats::as.formula(strFormula),
                       data=xwdata %>% mutate(w = as.factor(w)), estimand = "ATE",
                       treatATT = NULL,...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), trunc_fun(psmod$psList[[i]]$ps %>% pull(es.max.ATE)))
    }

  } else if (method == "RAMS-SL") {
    SL.library <- parent.frame()$SL.library
    weightit_superlearner <- WeightIt::weightit(w~., data = xwdata,
                                                method = "super", estimand = "ATE",SL.library = SL.library,...)
    for (i in 1:n_trt){
      assign(paste0("ps", i), 1/weightit_superlearner$weights)
    }

  } else if (method == "RAMS-SL-Trim") {
    SL.library <- parent.frame()$SL.library
    trim_perc <- parent.frame()$trim_perc
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
