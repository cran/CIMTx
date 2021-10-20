#' Bayesian Additive Regression Trees (BART) for ATE estimation
#'
#' This function implements the BART method when estimand is ATE. Please use our main function ce_estimate.R.
#'
#' @param y a numeric vector (0, 1) representing a binary outcome
#' @param x a dataframe, including all the covariates but not treatments
#' @param w a numeric vector representing the treatment groups
#' @param ndpost number of posterior samples from BART
#' @param discard "No" or "Yes" indicating whether to use the discarding rules for the BART based method. The default is "No"
#' @param ... Other parameters that can be passed through the BART::pbart() function
#'
#' @return a list with w*(w-1)/2 elements for ATE effect. Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR
#' @export
#' @examples
#' library(CIMTx)
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
#'ce_estimate_bart_ate(y = data$y, x = data$covariates , w = data$w, ndpost = 10)
#'
ce_estimate_bart_ate <- function(y, x, w, discard = "No",ndpost = 1000,... ) {

  n_trt <- length(unique(w))
  for (i in 1:n_trt){
      assign(paste0("n",i), sum(w==i))
    }
    xwdata = cbind(w,x)
    # Fit BART
    bart_mod = BART::pbart(x.train = xwdata, y.train = y, ndpost = ndpost,...)

    # Predict potential outcomes for w
    for (i in 1:n_trt){
      assign(paste0("xp",i), xwdata)
      for (j in 1:(n_trt)){
        assign(paste0("xp",j),  as.data.frame(eval(parse(text =paste0("xp",i)))) %>%
                 dplyr::mutate(w = j))
        assign(paste0("bart_pred",i,j), stats::predict(bart_mod, newdata = eval(parse(text = paste0("xp",j)))))
        assign(paste0("pred_prop",i,j), eval(parse(text = paste0("bart_pred",i,j)))[["prob.test"]])
      }
    }
    if (discard == "No") {
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("pred_prop",i,j), eval(parse(text = paste0("bart_pred",i,j)))[["prob.test"]])
        }
      }
      discard_all <- 0
    } else if(discard == "Yes"){
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("post.ind.sd",i,j), apply(eval(parse(text = paste0("pred_prop",i,j))), 2, stats::sd))
        }
        assign(paste0("threshold",i), max(eval(parse(text = paste0("post.ind.sd",i,i)))))
      }
   # Discarding rule
      for (i in 1:n_trt){
        n_trt_no_i <- unique(w)[unique(w)!=i]
        assign(paste0("eligible",i), TRUE)
        assign(paste0("criteria",i), TRUE)
        for (j in 1:length((n_trt_no_i))) {
          assign(paste0("criteria",i,n_trt_no_i[j]), eval(parse(text = paste0("post.ind.sd",i,n_trt_no_i[j]))) <= eval(parse(text = paste0("threshold",i))))
          assign(paste0("eligible",i), eval(parse(text = paste0("criteria",i,n_trt_no_i[j]))) & eval(parse(text = paste0("criteria",i))))
        }

        assign(paste0("n_",i, "_discard"), sum(eval(parse(text = paste0("eligible",i)))== FALSE))
      }
      discard_all <- NULL
      for (i in 1:n_trt){
        discard_all <- c(discard_all, eval(parse(text = paste0("n_",i, "_discard"))))
      }
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("pred_prop",i,j), eval(parse(text = paste0("pred_prop",i,j))) %>%
                   as.data.frame() %>%
                   dplyr::select(which(eval(parse(text = paste0("eligible",i))))) %>%
                   as.matrix())
        }
      }
    }

    for (i in 1:(n_trt-1)){
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j, "_est"), NULL)
        assign(paste0("RR",i,j, "_est"), NULL)
        assign(paste0("OR",i,j, "_est"), NULL)
      }
    }

    for (m in 1:ndpost) {
      for (i in 1:n_trt){
        assign(paste0("y",i), NULL)
        for (j in 1:n_trt){
          assign(paste0("y",j,i), c(stats::rbinom(eval(parse(text =(paste0("n",j)))), 1, eval(parse(text =(paste0("pred_prop",j,i)))) %>% as.data.frame %>% dplyr::slice(m) %>% as.numeric())))
          assign(paste0("y",i), c(eval(parse(text =(paste0("y",i)))), eval(parse(text =(paste0("y",j,i))))))
        }
        assign(paste0("y",i, "_pred_",m), mean(eval(parse(text =(paste0("y",i))))))
      }

      for (i in 1:(n_trt-1)){
        for (j in (i + 1):n_trt){
          assign(paste0("RD",i,j, "_est_m"), eval(parse(text =(paste0("y",i, "_pred_",m)))) - eval(parse(text =(paste0("y",j, "_pred_",m)))))
          assign(paste0("RR",i,j, "_est_m"), eval(parse(text =(paste0("y",i, "_pred_",m)))) / eval(parse(text =(paste0("y",j, "_pred_",m)))))
          assign(paste0("OR",i,j, "_est_m"), (eval(parse(text =(paste0("y",i, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",i, "_pred_",m)))))) / (eval(parse(text =(paste0("y",j, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",j, "_pred_",m)))))))
          assign(paste0("RD",i,j, "_est"), c(eval(parse(text =(paste0("RD",i,j, "_est")))), eval(parse(text =(paste0("RD",i,j, "_est_m"))))))

          assign(paste0("RR",i,j, "_est"), c(eval(parse(text =(paste0("RR",i,j, "_est")))), eval(parse(text =(paste0("RR",i,j, "_est_m"))))))

          assign(paste0("OR",i,j, "_est"), c(eval(parse(text =(paste0("OR",i,j, "_est")))), eval(parse(text =(paste0("OR",i,j, "_est_m"))))))
        }
      }
    }

result <- NULL
    for (i in 1:(n_trt-1)){
      for (j in (i + 1):n_trt){
        assign(paste0("ate",i,j), posterior_summary(eval(parse(text =(paste0("RD",i,j, "_est")))), eval(parse(text =(paste0("RR",i,j, "_est")))), eval(parse(text =(paste0("OR",i,j, "_est"))))))
        assign(paste0("ATE",i,j), list(round(eval(parse(text =(paste0("ate",i,j)))), digits = 2)))
        assign(paste0("ATE",i,j), stats::setNames(eval(parse(text =(paste0("ATE",i,j)))), paste0("ATE",i,j)))
        result <- c(result, (eval(parse(text =(paste0("ATE",i,j))))))
      }
    }

return(c(result, list(n_discard = discard_all)))
}
