# Bayesian Additive Regression Trees (BART) for ATE estimation
#'
#' This function implements the BART method when estimand is ATT. Please use our main function ce_estimate.R.
#'
#' @param y a numeric vector (0, 1) representing a binary outcome
#' @param x a dataframe, including all the covariates but not treatments.
#' @param w a numeric vector representing the treatment groups
#' @param discard "No" or "Yes" indicating whether to use the discarding rules for the BART based method. The default is "No"
#' @param reference_trt reference treatment group for ATT effect
#' @param ndpost number of posterior samples from BART
#' @param ... Other parameters that can be passed through the BART::pbart() function
#'
#' @return a list with w-1 elements for ATT effect; Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR
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
#'ce_estimate_bart_att(y = data$y, x = data$covariates , w = data$w, ndpost = 10,reference_trt=1)
#'
ce_estimate_bart_att <- function(y, x, w, discard = "No", ndpost = 1000,reference_trt,...) {

  n_trt <- length(unique(w))
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference_trt]

  for (i in 1:n_trt){
    assign(paste0("n",i), sum(w==i))
  }

  xwdata = cbind(w,x)

  # Fit BART
  bart_mod = BART::pbart(x.train = xwdata, y.train = y, ndpost = ndpost, ...)

  # Predict potential outcomes for among those in reference_trt group

  assign(paste0("xp",reference_trt), xwdata[w==reference_trt,])
  for (i in trt_indicator[trt_indicator!=reference_trt]){
    assign(paste0("xp",i), xwdata[w==reference_trt,])
    assign(paste0("xp",i), as.data.frame(eval(parse(text = paste0("xp",i)))) %>%
             dplyr::mutate(w = i))
  }

  for (j in 1:(n_trt)){
      assign(paste0("bart_pred",reference_trt,j), BART::pwbart(eval(parse(text = paste0("xp",j))), bart_mod$treedraws, mu=mean(y)))
      assign(paste0("pred_prop",reference_trt,j), stats::pnorm(eval(parse(text = paste0("bart_pred",reference_trt,j)))))
    }

  if (discard == "No") {

    for (j in 1:(n_trt)){
      assign(paste0("pred_prop",reference_trt,j), stats::pnorm(eval(parse(text = paste0("bart_pred",reference_trt,j)))))
    }
    n_discard_att <- 0
  } else if (discard == "Yes"){
    for (j in 1:(n_trt)){
      assign(paste0("post.ind.sd",j), apply(eval(parse(text = paste0("pred_prop",reference_trt,j))), 2, stats::sd))
    }
    threshold <- max(eval(parse(text = paste0("post.ind.sd",reference_trt))))

    for (j in 1:length((trt_indicator_no_reference))) {
      assign(paste0("eligible",trt_indicator_no_reference[j]), (eval(parse(text = paste0("post.ind.sd",trt_indicator_no_reference[j]))) <= threshold))
    }
    eligible <- rep(TRUE, dim(eval(parse(text = paste0("pred_prop",reference_trt,1))))[2])
    for (j in 1:length((trt_indicator_no_reference))){
      eligible <- (eval(parse(text = paste0("eligible",trt_indicator_no_reference[j]))) & eligible)

    }
    n_discard_att <- sum(eligible == FALSE)
    for (j in 1:(n_trt)){
      assign(paste0("pred_prop",reference_trt,j), (eval(parse(text = paste0("pred_prop",reference_trt,j))) %>%
               as.data.frame() %>%
               dplyr::select(which(eligible)) %>%
               as.matrix()))
    }

  }
  # Estimate causal effects
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j], "_est"), NULL)
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j], "_est"), NULL)
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j], "_est"), NULL)
    }

  for (m in 1:ndpost) {
      for (j in 1:n_trt){
        assign(paste0("y",j,"_pred"), mean(stats::rbinom(eval(parse(text =(paste0("n",reference_trt)))), 1, eval(parse(text =(paste0("pred_prop",reference_trt,j)))))))
      }
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j], "_est_m"), eval(parse(text =(paste0("y",reference_trt, "_pred")))) - eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred")))))
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j], "_est_m"), eval(parse(text =(paste0("y",reference_trt, "_pred")))) / eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred")))))
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j], "_est_m"), (eval(parse(text =(paste0("y",reference_trt, "_pred")))) / (1 - eval(parse(text =(paste0("y",reference_trt, "_pred")))))) / (eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred")))) / (1 - eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred")))))))
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("RD",reference_trt,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RD",reference_trt,trt_indicator_no_reference[j], "_est_m"))))))

      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("RR",reference_trt,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference_trt,trt_indicator_no_reference[j], "_est_m"))))))

      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("OR",reference_trt,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference_trt,trt_indicator_no_reference[j], "_est_m"))))))
    }

  }
  result <- NULL
  for (j in 1:length((trt_indicator_no_reference))){
    assign(paste0("att",reference_trt,trt_indicator_no_reference[j]), posterior_summary(eval(parse(text =(paste0("RD",reference_trt,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference_trt,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference_trt,trt_indicator_no_reference[j], "_est"))))))
    assign(paste0("ATT",reference_trt,trt_indicator_no_reference[j]), list(round(eval(parse(text =(paste0("att",reference_trt,trt_indicator_no_reference[j])))), digits = 2)))
    assign(paste0("ATT",reference_trt,trt_indicator_no_reference[j]), stats::setNames(eval(parse(text =(paste0("ATT",reference_trt,trt_indicator_no_reference[j])))), paste0("ATT",reference_trt,trt_indicator_no_reference[j])))
    result <- c(result, (eval(parse(text =(paste0("ATT",reference_trt,trt_indicator_no_reference[j]))))))
  }

  return(c(result, list(n_discard = n_discard_att)))
}
