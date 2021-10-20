#' Regression Adjustment (RA) for ATT estimation
#'
#'This function implements the RA method when estimand is ATE. Please use our main function ce_estimate.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param w numeric vector for the treatment indicator
#' @param ndpost number of independent simulation draws to create
#' @param reference_trt Reference group for ATT
#'
#' @return a list with w-1 elements for ATT effect; Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR
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
#'ce_estimate_ra_att(y = data$y, x = data$covariates ,
#' w = data$w, reference_trt = 1, ndpost = 10)
ce_estimate_ra_att <- function(y, x, w, ndpost, reference_trt) {
  n_trt <- length(unique(w))
  for (i in 1:n_trt){
    assign(paste0("n",i), sum(w==i))
  }
  n = length(w)
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference_trt]

  xwdata = cbind(w,x)
  xwydata = cbind(y,xwdata)

  # Fit Bayesian logistic regression
  reg_mod = arm::bayesglm(y ~., data = xwydata, family = stats::binomial(link="logit"), x = TRUE)

  mod_sims = arm::sim(reg_mod, n.sims = ndpost)
  sim_beta = as.matrix(stats::coef(mod_sims))
  x_tilde  = as.data.frame(stats::model.matrix(reg_mod))

  for (i in 1:n_trt){
    assign(paste0("x_tilde", reference_trt, i), x_tilde[x_tilde[["w"]] == reference_trt, ])
  }

  for (i in 1:length(trt_indicator_no_reference)){
    assign(paste0("x_tilde", reference_trt, trt_indicator_no_reference[i]), eval(parse(text = paste0("x_tilde", reference_trt, trt_indicator_no_reference[i]))) %>%
             dplyr::mutate(w = trt_indicator_no_reference[i]))

  }

  for (i in 1:n_trt){
    assign(paste0("y",reference_trt, i,"_tilde"), NULL)
  }

  # predictive simulation using the stats::binomial distribution
  # predict potential outcomes

  for (s in 1:ndpost) {
    for (i in 1:n_trt) {
      assign(paste0("p",reference_trt, i,"_tilde"), arm::invlogit(as.matrix(eval(parse(text = paste0("x_tilde",reference_trt, i)))) %*% sim_beta[s,]))
      assign(paste0("y",reference_trt, i,"_tilde_", s),stats::rbinom(eval(parse(text = paste0("n",reference_trt))), 1, eval(parse(text = paste0("p",reference_trt, i,"_tilde")))))
      assign(paste0("y",reference_trt, i,"_tilde"),rbind(eval(parse(text = paste0("y",reference_trt, i,"_tilde"))),
                                                     eval(parse(text = paste0("y",reference_trt, i,"_tilde_", s)))))
    }
  }

  # Estimate causal effects
  for (j in 1:length(trt_indicator_no_reference)){
    assign(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est"), NULL)
    assign(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est"), NULL)
    assign(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"), NULL)
  }

  # RD12_est = RR12_est = OR12_est = NULL
  # RD13_est = RR13_est = OR13_est = NULL

  for (m in 1:ndpost) {
    # Estimate E(Y1), E(Y2), E(Y3)
    for (i in 1:n_trt){
      assign(paste0("y", i, "_pred_",m), mean(eval(parse(text = paste0("y",reference_trt,i,"_tilde"))) %>% as.data.frame %>% dplyr::slice(m) %>% as.numeric()))
    }

    for (j in 1:length(trt_indicator_no_reference)){
      assign(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est_m"), eval(parse(text =(paste0("y",reference_trt, "_pred_",m)))) - eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred_",m)))))
      assign(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est_m"), eval(parse(text =(paste0("y",reference_trt, "_pred_",m)))) / eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred_",m)))))
      assign(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est_m"), (eval(parse(text =(paste0("y",reference_trt, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",reference_trt, "_pred_",m)))))) / (eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred_",m)))))))
      assign(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est_m"))))))

      assign(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est_m"))))))

      assign(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est_m"))))))
    }
  }

  result <- NULL

  for (j in 1:length(trt_indicator_no_reference)){
    assign(paste0("att",reference_trt, trt_indicator_no_reference[j]), posterior_summary(eval(parse(text =(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"))))))
    assign(paste0("ATT",reference_trt, trt_indicator_no_reference[j]), list(round(eval(parse(text =(paste0("att",reference_trt, trt_indicator_no_reference[j])))), digits = 2)))
    assign(paste0("ATT",reference_trt, trt_indicator_no_reference[j]), stats::setNames(eval(parse(text =(paste0("ATT",reference_trt, trt_indicator_no_reference[j])))), paste0("ATT",reference_trt, trt_indicator_no_reference[j])))
    result <- c(result, (eval(parse(text =(paste0("ATT",reference_trt, trt_indicator_no_reference[j]))))))
  }
  return(result)
}
