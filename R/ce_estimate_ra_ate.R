#' Regression Adjustment (RA) for ATE estimation
#'
#'This function implements the RA method when estimand is ATE. Please use our main function ce_estimate.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param w numeric vector for the treatment indicator
#' @param ndpost number of independent simulation draws to create
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
#'ce_estimate_ra_ate(y = data$y, x = data$covariates
#', w = data$w, ndpost = 10)
#'
ce_estimate_ra_ate <-  function(y, x, w, ndpost) {
  n_trt <- length(unique(w)) # Number of unique treatments
  for (i in 1:n_trt){
    assign(paste0("n",i), sum(w==i)) # Number of individuals receiving each treatment
  }
  n = length(w)

  xwdata = cbind(w,x)
  xwydata = cbind(y,xwdata)

  # Fit Bayesian logistic regression
  reg_mod = arm::bayesglm(y ~ ., data = xwydata, family = stats::binomial(link="logit"), x = TRUE)

  mod_sims = arm::sim(reg_mod, n.sims = ndpost)
  sim_beta = as.matrix(stats::coef(mod_sims))
  x_tilde  = stats::model.matrix(reg_mod)

  for (i in 1:n_trt){
    assign(paste0("x_tilde",i), as.data.frame(x_tilde) %>%
             dplyr::mutate(w = i))
  }

  # predictive simulation using the stats::binomial distribution
  # predict potential outcomes
  for (i in 1:n_trt){
    assign(paste0("y",i,"_tilde"), NULL)
  }

  for (s in 1:ndpost) {
    for (i in 1:n_trt) {
      assign(paste0("p",i,"_tilde"), arm::invlogit(as.matrix(eval(parse(text = paste0("x_tilde",i)))) %*% sim_beta[s,]))
      assign(paste0("y",i,"_tilde_", s),stats::rbinom(n, 1, eval(parse(text = paste0("p",i,"_tilde")))))
      assign(paste0("y",i,"_tilde"),rbind(eval(parse(text = paste0("y",i,"_tilde"))),
                                          eval(parse(text = paste0("y",i,"_tilde_", s)))))
    }
  }

  # Estimate causal effects
  for (i in 1:(n_trt-1)){
    for (j in (i+1):(n_trt)){
      assign(paste0("RD",i,j, "_est"), NULL)
      assign(paste0("RR",i,j, "_est"), NULL)
      assign(paste0("OR",i,j, "_est"), NULL)
    }
  }

  for (m in 1:ndpost) {
    # Estimate E(Y1), E(Y2), E(Y3)
    for (i in 1:n_trt){
      assign(paste0("y",i, "_pred_",m), mean(eval(parse(text = paste0("y",i,"_tilde"))) %>% as.data.frame %>% dplyr::slice(m) %>% as.numeric()))
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
  return(result)
}
