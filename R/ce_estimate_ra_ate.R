
ce_estimate_ra_ate <-  function(y, x, w, ndpost) {
  n_trt <- length(unique(w)) # Number of unique treatments
  for (i in 1:n_trt){
    assign(paste0("n",i), sum(w==i)) # Number of individuals receiving each treatment
  }
  n = length(w)

  xwdata = cbind(w,x)
  xwydata = cbind(y,xwdata)

  # Fit Bayesian logistic regression
  reg_mod = arm::bayesglm(y ~ ., data = as.data.frame(xwydata), family = stats::binomial(link="logit"), x = TRUE)

  mod_sims = arm::sim(reg_mod, n.sims = ndpost)
  sim_beta = as.matrix(stats::coef(mod_sims))
  x_tilde  = stats::model.matrix(reg_mod)

  for (i in 1:n_trt){
    assign(paste0("x_tilde",i), as.data.frame(x_tilde) %>%
             dplyr::mutate(w = i))
  }


  # predict potential outcomes

  for (i in 1:n_trt) {
      assign(paste0("p",i,"_tilde"), arm::invlogit(as.matrix(sim_beta %*% t(eval(parse(text = paste0("x_tilde",i)))))))
      assign(paste0("y",i,"_tilde"), matrix(stats::rbinom(ndpost * n, 1, eval(parse(text = paste0("p",i,"_tilde")))), nrow = ndpost))
  }


  # Estimate causal effects
  for (i in 1:(n_trt-1)){
    for (j in (i+1):(n_trt)){
      assign(paste0("RD",i,j, "_est"), NULL)
      assign(paste0("RR",i,j, "_est"), NULL)
      assign(paste0("OR",i,j, "_est"), NULL)
    }
  }

  for (i in 1:n_trt){
    assign(paste0("y",i, "_pred"), rowMeans(eval(parse(text = paste0("y",i,"_tilde")))))
  }
  result <- NULL
  for (i in 1:(n_trt-1)){
    for (j in (i + 1):n_trt){
      assign(paste0("RD",i,j, "_est"), list(eval(parse(text =(paste0("y",i, "_pred")))) - eval(parse(text =(paste0("y",j, "_pred"))))))
      assign(paste0("RR",i,j, "_est"), list(eval(parse(text =(paste0("y",i, "_pred")))) / eval(parse(text =(paste0("y",j, "_pred"))))))
      assign(paste0("OR",i,j, "_est"), list((eval(parse(text =(paste0("y",i, "_pred")))) / (1 - eval(parse(text =(paste0("y",i, "_pred")))))) / (eval(parse(text =(paste0("y",j, "_pred")))) / (1 - eval(parse(text =(paste0("y",j, "_pred"))))))))
      assign(paste0("RD",i,j, "_est"), stats::setNames(eval(parse(text =(paste0("RD",i,j, "_est")))), paste0("ATE_RD",i,j)))
      assign(paste0("RR",i,j, "_est"), stats::setNames(eval(parse(text =(paste0("RR",i,j, "_est")))), paste0("ATE_RR",i,j)))
      assign(paste0("OR",i,j, "_est"), stats::setNames(eval(parse(text =(paste0("OR",i,j, "_est")))), paste0("ATE_OR",i,j)))
      result <- c(result, (eval(parse(text =(paste0("RD",i,j, "_est"))))), (eval(parse(text =(paste0("RR",i,j, "_est"))))), (eval(parse(text =(paste0("OR",i,j, "_est"))))))
    }
  }

  class(result) <- "CIMTx_ATE_posterior"

  return(result)
}
