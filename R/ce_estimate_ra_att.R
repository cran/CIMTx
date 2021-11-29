
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


  for (i in 1:n_trt) {
    assign(paste0("p",reference_trt,i,"_tilde"), arm::invlogit(as.matrix(sim_beta %*% t(eval(parse(text = paste0("x_tilde", reference_trt,i)))))))
    assign(paste0("y", reference_trt,i,"_tilde"), matrix(stats::rbinom(ndpost * eval(parse(text = paste0("n",reference_trt))), 1, eval(parse(text = paste0("p",reference_trt,i,"_tilde")))), nrow = ndpost))
    assign(paste0("y",i,"_pred"), rowMeans(eval(parse(text = paste0("y",reference_trt, i,"_tilde")))))
  }

  # Estimate causal effects
  for (j in 1:length(trt_indicator_no_reference)){
    assign(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est"), NULL)
    assign(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est"), NULL)
    assign(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"), NULL)

  }

  result <- NULL
  for (i in 1:(n_trt-1)){
    for (j in 1:length(trt_indicator_no_reference)){
      assign(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est"), list(eval(parse(text =(paste0("y",reference_trt,"_pred")))) - eval(parse(text =(paste0("y", trt_indicator_no_reference[j], "_pred"))))))
      assign(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est"), list(eval(parse(text =(paste0("y",reference_trt, "_pred")))) / eval(parse(text =(paste0("y", trt_indicator_no_reference[j], "_pred"))))))
      assign(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"), list((eval(parse(text =(paste0("y",reference_trt, "_pred")))) / (1 - eval(parse(text =(paste0("y",reference_trt, "_pred")))))) / (eval(parse(text =(paste0("y", trt_indicator_no_reference[j], "_pred")))) / (1 - eval(parse(text =(paste0("y", trt_indicator_no_reference[j], "_pred"))))))))
      assign(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est"), stats::setNames(eval(parse(text =(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est")))), paste0("ATT_RD",reference_trt, trt_indicator_no_reference[j])))
      assign(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est"), stats::setNames(eval(parse(text =(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est")))), paste0("ATT_RR",reference_trt, trt_indicator_no_reference[j])))
      assign(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"), stats::setNames(eval(parse(text =(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est")))), paste0("ATT_OR",reference_trt, trt_indicator_no_reference[j])))
      result <- c(result, (eval(parse(text =(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est"))))), (eval(parse(text =(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est"))))), (eval(parse(text =(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"))))))
    }
  }
  class(result) <- "CIMTx_ATT_posterior"
  return(result)
}
