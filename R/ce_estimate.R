#' Causal inference with multiple treatments using observational data
#'
#' This function implements the 6 different methods for causal inference with multiple treatments using observational data.
#'
#' @param y a numeric vector (0, 1) representing a binary outcome
#' @param x a dataframe, including all the covariates but not treatments.
#' @param w a numeric vector representing the treatment groups
#' @param method a character string. Users can selected from the following methods including "RA", "VM", "BART", "TMLE", "IPTW-Multinomial", "IPTW-GBM", "IPTW-SL", "RAMS-Multinomial", "RAMS-GBM", "RAMS-SL"
#' @param discard "No" or "Yes" indicating whether to use the discarding rules for the BART based method. The default is "No"
#' @param estimand "ATT" or "ATE" representing the type of causal estimand. When the estimand = "ATT", users also need to specify the reference treatment group by setting the reference_trt argument.
#' @param trim_perc the percentile at which the inverse probability of treatment weights shouldbe trimmed
#' @param SL.library a character vector of prediction algorithms. A list of functions included in the SuperLearner package can be found with listWrappers().
#' @param reference_trt reference treatment group for ATT effect
#' @param boot is logical, indicating whether or not to use nonparametric bootstrap to calculate the 95\% confidence intervals of the causal effect estimates.
#' @param nboots only need to set up when the boot = TRUE; is a numeric value representing the number of bootstrap samples
#' @param ndpost is the number of posterior draws for the Bayesian methods
#' @param caliper only need to set up when the method is set to VM; is a numeric value denoting the caliper which should be used when matching on the logit of GPS within each cluster formed by K-means clustering. The caliper is in standardized units. For example, caliper = 0.25 means that all matches greater than 0.25 standard deviations of the logit of GPS are dropped. The default value is 0.25
#' @param n_cluster only need to set up when the method is set to VM; a numeric value denoting the number of clusters to form using K means clustering on the logit of GPS. The default value is 5.
#' @param ... Other parameters that can be passed through the functions
#'
#' @return a list with w-1 elements for ATT effect; a list with w*(w-1)/2 elements for ATE effect. Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR
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
#'ce_estimate(y = data$y, x = data$covariates, w = data$w,
#'ndpost=100, method = "RA", estimand = "ATE")
ce_estimate <- function(y, x, w, method, discard = "No", estimand, trim_perc, SL.library, reference_trt,  boot = FALSE,nboots, ndpost = 1000,caliper = 0.25,n_cluster  = 5,...){
  if (method == "RA" && estimand == "ATE") {
    result <- ce_estimate_ra_ate(
      y = y,
      x = x,
      w = w,
      ndpost = ndpost
    )
  } else if (method == "RA" && estimand == "ATT"){
    result <- ce_estimate_ra_att(
      y = y,
      x = x,
      w = w,
      ndpost = ndpost,
      reference_trt = reference_trt
    )
  } else if (method == "VM" && estimand == "ATT"&& boot == FALSE){
    result <- ce_estimate_vm_att(
      y = y,
      x = x,
      w = w,
      reference_trt = reference_trt,
      caliper = caliper,
      n_cluster = n_cluster
    )
  } else if (method == "VM" && estimand == "ATT"&& boot == TRUE){
    result <- ce_estimate_vm_att_boot(
      y = y,
      x = x,
      w = w,
      reference_trt = reference_trt,
      caliper = caliper,
      n_cluster = n_cluster,
      nboots = nboots
    )
  } else if (method == "BART" && estimand == "ATE") {
    result <- ce_estimate_bart_ate(
      y = y,
      x = x,
      w = w,
      ndpost = ndpost,
      discard = discard,...
    )
  } else if (method == "BART" && estimand == "ATT") {
    result <- ce_estimate_bart_att(
      y = y,
      x = x,
      w = w,
      ndpost = ndpost,
      reference_trt = reference_trt,
      discard = discard,...
    )
  } else if (method == "TMLE" && estimand == "ATE") {
    result <- ce_estimate_tmle_ate(
      y = y,
      x = x,
      w = w,
      SL.library = SL.library,...
    )
  } else if (method %in% c("RAMS-Multinomial", "RAMS-Multinomial-Trim", "RAMS-GBM", "RAMS-GBM-Trim", "RAMS-SL", "RAMS-SL-Trim") && estimand == "ATE"&& boot == FALSE) {
    result <- ce_estimate_rams_ate(
      y = y,
      x = x,
      w = w,
      method = method,...
    )
  } else if (method %in% c("RAMS-Multinomial", "RAMS-Multinomial-Trim", "RAMS-GBM", "RAMS-GBM-Trim", "RAMS-SL", "RAMS-SL-Trim") && estimand == "ATE"&& boot == TRUE) {
    result <- ce_estimate_rams_ate_boot(
      y = y,
      x = x,
      w = w,
      method = method,
      nboots = nboots,...
    )
  } else if (method %in% c("RAMS-Multinomial", "RAMS-Multinomial-Trim", "RAMS-GBM", "RAMS-GBM-Trim", "RAMS-SL", "RAMS-SL-Trim") && estimand == "ATT" && boot == FALSE) {
    result <- ce_estimate_rams_att(
      y = y,
      x = x,
      w = w,
      method = method,
      reference_trt = reference_trt,...
    )
  } else if (method %in% c("RAMS-Multinomial", "RAMS-Multinomial-Trim", "RAMS-GBM", "RAMS-GBM-Trim", "RAMS-SL", "RAMS-SL-Trim") && estimand == "ATT" && boot == TRUE) {
    result <- ce_estimate_rams_att_boot(
      y = y,
      x = x,
      w = w,
      method = method,
      nboots = nboots,
      reference_trt = reference_trt,...
    )
  } else if (method %in% c("IPTW-Multinomial", "IPTW-Multinomial-Trim", "IPTW-GBM", "IPTW-GBM-Trim", "IPTW-SL", "IPTW-SL-Trim") && estimand == "ATE"&& boot == FALSE) {
    result <- ce_estimate_iptw_ate(
      y = y,
      x = x,
      w = w,
      method = method,...
    )
  }  else if (method %in% c("IPTW-Multinomial", "IPTW-Multinomial-Trim", "IPTW-GBM", "IPTW-GBM-Trim", "IPTW-SL", "IPTW-SL-Trim") && estimand == "ATE"&& boot == TRUE) {
    result <- ce_estimate_iptw_ate_boot(
      y = y,
      x = x,
      w = w,
      method = method,
      nboots = nboots,...
    )
  } else if (method %in% c("IPTW-Multinomial", "IPTW-Multinomial-Trim", "IPTW-GBM", "IPTW-GBM-Trim", "IPTW-SL", "IPTW-SL-Trim") && estimand == "ATT"&& boot == FALSE) {
    result <- ce_estimate_iptw_att(
      y = y,
      x = x,
      w = w,
      method = method,
      reference_trt = reference_trt,...
    )
  } else if (method %in% c("IPTW-Multinomial", "IPTW-Multinomial-Trim", "IPTW-GBM", "IPTW-GBM-Trim", "IPTW-SL", "IPTW-SL-Trim") && estimand == "ATT"&& boot == TRUE) {
    result <- ce_estimate_iptw_att_boot(
      y = y,
      x = x,
      w = w,
      method = method,
      reference_trt = reference_trt,
      nboots = nboots,...
    )
  }
  return(result)
}
