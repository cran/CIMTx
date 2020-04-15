
#' OR/RR/RD estimation for ATT
#'
#' The function estimates OR/RR/RD for ATT. Please use our main function causal_multi_treat.R.
#'
#' @param wt12 weight for treatment group 2 in ATT
#' @param wt13 weight for treatment group 3 in ATT
#' @param y numeric vector for the binary outcome
#' @param trt_ind numeric vector for the treatment indicator
#'
#' @return list with 2 elements for ATT effect. It contains
#' \item{ATT12:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#'\item{ATT13:}{A dataframe containing the estimation,
#'standard error, lower and upper 95\% CI for RD/RR/OR}
#' @export
#' @examples
#'library(CIMTx)
#'att_fun(wt12 = 1:100, wt13 = 1:100, y  =1:100, trt_ind =  rep(1:3, c(32,32,36)))
att_fun <- function(wt12, wt13,y, trt_ind) {
  mu_1_1_hat_iptw = mean(y[trt_ind == 1])
  mu_1_2_hat_iptw = sum(y[trt_ind == 2] * wt12[trt_ind == 2]) / sum(wt12[trt_ind == 2])
  mu_1_3_hat_iptw = sum(y[trt_ind == 3] * wt13[trt_ind == 3]) / sum(wt13[trt_ind ==  3])
  RD12 = mu_1_1_hat_iptw - mu_1_2_hat_iptw
  RD13  = mu_1_1_hat_iptw - mu_1_3_hat_iptw
  RR12 = mu_1_1_hat_iptw / mu_1_2_hat_iptw
  RR13 = mu_1_1_hat_iptw / mu_1_3_hat_iptw
  OR12 = (mu_1_1_hat_iptw / (1 - mu_1_1_hat_iptw)) / (mu_1_2_hat_iptw / (1 - mu_1_2_hat_iptw))
  OR13 = (mu_1_1_hat_iptw / (1 - mu_1_1_hat_iptw)) / (mu_1_3_hat_iptw / (1 - mu_1_3_hat_iptw))
  res = list(RD12, RD13, RR12, RR13, OR12, OR13)
  return (res)
}
