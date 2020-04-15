
#' OR/RR/RD estimation for ATE
#'
#' The function estimates OR/RR/RD for ATE. Please use our main function causal_multi_treat.R.
#' @param wt1 weight for treatment group 1 in ATE
#' @param wt2 weight for treatment group 2 in ATE
#' @param wt3 weight for treatment group 3 in ATE
#' @param y numeric vector for the binary outcome
#' @param trt_ind numeric vector for the treatment indicator
#'
#' @return list with 3 elements for ATE effect. It contains
#' \item{ATE12:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#' \item{ATE13:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#'\item{ATE23:}{A dataframe containing the estimation,
#'standard error, lower and upper 95\% CI for RD/RR/OR}
#' @export
#' @examples
#'library(CIMTx)
#'ate_fun(wt1 = 1:100, wt2 = 1:100, wt3 = 1:100, y  =1:100, trt_ind =  rep(1:3, c(32,32,36)))
ate_fun <- function(wt1,
                    wt2,
                    wt3,
                    y,
                    trt_ind) {
  mu_1_hat_iptw = sum(y[trt_ind == 1] * wt1[trt_ind == 1]) / sum(wt1[trt_ind == 1])
  mu_2_hat_iptw = sum(y[trt_ind == 2] * wt2[trt_ind == 2]) / sum(wt2[trt_ind == 2])
  mu_3_hat_iptw = sum(y[trt_ind == 3] * wt3[trt_ind == 3]) / sum(wt3[trt_ind == 3])
  RD12 = mu_1_hat_iptw - mu_2_hat_iptw
  RD13 = mu_1_hat_iptw - mu_3_hat_iptw
  RD23 = mu_2_hat_iptw - mu_3_hat_iptw
  RR12 = mu_1_hat_iptw / mu_2_hat_iptw
  RR13 = mu_1_hat_iptw / mu_3_hat_iptw
  RR23 = mu_2_hat_iptw / mu_3_hat_iptw
  OR12 = (mu_1_hat_iptw / (1 - mu_1_hat_iptw)) / (mu_2_hat_iptw / (1 - mu_2_hat_iptw))
  OR13 = (mu_1_hat_iptw / (1 - mu_1_hat_iptw)) / (mu_3_hat_iptw / (1 - mu_3_hat_iptw))
  OR23 = (mu_2_hat_iptw / (1 - mu_2_hat_iptw)) / (mu_3_hat_iptw / (1 - mu_3_hat_iptw))
  res = list(RD12, RD13, RD23, RR12, RR13, RR23, OR12, OR13, OR23)
  return (res)
}
