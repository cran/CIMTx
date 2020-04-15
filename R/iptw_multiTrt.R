
#' Inverse probability of treatment weighting (IPTW)
#'
#' This function implements the IPTW method. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param trt numeric vector for the treatment indicator
#' @param psdat data frame containing the treatment indicator and covariates
#' @param estimand causal estimands, "ATT" or "ATE"
#' @param method methods for causal inference with multiple treatments, inherited from causal_multi_treat.R
#' @param wt1 weight for treatment group 1 in ATE
#' @param wt2 weight for treatment group 2 in ATE
#' @param wt3 weight for treatment group 3 in ATE
#' @param wt12 weight for treatment group 2 in ATT
#' @param wt13 weight for treatment group 3 in ATT
#' @param trim_alpha alpha values for IPTW weight trimming, inherited from causal_multi_treat.R
#' @param SL.library methods specified with SL.library in Superlearner package, inherited from causal_multi_treat.R
#'
#' @return list with 2 elements for ATT effect. It contains
#' \item{ATT12:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#'\item{ATT13:}{A dataframe containing the estimation,
#'standard error, lower and upper 95\% CI for RD/RR/OR}
#' list with 3 elements for ATE effect. It contains
#' \item{ATE12:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#' \item{ATE13:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#'\item{ATE23:}{A dataframe containing the estimation,
#'standard error, lower and upper 95\% CI for RD/RR/OR}
#' @export
#' @examples
#'library(CIMTx)
#'set.seed(1)
#'idata = data_gen(n = 50, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'iptw_multiTrt(y=y, trt = trt_ind,SL.library = c("SL.glm"),
#'trim_alpha = 0.05, method = "IPTW-GBM", estimand = "ATT")

iptw_multiTrt = function(y, trt, psdat, estimand = "ATE", method, wt1, wt2, wt3,wt12, wt13, trim_alpha,SL.library) {
  if (estimand == "ATE") {
    iptw_est = iptw_multiTrt_ate(y, trt, psdat, wt1, wt2, wt3, method, trim_alpha,SL.library)
  }
  if (estimand == "ATT") {
    iptw_est = iptw_multiTrt_att(y, trt, psdat, wt12, wt13,method, trim_alpha,SL.library)
  }
  return(iptw_est)
}
