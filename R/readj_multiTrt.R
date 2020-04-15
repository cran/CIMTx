
#' Regression Adjustment
#'
#' This function implements the regression adjustment method. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param trt numeric vector for the treatment indicator
#' @param estimand causal estimands. Please select "ATT" or "ATE"
#' @param ndpost number of independent simulation draws to create
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
#' library(CIMTx)
#' set.seed(3242019)
#' idata = data_gen(n = 12, ratio =1,scenario = 1)
#' trt_ind <- as.numeric(idata$trtdat$trt_ind)
#' all_vars <- idata$trtdat[, -1]
#' y <- idata$Yobs
#' regadj_multiTrt(y = y, x = idata$trtdat,trt = trt_ind, estimand="ATE")

regadj_multiTrt = function(y, x, trt, estimand="ATE", ndpost=1000) {
  x <- x[, -1]
  # Data structure
  #        Y(1) Y(2) Y(3)
  # trt=1   *    ?    ?
  # trt=2   ?    *    ?
  # trt=3   ?    ?    *

  #        Y(1) Y(2) Y(3)
  # trt=1  y11  y12  y13
  # trt=2  y21  y22  y23
  # trt=3  y31  y32  y33

  if (estimand=="ATE") {
    regadj_est = regadj_multiTrt_ate(y, x, trt, ndpost=1000)
  }

  if (estimand=="ATT") {
    regadj_est = regadj_multiTrt_att(y, x, trt, ndpost=1000)
  }

  return(regadj_est)
}
