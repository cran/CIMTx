% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/causal_multi_treat.R
\name{causal_multi_treat}
\alias{causal_multi_treat}
\title{Estimation of causal effects of multiple treatments}
\usage{
causal_multi_treat(
  y,
  x,
  trt,
  method,
  discard = "No",
  estimand,
  trim_alpha = 0.05,
  SL.library = c("SL.glm", "SL.gam", "SL.knn"),
  reference_trt = 1,
  ndpost = 1000
)
}
\arguments{
\item{y}{numeric vector for the binary outcome}

\item{x}{dataframe including the treatment indicator and the covariates}

\item{trt}{numeric vector for the treatment indicator}

\item{method}{methods for causal inference with multiple treatments. Please select one of the following methods:
\enumerate{
\item {Regression Adjustment: }{Logistics regression to impute missing outcomes}
\item {VM Matching: }{vector matching}
\item {BART: }{Bayesian Additive Regression Trees}
\item {TMLE: }{Targeted maximum likelihood}
\item {IPTW-Logistics: }{Inverse probability of treatment weighting (IPTW) with weights from logistics regression}
\item {IPTW-Logistics-Trim: }{IPTW with trimmed weights from logistics regression}
\item {IPTW-GBM: }{IPTW with weights from generalized boosted method}
\item {IPTW-GBM-Trim: }{IPTW with trimmed weights from generalized boosted method}
\item {IPTW-Superlearner: }{IPTW with weights from superlearner}
\item {IPTW-Superlearner-Trim: }{IPTW with trimmed weights from superlearner}
}}

\item{discard}{discarding rules for BART method. Please select "No", "Lenient" or "Stringent". The default is "No".}

\item{estimand}{causal estimands. Please select "ATT" or "ATE"}

\item{trim_alpha}{alpha values for IPTW weight trimming. The default is 0.05, which means we truncate upper 95\% and lower 5\% of the weights for further IPTW estimation. The default is a combination of SL.glm, SL.gam and SL.knn.}

\item{SL.library}{methods specified with SL.library in Superlearner package}

\item{reference_trt}{Reference group for ATT}

\item{ndpost}{number of independent simulation draws to create}
}
\value{
list with 2 elements for ATT effect. It contains
\item{ATT12:}{A dataframe containing the estimation,
standard error, lower and upper 95\% CI for RD/RR/OR}
\item{ATT13:}{A dataframe containing the estimation,
standard error, lower and upper 95\% CI for RD/RR/OR}
list with 3 elements for ATE effect. It contains
\item{ATE12:}{A dataframe containing the estimation,
standard error, lower and upper 95\% CI for RD/RR/OR}
\item{ATE13:}{A dataframe containing the estimation,
standard error, lower and upper 95\% CI for RD/RR/OR}
\item{ATE23:}{A dataframe containing the estimation,
standard error, lower and upper 95\% CI for RD/RR/OR}
}
\description{
This function estimates the causal effects of multiple treatments with a binary outcome.
}
\examples{
library(CIMTx)
set.seed(3242019)
idata = data_gen(n = 12, ratio =1,scenario = 1)
trt_ind <- as.numeric(idata$trtdat$trt_ind)
all_vars <- idata$trtdat[, -1] #exclude treatment indicator
y <- idata$Yobs

# Regression Adjustment
causal_multi_treat(y = y, x = idata$trtdat,ndpost = 10,
trt = trt_ind, method ="Regression Adjustment", estimand = "ATT", reference_trt = 3)
causal_multi_treat(y = y, x = idata$trtdat,ndpost = 10,
trt = trt_ind, method ="Regression Adjustment",
estimand = "ATE")



# BART with and without discarding
\dontrun{
causal_multi_treat(y = y, x = idata$trtdat,
trt = trt_ind, method = "BART", estimand = "ATE", discard = "No")
causal_multi_treat(y = y, x = idata$trtdat,
trt = trt_ind, method = "BART", estimand = "ATE", discard = "No")
causal_multi_treat(y = y, x = idata$trtdat,
trt = trt_ind, method = "BART", estimand = "ATT", discard = "Stringent")
causal_multi_treat(y = y, x = idata$trtdat,
trt = trt_ind, method = "BART", estimand = "ATT", discard = "Stringent")
causal_multi_treat(y = y, x = idata$trtdat,
trt = trt_ind, method = "BART", estimand = "ATT", discard = "Lenient")
causal_multi_treat(y = y, x = idata$trtdat,
trt = trt_ind, method = "BART", estimand = "ATT", discard = "Lenient")

# VM Matching
causal_multi_treat(y = y, x = idata$trtdat,
trt = trt_ind,method = "VM Matching", estimand = "ATT")

# IPTW-related methods
causal_multi_treat(y = y,trt = trt_ind,
method = "IPTW-Logistics", estimand = "ATT")
causal_multi_treat(y = y,trt = trt_ind,
method = "IPTW-Logistics", estimand = "ATE")
causal_multiple_treatment(y = y, x = idata$trtdat,
trt = trt_ind, method = "IPTW-GBM", estimand = "ATE")
causal_multiple_treatment(y = y, x = idata$trtdat,
trt = trt_ind, method = "IPTW-GBM-Trim", estimand = "ATE")
causal_multiple_treatment(y = y, x = idata$trtdat,
trt = trt_ind, method = "IPTW-Superlearner", estimand = "ATE")
causal_multiple_treatment(y = y, x = idata$trtdat,
trt = trt_ind, method = "IPTW-Superlearner-Trim", estimand = "ATE")
causal_multiple_treatment(y = y, x = idata$trtdat,
trt = trt_ind,method = "IPTW-Superlearner", estimand = "ATT")
causal_multiple_treatment(y = y, x = idata$trtdat,
trt = trt_ind, method = "IPTW-Superlearner-Trim", estimand = "ATT")

}
}
