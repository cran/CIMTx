% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/true_c_fun_cal.R
\name{true_c_fun_cal}
\alias{true_c_fun_cal}
\title{Calculate the true c functions with 3 treatments and a binary predictor}
\usage{
true_c_fun_cal(x, w)
}
\arguments{
\item{x}{A matrix with one column for the binary predictor
with values 0 and 1}

\item{w}{A treatment indicator}
}
\value{
A matrix with 2 rows and 6 columns
}
\description{
This function calculates the true confounding functions
with 3 treatments and a binary predictor for simulated data.
}
\examples{
set.seed(111)
data_SA <- data_sim(
  sample_size = 100,
  n_trt = 3,
  x = c(
    "rbinom(1, .5)", # x1:measured confounder
    "rbinom(1, .4)"
  ), # x2:unmeasured confounder
  lp_y = rep(".2*x1+2.3*x2", 3), # parallel response surfaces
  nlp_y = NULL,
  align = FALSE, # w model is not the same as the y model
  lp_w = c(
    "0.2 * x1 + 2.4 * x2", # w = 1
    "-0.3 * x1 - 2.8 * x2"
  ),
  nlp_w = NULL,
  tau = c(-2, 0, 2),
  delta = c(0, 0),
  psi = 1
)
x1 <- data_SA$covariates[, 1, drop = FALSE]
w <- data_SA$w
Y1 <- data_SA$Y_true[, 1]
Y2 <- data_SA$Y_true[, 2]
Y3 <- data_SA$Y_true[, 3]
true_c_fun <- true_c_fun_cal(x = x1, w = w)
}
