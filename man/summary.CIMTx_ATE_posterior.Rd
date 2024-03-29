% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.CIMTx_ATE_posterior.R
\name{summary.CIMTx_ATE_posterior}
\alias{summary.CIMTx_ATE_posterior}
\title{Summarize a CIMTx_ATE_posterior object}
\usage{
\method{summary}{CIMTx_ATE_posterior}(object, ...)
}
\arguments{
\item{object}{a \code{CIMTx_ATE_posterior} object
obtained with \code{\link{ce_estimate}} function.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
a list with w*(w-1)/2 elements for ATE effect.
Each element of the list contains the estimation,
standard error, lower and upper 95\% CI for RD/RR/OR.
}
\description{
Summarize a CIMTx_ATE_posterior object
}
\examples{
library(CIMTx)
lp_w_all <-
  c(
    ".4*x1 + .1*x2  - .1*x4 + .1*x5", # w = 1
    ".2 * x1 + .2 * x2  - .2 * x4 - .3 * x5"
  ) # w = 2
nlp_w_all <-
  c(
    "-.5*x1*x4  - .1*x2*x5", # w = 1
    "-.3*x1*x4 + .2*x2*x5"
  ) # w = 2
lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - .1*x4 - .2*x5", 3)
nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
X_all <- c(
  "rnorm(0, 0.5)", # x1
  "rbeta(2, .4)", # x2
  "runif(0, 0.5)", # x3
  "rweibull(1,2)", # x4
  "rbinom(1, .4)" # x5
)

set.seed(111111)
data <- data_sim(
  sample_size = 300,
  n_trt = 3,
  x = X_all,
  lp_y = lp_y_all,
  nlp_y = nlp_y_all,
  align = FALSE,
  lp_w = lp_w_all,
  nlp_w = nlp_w_all,
  tau = c(-1.5, 0, 1.5),
  delta = c(0.5, 0.5),
  psi = 1
)
ce_estimate_ra_ate_result <- ce_estimate(
  y = data$y, x = data$covariates,
  w = data$w, ndpost = 10, method = "RA", estimand = "ATE"
)
summary(ce_estimate_ra_ate_result)
}
\references{
Hadley Wickham (2019).
\emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
R package version 1.4.0.
URL:\url{https://CRAN.R-project.org/package=stringr}
}
