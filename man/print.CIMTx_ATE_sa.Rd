% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.CIMTx_ATE_sa.R
\name{print.CIMTx_ATE_sa}
\alias{print.CIMTx_ATE_sa}
\title{Print the ATE results for from sensitivity analysis}
\usage{
\method{print}{CIMTx_ATE_sa}(x, ...)
}
\arguments{
\item{x}{a \code{CIMTx_ATE_sa} object obtained
from \code{\link{sa}} function}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The output from \code{\link{print}}
}
\description{
The \code{\link{print}} method for class "CIMTx_ATE_sa"
}
\examples{
result <- list(estimand = "ATE")
class(result) <- "CIMTx_ATE_sa"
print(result)
}
