% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.CIMTx_sa_grid.R
\name{print.CIMTx_sa_grid}
\alias{print.CIMTx_sa_grid}
\title{Print the ATT results for from sensitivity analysis}
\usage{
\method{print}{CIMTx_sa_grid}(x, ...)
}
\arguments{
\item{x}{a \code{CIMTx_sa_grid} object obtained from
\code{\link{sa}} function}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The output from \code{\link{print}}
}
\description{
The \code{\link{print}} method for class "CIMTx_sa_grid"
}
\examples{
result <- list(estimand = "ATE")
class(result) <- "CIMTx_sa_grid"
print(result)
}
