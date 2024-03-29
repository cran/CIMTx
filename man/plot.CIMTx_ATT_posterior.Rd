% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.CIMTx_ATT_posterior.R
\name{plot.CIMTx_ATT_posterior}
\alias{plot.CIMTx_ATT_posterior}
\title{Plot for non-IPTW estimation methods for ATT effect}
\usage{
\method{plot}{CIMTx_ATT_posterior}(x, ...)
}
\arguments{
\item{x}{a \code{CIMTx_ATT_posterior} object obtained}

\item{...}{further arguments passed to or from other methods.}
}
\value{
an error message
}
\description{
Plot for non-IPTW estimation methods for ATT effect
}
\examples{
\dontrun{
result <- list(method = "RA")
class(result) <- "CIMTx_ATT_posterior"
plot(result)
}
}
