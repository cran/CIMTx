% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.CIMTx_sa_grid.R
\name{plot.CIMTx_sa_grid}
\alias{plot.CIMTx_sa_grid}
\title{Contour plot for the grid specification of sensitivity analysis}
\usage{
\method{plot}{CIMTx_sa_grid}(x, ate = NULL, att = NULL, ...)
}
\arguments{
\item{x}{Object from sa function}

\item{ate}{a character indicating the ATE effect to plot,
eg, "1,3" or "2,3"}

\item{att}{a character indicating the ATT effect to plot,
eg, "1,3" or "1,2"}

\item{...}{further arguments passed to or from other methods.}
}
\value{
A ggplot figure
}
\description{
This function make the countor plot after the grid specification
of sensitivity analysis.
The input of the function is from the output of the sa.R function.
}
\examples{
sa_object_example <- list(
  ATE13 = seq(0, 1, length.out = 25), grid_index = c(4, 5),
  c_functions = data.frame(
    c4 = rep(seq(-0.6, 0, 0.15), each = 5),
    c5 = rep(seq(0, 0.6, 0.15), 5)
  )
)
class(sa_object_example) <- "CIMTx_sa_grid"
plot(sa_object_example, ate = "1,3")
}
\references{
Hadley Wickham (2019).
\emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
R package version 1.4.0.
URL:\url{https://CRAN.R-project.org/package=stringr}

Elio Campitelli (2021).
 \emph{metR: Tools for Easier Analysis of Meteorological Fields}.
  R package version 0.11.0.
  URL:\url{https://github.com/eliocamp/metR}
}
