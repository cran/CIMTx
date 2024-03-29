% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ce_estimate_vm_att.R
\name{ce_estimate_vm_att}
\alias{ce_estimate_vm_att}
\title{Causal inference with multiple treatments using VM for ATT effects}
\usage{
ce_estimate_vm_att(y, x, w, reference_trt, caliper, n_cluster)
}
\arguments{
\item{y}{A numeric vector (0, 1) representing a binary outcome.}

\item{x}{A dataframe, including all the covariates but not treatments.}

\item{w}{A numeric vector representing the treatment groups.}

\item{reference_trt}{A numeric value indicating reference treatment group
for ATT effect.}

\item{caliper}{A numeric value denoting the caliper on the logit of
GPS within each cluster formed by K-means clustering.
The caliper is in standardized units.
For example, \code{caliper = 0.25} means that
all matches greater than 0.25 standard deviations of the
logit of GPS are dropped. The default value is 0.25.}

\item{n_cluster}{A numeric value denoting the number of clusters to
form using K means clustering on the logit of GPS.}
}
\value{
A summary of the effect estimates can be obtained
with \code{summary} function. The output also contains the number
of matched individuals.
}
\description{
The function \code{ce_estimate_vm_att} implements
VM to estimate ATT effect with
multiple treatments using observational data.
}
\references{
Venables, W. N. & Ripley, B. D. (2002)
\emph{Modern Applied Statistics with S}.
Fourth Edition. Springer, New York. ISBN 0-387-95457-0

Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021).
\emph{dplyr: A Grammar of Data Manipulation}.
R package version 1.0.7.
URL: \url{https://CRAN.R-project.org/package=dplyr}

Jasjeet S. Sekhon (2011).
Multivariate and Propensity Score Matching Software with A
utomated Balance Optimization: The Matching Package for R.
\emph{Journal of Statistical Software}, \strong{42}(7), 1-52
}
