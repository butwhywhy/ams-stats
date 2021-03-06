% Generated by roxygen2 (4.0.1): do not edit by hand
\name{NormalErrorGenerator}
\alias{NormalErrorGenerator}
\title{Generates a normal noise generator for use in \code{FakeMeasures} function.}
\usage{
NormalErrorGenerator(sigma)
}
\arguments{
\item{sigma}{Numeric. The standard deviation of the simulated normal noise.}
}
\description{
Generates a normal noise generator for use in \code{FakeMeasures} function.
}
\examples{
# Normal error distribution function for synthetic AMS measures with
# standard debiation 0.1
error.dist <- NormalErrorGenerator(0.1)
}

