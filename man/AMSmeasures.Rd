% Generated by roxygen2 (4.0.1): do not edit by hand
\name{AMSmeasures}
\alias{AMSmeasures}
\title{Constructs an AMS measurements object}
\usage{
AMSmeasures(repetitions, positions, values)
}
\arguments{
\item{repetitions}{Integer. It must have the same cardinality as the
\code{values} and represents the number of repetition of the
correponding measurement.}

\item{positions}{Integer vector or factor. It must have the same
cardinality as the \code{vaues} parameter, and represents
the order of the corresponding measurement in the experimental
setup.}

\item{values}{Numberic vector.The measured values.}
}
\description{
This method constructs an object with class \code{AMSmeasures}, which uses
a \code{data.frame} internally.
}
\examples{
reps <- rep(1, 15)
poss <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
values <- c(0.86,1.26,1.10,0.86,1.26,1.26,1.26,1.02,1.26,1.26,0.80,1.80,
            1.50,0.80,1.80)
# Construct AMS measures object
AMSmeasures(reps, poss, values)
}
\seealso{
\code{\link{as.AMSmeasures}}; \code{\link{reps.AMSmeasures}},
    \code{\link{poss.AMSmeasures}}, \code{\link{values.AMSmeasures}}
}

