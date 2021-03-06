% Generated by roxygen2 (4.0.1): do not edit by hand
\name{SpherEllipse}
\alias{SpherEllipse}
\title{Constructs an spherical ellipse}
\usage{
SpherEllipse(center.dir, axis1.dir, axis2.dir, semiangle1, semiangle2)
}
\arguments{
\item{center.dir}{Numeric vector, representing longitud and latitude of the
center of the ellipse}

\item{axis1.dir}{Numeric vector, representing longitud and latitude of the
first semiaxis of the ellipse}

\item{axis2.dir}{Numeric vector, representing longitud and latitude of the
second semiaxis of the ellipse}

\item{semiangle1}{Numberic value, angular magnitud of the first semiaxis}

\item{semiangle2}{Numberic value, angular magnitud of the second semiaxis}
}
\description{
An ellipse in the sphere is described by a list object 'elip' with constituent
elements 'elip$center', 'elip$eta', 'elip$zeta', 'elip$edir', 'elip$zdir'.
Of these, 'elip$center', 'elip$edir', 'elip$zdir', are vectors containing the
longitud, latitud and radial coordinates of the center direction and the two
semiangles, called eta and zeta. 'elip$eta' and 'elip$zeta' are the value of
the two semiangles. All angular magnitudes in radians.
}
\examples{
center.dir <- c(0, pi/2)
axis1.dir <- c(pi/2, 0)
axis2.dir <- c(pi/2, pi/2)
ellipse <- SpherEllipse(center.dir, axis1.dir, axis2.dir, 1/2, 1)
plot(ellipse)
}
\seealso{
\code{\link{plot.SpherEllipse}}
}

