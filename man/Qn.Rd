\name{qn}
\alias{qn}
\title{ scale estimation using the robust Qn estimator }
\description{
  Returns a scale estimation as calculated by the (robust) Qn estimator.
}
\usage{
qn(x, corrFact)
}
\arguments{
  \item{x}{ a vector of data }
  \item{corrFact}{ the finite sample bias correction factor. By default a value of ~ 2.219144 is used (assuming normality). }
}
\details{
  The Qn estimator computes the first quartile of the pairwise absolute 
  differences of all data values.
}
\value{
  The estimated scale of the data.
}
\references{ 
  P.J. Rousseeuw, C. Croux (1993)
  Alternatives to the Median Absolute Deviation,
  \emph{JASA}, \strong{88}, 1273-1283.
}
\author{ Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}

\section{Warning }{
Earlier implementations used a wrong correction factor for the final result. Thus qn estimations computed with package pcaPP version > 1.8-1 differ about 0.12\% from earlier estimations (version <= 1.8-1).
}

\seealso{ \code{\link{mad}} }
\note{See the vignette "Compiling pcaPP for Matlab" which comes with this package to compile and use this function in Matlab.}
\examples{
  # data with outliers
  x <- c(rnorm(100), rnorm(10, 10))
  qn(x)
}
\keyword{multivariate}
\keyword{robust}
