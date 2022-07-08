\name{l1median}
\alias{l1median}
\title{ Multivariate L1 Median }
\description{
  Computes the multivariate L1 median (also called spatial median) of a data matrix.  
}
\usage{
l1median(X, MaxStep = 200, ItTol = 10^-8, trace = 0, m.init = .colMedians (X))
}
\arguments{
  \item{X}{ A matrix containing the values whose multivariate L1 median is to be computed. }
  \item{MaxStep}{ The maximum number of iterations. }
  \item{ItTol}{ Tolerance for convergence of the algorithm.}
  \item{trace}{ The tracing level. }
  \item{m.init}{ An initial estimate. }
}
\value{
  returns the vector of the coordinates of the L1 median.
}
\references{
  C. Croux, P. Filzmoser, M. Oliveira, (2007).
  Algorithms for Projection-Pursuit Robust Principal Component Analysis,
  \emph{Chemometrics and Intelligent Laboratory Systems}, Vol. 87, pp. 218-225.
}
\author{Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}

\seealso{ \code{\link{median}}}
\examples{
  l1median(rnorm(100), trace = -1) # this returns the median of the sample

  # multivariate data with outliers
  library(mvtnorm)
  x <- rbind(rmvnorm(200, rep(0, 4), diag(c(1, 1, 2, 2))), 
             rmvnorm( 50, rep(3, 4), diag(rep(2, 4))))
  l1median(x, trace = -1)
  # compare with coordinate-wise median:
  apply(x,2,median)
}
\keyword{ multivariate }
\keyword{ robust }
