\name{covPC}
\alias{covPC}
\title{ Covariance Matrix Estimation from princomp Object}
\description{
  computes the covariance matrix from a princomp object. 
  The number of components k can be given as input.
}
\usage{
covPC(x, k, method)
}
\arguments{
  \item{x}{ an object of class princomp. }
  \item{k}{ number of PCs to use for covariance estimation (optional). }
  \item{method}{ method how the PCs have been estimated (optional). }
}
\details{
  There are several possibilities to estimate the principal components (PCs)
  from an input data matrix, including the functions \code{\link{PCAproj}} and 
  \code{\link{PCAgrid}}. This function uses the estimated PCs to reconstruct
  the covariance matrix. Not all PCs have to be used, the number k of
  PCs (first k PCs) can be given as input to the function.
}
\value{
  \item{cov}{ the estimated covariance matrix}
  \item{center}{ the center of the data, as provided from the princomp object.}
  \item{method}{ a string describing the method that was used to calculate the
  PCs.}
}
\references{
  C. Croux, P. Filzmoser, M. Oliveira, (2007).
  Algorithms for Projection-Pursuit Robust Principal Component Analysis,
  \emph{Chemometrics and Intelligent Laboratory Systems}, Vol. 87, pp. 218-225.
}
\author{ Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}> }
\seealso{ \code{\link{PCAgrid}}, \code{\link{PCAproj}}, \code{\link{princomp}} }
\examples{
  # multivariate data with outliers
  library(mvtnorm)
  x <- rbind(rmvnorm(200, rep(0, 6), diag(c(5, rep(1,5)))),
             rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))
  pc <- princomp(x)
  covPC(pc, k=2)
}
\keyword{ multivariate }
