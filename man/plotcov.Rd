\name{plotcov}
\alias{plotcov}
\title{ Compare two Covariance Matrices in Plots}
\description{
  allows a direct comparison of two estimations of the covariance 
  matrix (e.g. resulting from covPC) in a plot.
}
\usage{
plotcov(cov1, cov2, method1, labels1, method2, labels2, ndigits, ...)
}
\arguments{
  \item{cov1}{ a covariance matrix (from cov, covMcd, covPC, 
     covPCAgrid, covPCAproj, etc. }
  \item{cov2}{ a covariance matrix (from cov, covMcd, covPC, 
     covPCAgrid, covPCAproj, etc. }
  \item{method1}{ legend for ellipses of estimation method1}
  \item{method2}{ legend for ellipses of estimation method2}
  \item{labels1}{ legend for numbers of estimation method1}
  \item{labels2}{ legend for numbers of estimation method2}
  \item{ndigits}{ number of digits to use for printing covariances, by default 
    ndigits=4}
  \item{...}{ additional arguments for text or plot}
}
\details{
  Since (robust) PCA can be used to re-compute the (robust) covariance matrix,
  one might be interested to compare two different methods of covariance
  estimation visually. This routine takes as input objects for the covariances
  to compare the output of \code{\link{cov}}, but also the return objects 
  from \code{\link{covPCAgrid}}, \code{\link{covPCAproj}}, \code{\link{covPC}},
  and \code{\link[robustbase]{covMcd}}. The comparison of the two covariance matrices
  is done by numbers (the covariances) and by ellipses.
}
\value{
  only the plot is generated
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
  plotcov(covPCAproj(x),covPCAgrid(x))

}
\keyword{ multivariate }
