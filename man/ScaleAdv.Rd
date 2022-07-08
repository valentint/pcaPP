\name{ScaleAdv}
\alias{ScaleAdv}
\title{ centers and rescales data }
\description{
  Data is centered and rescaled (to have mean 0 and a standard deviation of 1).
}
\usage{
ScaleAdv(x, center = mean, scale = sd)
}
\arguments{
  \item{x}{ matrix containing the observations. If this is not a matrix, but
  a data frame, it is automatically converted into a matrix using the function
  \code{\link{as.matrix}}. In any other case, (eg. a vector) it is converted 
  into a matrix with one single column. }
  \item{center}{ this argument indicates how the data is to be centered. It 
  can be a function like \code{\link{mean}} or \code{\link{median}} or a vector
  of length \code{ncol(x)} containing the center value of each column. }
  \item{scale}{ this argument indicates how the data is to be rescaled. It 
  can be a function like \code{\link{sd}} or \code{\link{mad}} or a vector
  of length \code{ncol(x)} containing the scale value of each column. }
}
\value{
  The function returns a list containing
  \item{ x }{centered and rescaled data matrix.}
  \item{ center }{ a vector of the centers of each column x. If you add to 
  each column of \code{x} the appropriate value from \code{center}, you will obtain
  the data with the original location of the observations. }
  \item{ scale }{ a vector of the scale factors of each column x. If you multiply 
  each column of \code{x} by the appropriate value from \code{scale}, you will obtain
  the data with the original scales. }
}
\details{
  The default \code{scale} being \code{NULL} means that no rescaling is done.
}
\references{
  C. Croux, P. Filzmoser, M. Oliveira, (2007).
  Algorithms for Projection-Pursuit Robust Principal Component Analysis,
  \emph{Chemometrics and Intelligent Laboratory Systems}, Vol. 87, pp. 218-225.
}
\author{ Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}
\examples{
  x <- rnorm(100, 10, 5)
  x <- ScaleAdv(x)$x

  # can be used with multivariate data too
  library(mvtnorm)
  x <- rmvnorm(100, 3:7, diag((7:3)^2))
  res <- ScaleAdv(x, center = l1median, scale = mad)
  res

  # instead of using an estimator, you could specify the center and scale yourself too
  x <- rmvnorm(100, 3:7, diag((7:3)^2))
  res <- ScaleAdv(x, 3:7, 7:3)
  res
}
\keyword{ multivariate }
