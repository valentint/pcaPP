\name{objplot}
\alias{objplot}
\title{ Objective Function Plot for Sparse PCs }
\description{
Plots an objective function (TPO or BIC) of one or more sparse PCs for a 
series of lambdas.
}
\usage{
  objplot (x, k, ...)
}
\arguments{
  \item{x}{ An \code{\link{opt.TPO}} or \code{\link{opt.BIC}} object. }
  \item{k}{ This function displays the objective function's values of the 
    \code{k}-th component for \code{opt.TPO}-objects, or the first 
    \code{k} components for \code{opt.BIC}-objects.
  }
  \item{...}{ Further arguments passed to or from other functions. }
}
\details{
This function operates on the return object of function
\code{\link{opt.TPO}} or \code{\link{opt.BIC}}.
The model (\code{lambda}) selected by the minimization of the corresponding 
criterion is highlighted by a dashed vertical line.

The component the argument \code{k} refers to, corresponds to the 
\code{$pc.noord} item of argument \code{x}.
For more info on the order of sparse PCs see the details section of 
\code{\link{opt.TPO}}.
}
%\value{none }
\references{
C. Croux, P. Filzmoser, H. Fritz (2011).
Robust Sparse Principal Component Analysis Based on Projection-Pursuit,
\emph{??} To appear.
}
\author{Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}
\seealso{ \code{\link{sPCAgrid}}, \code{\link{princomp}} }
\examples{

  set.seed (0)
                      ##  generate test data
  x <- data.Zou (n = 250)

  k.max <- 3          ##  max number of considered sparse PCs

                      ##  arguments for the sPCAgrid algorithm
  maxiter <- 25       ##    the maximum number of iterations
  method <- "sd"      ##    using classical estimations

                      ##  Optimizing the TPO criterion
  oTPO <- opt.TPO (x, k.max = k.max, method = method, maxiter = maxiter)

                      ##  Optimizing the BIC criterion
  oBIC <- opt.BIC (x, k.max = k.max, method = method, maxiter = maxiter)

          ##  Objective function vs. lambda
  par (mfrow = c (2, k.max))
  for (i in 1:k.max)        objplot (oTPO, k = i)
  for (i in 1:k.max)        objplot (oBIC, k = i)
}
\keyword{multivariate}
\keyword{robust}
