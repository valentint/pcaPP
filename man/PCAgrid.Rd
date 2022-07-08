\name{PCAgrid}
\alias{PCAgrid}
\alias{sPCAgrid}
\title{ (Sparse) Robust Principal Components using the Grid search algorithm }
\description{
  Computes a desired number of (sparse) (robust) principal components using 
  the grid search algorithm in the plane.
  The global optimum of the objective function is searched in planes, not in
  the p-dimensional space, using regular grids in these planes.
}
\usage{
PCAgrid (x, k = 2, method = c ("mad", "sd", "qn"),
         maxiter = 10, splitcircle = 25, scores = TRUE, zero.tol = 1e-16, 
	 center = l1median, scale, trace = 0, store.call = TRUE, control, ...)

sPCAgrid (x, k = 2, method = c ("mad", "sd", "qn"), lambda = 1,
          maxiter = 10, splitcircle = 25, scores = TRUE, zero.tol = 1e-16, 
	  center = l1median, scale, trace = 0, store.call = TRUE, control, ...)
}
\arguments{
  \item{x}{ a numerical matrix or data frame of dimension (\code{n x p})which 
  provides the data for the principal components analysis. }
  \item{k}{ the desired number of components to compute }
  \item{method}{ the scale estimator used to detect the direction with the 
  largest variance. Possible values are \code{"sd"}, \code{"mad"} and 
  \code{"qn"}, the  latter can be called \code{"Qn"} too. \code{"mad"} is the 
  default value.}
  \item{lambda}{ the sparseness constraint's strength(\code{sPCAgrid} only).
    A single value for all components, or a vector of length \code{k} with 
    different values for each component can be specified.
    See \code{\link{opt.TPO}} for the choice of this argument. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{splitcircle}{ the number of directions in which the algorithm should
  search for the largest variance. The direction with the largest variance
  is searched for in the directions defined by a number of equally spaced points
  on the unit circle. This argument determines, how many such points are used to
  split the unit circle. }
  \item{scores}{ A logical value indicating whether the scores of the
  principal component should be calculated. }
  \item{zero.tol}{ the zero tolerance used internally for checking 
  convergence, etc. }
%  \item{anglehalving}{ boolean stating whether angle halving is to be used or not. 
%  Angle halving will usually improve the solution quite a lot.}
%  \item{fact2dim}{ an integer that is multiplied to splitcircle if x is only 
%  two-dimensional. In higher dimensions, fewer search directions are needed to allow
%  for faster computation. In two dimensions, more search directions are required to 
%  grant higher precision. \code{fact2dim} is used to take account of this.}
  \item{center}{ this argument indicates how the data is to be centered. It
  can be a function like \code{\link{mean}} or \code{\link{median}} or a vector
  of length \code{ncol(x)} containing the center value of each column. }
  \item{scale}{ this argument indicates how the data is to be rescaled. It
  can be a function like \code{\link{sd}} or \code{\link{mad}} or a vector
  of length \code{ncol(x)} containing the scale value of each column. }
  \item{trace}{ an integer value >= 0, specifying the tracing level. }
%  \item{cut.pc}{ a logical value, specifying whether only the first \code{k} 
%  columns of the resulting loadings and scores matrix shall be returned 
%  (\code{TRUE}). If this value is \cpode{FALSE}, the algorithm returns 
%  an \code{p x p} loadings- and an \code{n x p} scores matrix, whereas the
%  last \code{p - k} components form an arbitrary basis of the complementary
%  space of the first \code{k} found components.}
%  \item{pc.ini}{ an optional pre calculated \code{princomp} (S3) object. If
%  provided, the algorithm searches for additional PCs in the orthogonal space
%  of the components provided by this object. }
%  \item{k.ini}{ an optional integer value, specifying how many components of
%  \code{pc.ini} shall be considered. }
%  \item{ord.all}{ a logical value, specifying wheter}
%  \item{HDred}{ }
  \item{store.call}{ a logical variable, specifying whether the function call
    shall be stored in the result structure. }
  \item{control}{ a list which elements must be the same as (or a subset of) 
  the parameters above. If the control object is supplied, the parameters from 
  it will be used and any other given parameters are overridden. }
  \item{...}{ further arguments passed to or from other functions. }
}
\details{
  In contrast to \code{PCAgrid}, the function \code{sPCAgrid} computes sparse
  principal components. The strength of the applied sparseness constraint is
  specified by argument \code{lambda}.
  %Setting \code{lambda = 0} yields the same estimation for both functions 
  %\code{PCAgrid} and\code{sPCAgrid}.

  Similar to the function \code{\link{princomp}}, there is a \code{print} method 
  for the these objects that prints the results in a nice format and the 
  \code{plot} method produces a scree plot (\code{\link{screeplot}}). There is 
  also a \code{biplot} method.

  Angle halving is an extension of the original algorithm. In the original 
  algorithm, the search directions are determined by a number of points on the
  unit circle in the interval [-pi/2 ; pi/2). Angle halving means this angle is
  halved in each iteration, eg. for the first approximation, the above mentioned
  angle is used, for the second approximation, the angle is halved to 
  [-pi/4 ; pi/4) and so on. This usually gives better results with less 
  iterations needed. \cr
  NOTE: in previous implementations angle halving could be suppressed by the 
  former argument "\code{anglehalving}". This still can be done by setting 
  argument \code{maxiter = 0}.
}
\value{
  The function returns an object of class \code{"princomp"}, i.e. a list 
  similar to the output of the function \code{\link{princomp}}.
  \item{sdev}{the (robust) standard deviations of the principal components.}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix whose columns
  contain the eigenvectors).  This is of class \code{"loadings"}:
  see \code{\link{loadings}} for its \code{\link{print}} method.}
  \item{center}{the means that were subtracted.}
  \item{scale}{the scalings applied to each variable.}
  \item{n.obs}{the number of observations.}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied data on the
  principal components.}
  \item{call}{the matched call.}
  \item{obj}{A vector containing the objective functions values. For function
  \code{PCAgrid} this is the same as \code{sdev}. }
  \item{lambda}{The lambda each component has been calculated with 
  (\code{\link{sPCAgrid}} only).}
}
\references{
  C. Croux, P. Filzmoser, M. Oliveira, (2007).
  Algorithms for Projection-Pursuit Robust Principal Component Analysis,
  \emph{Chemometrics and Intelligent Laboratory Systems}, Vol. 87, pp. 218-225.

  C. Croux, P. Filzmoser, H. Fritz (2011).
  Robust Sparse Principal Component Analysis Based on Projection-Pursuit,
  \emph{??} To appear.
}
\author{Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}
\seealso{ \code{\link{PCAproj}}, \code{\link{princomp}} }
\note{See the vignette "Compiling pcaPP for Matlab" which comes with this package to compile and use these functions in Matlab.}
\examples{
  # multivariate data with outliers
  library(mvtnorm)
  x <- rbind(rmvnorm(200, rep(0, 6), diag(c(5, rep(1,5)))),
             rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))
  # Here we calculate the principal components with PCAgrid
  pc <- PCAgrid(x)
  # we could draw a biplot too:
  biplot(pc)
  # now we want to compare the results with the non-robust principal components
  pc <- princomp(x)
  # again, a biplot for comparison:
  biplot(pc)

  ##  Sparse loadings
  set.seed (0)
  x <- data.Zou ()

                   ##  applying PCA
  pc <-  princomp (x)
                   ##  the corresponding non-sparse loadings
  unclass (pc$load[,1:3])
  pc$sdev[1:3]

                   ##  lambda as calculated in the opt.TPO - example
  lambda <- c (0.23, 0.34, 0.005)
                   ##  applying sparse PCA
  spc <- sPCAgrid (x, k = 3, lambda = lambda, method = "sd")
  unclass (spc$load)
  spc$sdev[1:3]

                   ## comparing the non-sparse and sparse biplot
  par (mfrow = 1:2)
  biplot (pc, main = "non-sparse PCs")
  biplot (spc, main = "sparse PCs")
}
\keyword{multivariate}
\keyword{robust}
