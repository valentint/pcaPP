\name{PCAproj}
\alias{PCAproj}
\title{ Robust Principal Components using the algorithm of Croux and Ruiz-Gazen (2005) }
\description{
  Computes a desired number of (robust) principal components using the algorithm of 
  Croux and Ruiz-Gazen (JMVA, 2005). 
}
\usage{
PCAproj(x, k = 2, method = c("mad", "sd", "qn"), CalcMethod = c("eachobs",
"lincomb", "sphere"), nmax = 1000, update = TRUE, scores = TRUE, maxit = 5, 
maxhalf = 5, scale = NULL, center = l1median_NLM, zero.tol = 1e-16, control)
}
\arguments{
  \item{x}{ a numeric matrix or data frame which provides the data for the 
  principal components analysis. }
  \item{k}{ desired number of components to compute }
  \item{method}{ scale estimator used to detect the direction with the largest 
  variance. Possible values are \code{"sd"}, \code{"mad"} and \code{"qn"}, the 
  latter can be called \code{"Qn"} too. \code{"mad"} is the default value.}
  \item{CalcMethod}{ the variant of the algorithm to be used. Possible values are
  \code{"eachobs"}, \code{"lincomb"} and \code{"sphere"}, with \code{"eachobs"} being 
  the default. }
  \item{nmax}{ maximum number of directions to search in each step (only when
  using \code{"sphere"} or \code{"lincomb"} as the \code{CalcMethod}). }
  \item{update}{ a logical value indicating whether an update algorithm should be 
  used.}
  \item{scores}{ a logical value indicating whether the scores of the
  principal component should be calculated. }
  \item{maxit}{ maximim number of iterations. }
  \item{maxhalf}{ maximum number of steps for angle halving. }
  \item{scale}{ this argument indicates how the data is to be rescaled. It
  can be a function like \code{\link{sd}} or \code{\link{mad}} or a vector
  of length \code{ncol(x)} containing the scale value of each column. }
  \item{center}{ this argument indicates how the data is to be centered. It
  can be a function like \code{\link{mean}} or \code{\link{median}} or a vector
  of length \code{ncol(x)} containing the center value of each column. }
  \item{zero.tol}{ the zero tolerance used internally for checking 
  convergence, etc. }
  \item{control}{ a list which elements must be the same as (or a subset of) 
  the parameters above. If the control object is supplied, the parameters from 
  it will be used and any other given parameters are overridden. }
}
\details{
  Basically, this algrithm considers the directions of each observation 
  through the origin of the centered data as possible projection directions.
  As this algorithm has some drawbacks, especially if \code{ncol(x) > nrow(x)}
  in the data matrix, there are several improvements that can be used with this
  algorithm.
\itemize{
  \item update - An updating step basing on the algorithm for finding the
  eigenvectors is added to the algorithm. This can be used with any 
  \code{CalcMethod}
  \item sphere - Additional search directions are added using random directions.
  The random directions are determined using random data points generated from 
  a p-dimensional multivariate standard normal distribution. These new data 
  points are projected to the unit sphere, giving the new search directions.
  \item lincomb - Additional search directions are added using linear 
  combinations of the observations. It is similar to the 
  \code{"sphere"} - algorithm, but the new data points are generated using linear 
  combinations of the original data \code{b_1*x_1 + ... + b_n*x_n} where the 
  coefficients \code{b_i} come from a uniform distribution in the interval 
  \code{[0, 1]}.
}
  Similar to the function \code{\link{princomp}}, there is a \code{print} method 
  for the these objects that prints the results in a nice format and the \code{plot} 
  method produces a scree plot (\code{\link{screeplot}}). There is also a 
  \code{\link{biplot}} method.
}
\value{
  The function returns a list of class \code{"princomp"}, i.e. a list similar to the 
  output of the function \code{\link{princomp}}.
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
}
\references{
  C. Croux, P. Filzmoser, M. Oliveira, (2007).
  Algorithms for Projection-Pursuit Robust Principal Component Analysis,
  \emph{Chemometrics and Intelligent Laboratory Systems}, Vol. 87, pp. 218-225.
}
\author{ Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}
\seealso{ \code{\link{PCAgrid}}, \code{\link{ScaleAdv}}, \code{\link{princomp}} }
\examples{
  # multivariate data with outliers
  library(mvtnorm)
  x <- rbind(rmvnorm(200, rep(0, 6), diag(c(5, rep(1,5)))),
             rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))
  # Here we calculate the principal components with PCAgrid
  pc <- PCAproj(x, 6)
  # we could draw a biplot too:
  biplot(pc)

  # we could use another calculation method and another objective function, and 
  # maybe only calculate the first three principal components:
  pc <- PCAproj(x, 3, "qn", "sphere")
  biplot(pc)

  # now we want to compare the results with the non-robust principal components
  pc <- princomp(x)
  # again, a biplot for comparision:
  biplot(pc)
}
\keyword{ robust }
\keyword{ multivariate }
