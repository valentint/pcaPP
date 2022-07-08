\name{opt.TPO}
\alias{opt.TPO}
\alias{opt.BIC}
\title{ Model Selection for Sparse (Robust) Principal Components }
\description{
  These functions compute a suggestion for the sparseness parameter 
  \code{lambda} which is required by function \code{\link{sPCAgrid}}.
  A range of different values for lambda is tested and 
  according to an objective function, the best solution is selected.
  Two different approaches (TPO and BIC) are available, which is further 
  discussed in the details section.
  A graphical summary of the optimization can be obtained by plotting 
  the function's return value (\code{\link{plot.opt.TPO}}, 
  \code{\link{plot.opt.BIC}} for tradeoff curves or \code{\link{objplot}}
  for an objective function plot).
}
\usage{
  opt.TPO (x, k.max = ncol (x), n.lambda = 30, lambda.max, ...)
  opt.BIC (x, k.max = ncol (x), n.lambda = 30, lambda.max, ...)
}
\arguments{
  \item{x}{ a numerical matrix or data frame of dimension (\code{n x p}), 
  which provides the data for the principal components analysis. }
  \item{k.max}{ the maximum number of components which shall be considered
  for optimizing an objective function (optional). }
  \item{n.lambda}{ the number of lambdas to be checked for each component 
  (optional). }
  \item{lambda.max}{ the maximum value of lambda to be checked (optional).
  If omitted, the lambda which yields "full sparseness" (i.e. loadings of 
  only zeros and ones) is computed and used as default value. }
  \item{...}{ further arguments passed to \code{\link{sPCAgrid}} }
}
\details{
% explain difference between TPO and BIC
The choice for a particular lambda is done by optimizing an objective function, 
which is calculated for a set of \code{n.lambda} models with different 
lambdas, ranging from 0 to \code{lambda.max}. If \code{lambda.max} is not 
specified, the minimum lambda yielding "full sparseness" is used.
"Full sparseness" refers to a model with minimum possible absolute sum of 
loadings, which in general implies only zeros and ones in the loadings matrix.

The user can choose between two optimization methods: 
TPO (Tradeoff Product Optimization; see below), or the 
BIC (see Guo et al., 2011; Croux et al., 2011).
The main difference is, that optimization based on the BIC always chooses the 
same lambda for all PCs, and refers to a particular choice of \code{k}, 
the number of considered components.
TPO however is optimized separately for each component, and so delivers 
different lambdas within a model and does not depend on a decision on \code{k}.
\cr
This characteristic can be noticed in the return value of the function:
\code{\link{opt.TPO}} returns a single model with \code{k.max} PCs and 
different values of \code{lambda} for each PC. 
On the contrary \code{\link{opt.BIC}} returns \code{k.max} distinct models 
with \code{k.max} different lambdas, whereas for each model a different
number of components \code{k} has been considered for the optimization.
Applying the latter method, the user finally has to select one of these 
\code{k.max} models manually, 
e.g. by considering the cumulated explained variance,
whereas the TPO method does not require any further decisions.

% TPO 
The tradeoff made in the context of sparse PCA refers to the loss of explained
variance vs. the gain of sparseness. TPO (Tradeoff Product Optimization) 
maximizes the area under the tradeoff curve (see \code{\link{plot.opt.TPO}}),
in particular it maximizes the explained variance multiplied by the number of 
zero loadings of a particular component. As in this context the according 
criterion is minimized, the negative product is considered.

% explain problem of ordered / not ordered PCs
Note that in the context of sparse PCA, there are two sorting orders of PCs, 
which must be considered: Either according to the objective function's value, 
(item \code{$pc.noord})or the variance of each PC(item \code{$pc}).
As in none-sparse PCA the objective function is identical to the PCs' 
variance, this is not an issue there.\cr
The sPCAgrid algorithm delivers the components in decreasing order, according 
to the objective function (which apart from the variance also includes sparseness
terms), whereas the method \code{\link{sPCAgrid}} subsequently re-orders the 
components according to their explained variance.
}
\value{
The functions return an S3 object of type "opt.TPO" or "opt.BIC" respectively, 
containing the following items:
\item{pc}{
An S3 object of type \code{princomp} (\code{\link{opt.TPO}}), or a list of 
  such objects of length \code{k.max} (\code{\link{opt.BIC}}), 
  as returned by \code{\link{sPCAgrid}}.}
\item{pc.noord}{An S3 object of type \code{princomp} (\code{\link{opt.TPO}}), 
  or a list of such objects of length \code{k.max} (\code{\link{opt.BIC}}), 
  as returned by \code{\link{sPCAgrid}}.\cr
  Here the PCs have not been re-ordered according 
  to their variance, but are still ordered according to their objective 
  function's value
  as returned by the \code{\link{sPCAgrid}} - algorithm.
  This information is used for according tradeoff curves and the objective 
  function plot.}
\item{x}{ The input data matrix as provided by the user. }
\item{k.ini, opt}{ These items contain optimization information, as used in 
  functions \code{\link{plot.opt.TPO}}, \code{\link{plot.opt.BIC}} and 
  \code{\link{objplot}}.}
}
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

  oTPO$pc             ##  the model selected by opt.TPO
  oTPO$pc$load        ##  and the according sparse loadings.

                      ##  Optimizing the BIC criterion
  oBIC <- opt.BIC (x, k.max = k.max, method = method, maxiter = maxiter)

  oBIC$pc[[1]]        ##  the first model selected by opt.BIC (k = 1)

          ##  Tradeoff Curves: Explained Variance vs. sparseness
  par (mfrow = c (2, k.max))
  for (i in 1:k.max)        plot (oTPO, k = i)
  for (i in 1:k.max)        plot (oBIC, k = i)

          ##  Tradeoff Curves: Explained Variance vs. lambda
  par (mfrow = c (2, k.max))
  for (i in 1:k.max)        plot (oTPO, k = i, f.x = "lambda")
  for (i in 1:k.max)        plot (oBIC, k = i, f.x = "lambda")

          ##  Objective function vs. lambda
  par (mfrow = c (2, k.max))
  for (i in 1:k.max)        objplot (oTPO, k = i)
  for (i in 1:k.max)        objplot (oBIC, k = i)
}
\keyword{multivariate}
\keyword{robust}
