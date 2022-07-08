\name{plot.opt.TPO}
\alias{plot.opt.TPO}
\alias{plot.opt.BIC}
\title{ Tradeoff Curves for Sparse PCs }
\description{
Tradeoff curves of one or more sparse PCs for a series of lambdas, which 
contrast the loss of explained variance and the gain of sparseness.
}
\usage{
\method{plot}{opt.TPO} (x, k, f.x = c ("l0", "pl0", "l1", "pl1", "lambda"),
                        f.y = c ("var", "pvar"), ...)
\method{plot}{opt.BIC} (x, k, f.x = c ("l0", "pl0", "l1", "pl1", "lambda"),
                        f.y = c ("var", "pvar"), ...)
}
\arguments{
  \item{x}{ An \code{\link{opt.TPO}} or \code{\link{opt.BIC}} object. }
  \item{k}{ This function plots the tradeoff curve of the 
    \code{k}-th component for \code{opt.TPO}-objects, or the first 
    \code{k} components for \code{opt.BIC}-objects.
  }
  \item{f.x, f.y}{ A string, specifying which information shall be plotted on the 
  x and y - axis. See the details section for more information.}
  \item{...}{ Further arguments passed to or from other functions. }
}
\details{
The argument \code{f.x} can obtain the following values:
\itemize{
  \item \code{"l0"}: l0 - sparseness, which corresponds to the number of 
    zero loadings of the considered component(s).
  \item \code{"pl0"}: l0 - sparseness in percent (l0 - sparseness 
    ranges from \code{0} to \code{p-1} for each component).
  \item \code{"l1"}: l1 - sparseness, which corresponds to 
    the negative sum of absolute 
    loadings of the considered component(s).\cr
    (The exact value displayed for a single component is
    \code{sqrt (p) - S}, with \code{S} as the the absolute sum of loadings.)
    \cr
    As this value is a part of the objective function which selects 
    the candidate directions within the \code{\link{sPCAgrid}} function,
    this option is provided here.
    \item \code{"pl1"} The "l1 - sparseness" in percent (l1 - sparseness 
    ranges from \code{0} to \code{sqrt (p-1)} for each component).
  \item \code{"lambda"}: The lambda used for computing a particular model.
}
The argument \code{f.y} can obtain the following values:
\itemize{
  \item \code{"var"}: The (cumulated) explained variance of the considered
  component(s). The value shown here is calculated using the variance 
  estimator specified via the \code{method} argument of function 
  \code{\link{sPCAgrid}}.
  \item \code{"pvar"}: The (cumulated) explained variance of the considered
  component(s) in percent. The 100\%-level is assumed as the sum of variances 
  of all columns of argument \code{x}.\cr
  Again the same variance estimator is 
  used as specified via the \code{method} argument of function 
  \code{\link{sPCAgrid}}.\cr
}

The subtitle summarizes the result of the applied criterion for selecting a
value of lambda:
\itemize{
\item The name of the applied method (TPO/BIC).
\item The selected value of \code{lambda} for the \code{k}-th component 
  (\code{\link{opt.TPO}}) or all computed components (\code{\link{opt.BIC}}).
\item The empirical cumulated variance (ECV) of the first \code{k} components 
  in percent.
\item The obtained l0-sparseness of the first \code{k} components.
}

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

          ##  Tradeoff Curves: Explained Variance vs. sparseness
  par (mfrow = c (2, k.max))
  for (i in 1:k.max)        plot (oTPO, k = i)
  for (i in 1:k.max)        plot (oBIC, k = i)

          ##  Explained Variance vs. lambda
  par (mfrow = c (2, k.max))
  for (i in 1:k.max)        plot (oTPO, k = i, f.x = "lambda")
  for (i in 1:k.max)        plot (oBIC, k = i, f.x = "lambda")
}
\keyword{multivariate}
\keyword{robust}
