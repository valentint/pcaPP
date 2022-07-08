\name{l1median_NLM}
\alias{l1median_NM}
\alias{l1median_CG}
\alias{l1median_BFGS}
%\alias{l1median_SA}
\alias{l1median_NLM}
\alias{l1median_HoCr}
\alias{l1median_VaZh}
\title{ Multivariate L1 Median }
\description{
  Computes the multivariate L1 median (also called spatial median) of a data matrix \code{X}.  
}
\usage{
l1median_NM (X, maxit = 200, tol = 10^-8, trace = 0,
             m.init = .colMedians (X), ...)
l1median_CG (X, maxit = 200, tol = 10^-8, trace = 0,
             m.init = .colMedians (X), ...)%, type = 1)
l1median_BFGS (X, maxit = 200, tol = 10^-8, trace = 0,
               m.init = .colMedians (X), REPORT = 10, ...)
l1median_NLM (X, maxit = 200, tol = 10^-8, trace = 0,
              m.init = .colMedians (X), ...)
l1median_HoCr (X, maxit = 200, tol = 10^-8, zero.tol = 1e-15, trace = 0,
               m.init = .colMedians (X), ...)
l1median_VaZh (X, maxit = 200, tol = 10^-8, zero.tol = 1e-15, trace = 0,
	       m.init = .colMedians (X), ...)
}
\arguments{
  \item{X}{ a matrix of dimension \code{n} x \code{p}.}
  \item{maxit}{ The maximum number of iterations to be performed. }
  \item{tol}{ The convergence tolerance. }
  \item{trace}{ The tracing level. Set \code{trace > 0} to retrieve additional information on the single iterations. }
  \item{m.init}{ A vector of length \code{p} containing the initial value of the L1-median. }
%  \item{pscale}{ A vector of length \code{p} containing the variables scale to be used. }
  \item{REPORT}{ A parameter internally passed to \code{\link[stats]{optim}}. }
  \item{zero.tol}{ The zero-tolerance level used in \code{l1median_VaZh} and \code{l1median_HoCr} for determining the equality of two observations (i.e. an observation and a current center estimate). }
  \item{\dots}{ Further parameters passed from other functions. }
%  \item{TMax}{  }
%  \item{TempInit}{  }
}
\value{
%  returns the vector of the coordinates of the L1 median.
%  return (list (par = ret$med, value = ret$dpar.out[1], code = ret$npar.out [1], iterations = ret$npar.out [2], iterations_gr = ret$npar.out [3]))
  \item{par}{ A vector of length \code{p} containing the L1-median. }
  \item{value}{ The value of the objective function \code{||X - l1median||} which is minimized for finding the L1-median. }
  \item{code}{ The return code of the optimization algorithm. See \code{\link[stats]{optim}} and \code{\link[stats]{nlm}} for further information. }
  \item{iterations}{ The number of iterations performed. }
  \item{iterations_gr}{ When using a gradient function this value holds the number of times the gradient had to be computed. }
  \item{time}{The algorithms runtime in milliseconds.}
}
\details{The L1-median is computed using the built-in functions 
    \code{\link[stats]{optim}} and \code{\link[stats]{nlm}}. 
    These functions are a transcript of the \code{L1median} 
    method of package \code{robustX}, porting as much code as possible into C++.}
\author{Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}

\seealso{ \code{\link{median}}}
\note{See the vignette "Compiling pcaPP for Matlab" which comes with this package to compile and use some of these functions in Matlab.}
\examples{

  # multivariate data with outliers
  library(mvtnorm)
  x <- rbind(rmvnorm(200, rep(0, 4), diag(c(1, 1, 2, 2))), 
             rmvnorm( 50, rep(3, 4), diag(rep(2, 4))))

  l1median_NM (x)$par
  l1median_CG (x)$par
  l1median_BFGS (x)$par
  l1median_NLM (x)$par
  l1median_HoCr (x)$par
  l1median_VaZh (x)$par

  # compare with coordinate-wise median:
  apply(x,2,median)
}
\keyword{ multivariate }
\keyword{ robust }
