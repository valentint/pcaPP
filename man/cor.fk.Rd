\name{cor.fk}
\alias{cor.fk}
\title{ Fast estimation of Kendall's tau rank correlation coefficient}
\description{
  Calculates Kendall's tau rank correlation coefficient in O (n log (n)) 
  rather than O (n^2) as in the current R implementation.
}
\usage{
%cor.fk (x, y = NULL, cor = TRUE)
cor.fk (x, y = NULL)
}
\arguments{
  \item{x}{ A vector, a matrix or a data frame of data. }
  \item{y}{ A vector of data. }
%  \item{cor}{ A logical value, specifying whether the correlation (\code{TRUE}) or the covariance (\code{FALSE} shall be estimated. )}
}
\value{
  The estimated correlation coefficient.
}
\details{
The code of this implementation of the fast Kendall's tau correlation algorithm has 
originally been published by David Simcha. 
Due to it's runtime (\code{O(n log n)}) it's essentially faster than the 
current R implementation (\code{O (n\^2)}), especially for large numbers of 
observations.
The algorithm goes back to Knight (1966) and has been described more detailed 
by Abrevaya (1999) and Christensen (2005).
}
\references{ 
  Knight, W. R. (1966). A Computer Method for Calculating Kendall's Tau with Ungrouped Data.
    Journal of the American Statistical Association, \strong{314}(61) Part 1, 436-439.\cr
  Christensen D. (2005). Fast algorithms for the calculation of Kendall's Tau.
    Journal of Computational Statistics \strong{20}, 51-62.\cr
  Abrevaya J. (1999). Computation of the Maximum Rank Correlation Estimator.
    Economic Letters \strong{62}, 279-285.
}
\author{ David Simcha, Heinrich Fritz, Christophe Croux, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}

\seealso{ \code{\link{cor}} }
\examples{

  set.seed (100)		## creating test data
  n <- 1000
  x <- rnorm (n)
  y <- x+  rnorm (n)

  tim <- proc.time ()[1]	## applying cor.fk
  cor.fk (x, y)
  cat ("cor.fk runtime [s]:", proc.time ()[1] - tim, "(n =", length (x), ")\n")

  tim <- proc.time ()[1]	## applying cor (standard R implementation)
  cor (x, y, method = "kendall")
  cat ("cor runtime [s]:", proc.time ()[1] - tim, "(n =", length (x), ")\n")

		##	applying cor and cor.fk on data containing				

  Xt <- cbind (c (x, as.integer (x)), c (y, as.integer (y)))

  tim <- proc.time ()[1]	## applying cor.fk
  cor.fk (Xt)
  cat ("cor.fk runtime [s]:", proc.time ()[1] - tim, "(n =", nrow (Xt), ")\n")

  tim <- proc.time ()[1]	## applying cor (standard R implementation)
  cor (Xt, method = "kendall")
  cat ("cor runtime [s]:", proc.time ()[1] - tim, "(n =", nrow (Xt), ")\n")
}
\keyword{multivariate}
\keyword{robust}
