\name{data.Zou}
\alias{data.Zou}
\title{ Test Data Generation for Sparse PCA examples }
\description{
Draws a sample data set, as introduced by Zou et al. (2006).
}
\usage{
data.Zou (n = 250, p =  c(4, 4, 2), ...)
}
\arguments{
  \item{n}{ The required number of observations. }
  \item{p}{ A vector of length 3, specifying how many variables shall be 
  constructed using the three factors V1, V2 and V3. }
  \item{...}{ Further arguments passed to or from other functions. }
}
\details{
This data set has been introduced by Zou et al. (2006), and then been referred 
to several times, e.g. by Farcomeni (2009), Guo et al. (2010) and 
Croux et al. (2011).

The data set contains two latent factors V1 ~ N(0, 290) and V2 ~ N(0, 300) and 
a third mixed component V3 = -0.3 V1 + 0.925V2 + e; e ~ N(0, 1).\cr
The ten variables Xi of the original data set are constructed the following 
way:\cr
Xi = V1 + ei; i = 1, 2, 3, 4\cr
Xi = V2 + ei; i = 5, 6, 7, 8\cr
Xi = V3 + ei; i = 9, 10\cr
whereas ei ~ N(0, 1) is indepependent for i = 1 , ..., 10
}
\value{
 A matrix of dimension \code{n x sum (p)} containing the generated sample data 
 set.
 }
\references{
C. Croux, P. Filzmoser, H. Fritz (2011).
Robust Sparse Principal Component Analysis Based on Projection-Pursuit,
\emph{??} To appear.

A. Farcomeni (2009).
An exact approach to sparse principal component analysis,
\emph{Computational Statistics}, Vol. 24(4), pp. 583-604.

J. Guo, G.  James, E. Levina, F. Michailidis, and J. Zhu (2010).
Principal component analysis with sparse fused loadings,
\emph{Journal of Computational and Graphical Statistics.} To appear.

H. Zou, T. Hastie, R. Tibshirani (2006). 
Sparse principal component analysis,
\emph{Journal of Computational and Graphical Statistics}, Vol. 15(2), pp. 265-286.
}
\author{Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}
\seealso{ \code{\link{sPCAgrid}}, \code{\link{princomp}} }
\examples{
                   ##  data generation
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
