\name{PCdiagplot}
\alias{PCdiagplot}
\title{Diagnostic plot for principal components}
\description{
  Computes Orthogonal Distances (OD) and Score Distances (SD) for already computed principal components
  using the projection pursuit technique.
}
\usage{
PCdiagplot(x, PCobj, crit = c(0.975, 0.99, 0.999), ksel = NULL, plot = TRUE,
           plotbw = TRUE, raw = FALSE, colgrid = "black", ...)
}
\arguments{
  \item{x}{ a numeric matrix or data frame which provides the data for the
  principal components analysis. }
  \item{PCobj}{ a PCA object resulting from \code{\link{PCAproj}} or \code{\link{PCAgrid}} }
  \item{crit}{ quantile(s) used for the critical value(s) for OD and SD }
  \item{ksel}{ range for the number of PCs to be used in the plot; if NULL all PCs provided are used }
  \item{plot}{ if TRUE a plot is generated, otherwise only the values are returned }
  \item{plotbw}{ if TRUE the plot uses gray, otherwise color representation }
  \item{raw}{ if FALSE, the distribution of the SD will be transformed to approach chisquare 
  distribution, otherwise the raw values are reported and used for plotting }
  \item{colgrid}{ the color used for the grid lines in the plot }
  \item{\dots}{ additional graphics parameters as used in \code{\link{par}} }
}
\details{
  Based on (robust) principal components, a diagnostics plot is made using Orthogonal Distance (OD)
  and Score Distance (SD). This plot can provide important information about the multivariate
  data structure.
}
\value{
  \item{ODist}{matrix with OD for each observation (rows) and each selected PC (cols) }
  \item{SDist}{matrix with SD for each observation (rows) and each selected PC (cols) }
  \item{critOD}{matrix with critical values for OD for each selected PC (rows) and each
        critical value (cols) }
  \item{critSD}{matrix with critical values for SD for each selected PC (rows) and each
        critical value (cols) }
}
\references{ 
  P. Filzmoser and H. Fritz (2007).
  Exploring high-dimensional data with robust principal components.
  In  S. Aivazian, P. Filzmoser, and Yu. Kharin, editors, Proceedings of the Eighth 
  International Conference on Computer Data Analysis and Modeling, volume 1, pp. 43-50, 
  Belarusian State University, Minsk.

  M. Hubert, P.J. Rousseeuwm, K. Vanden Branden (2005).
  ROBCA: a new approach to robust principal component analysis
  Technometrics 47, pp. 64-79.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}
\seealso{\code{\link{PCAproj}}, \code{\link{PCAgrid}}}
\examples{
  # multivariate data with outliers
  library(mvtnorm)
  x <- rbind(rmvnorm(85, rep(0, 6), diag(c(5, rep(1,5)))),
             rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))
  # Here we calculate the principal components with PCAgrid
  pcrob <- PCAgrid(x, k=6)
  resrob <- PCdiagplot(x,pcrob,plotbw=FALSE)

  # compare with classical method:
  pcclass <- PCAgrid(x, k=6, method="sd")
  resclass <- PCdiagplot(x,pcclass,plotbw=FALSE)
}
\keyword{robust}
