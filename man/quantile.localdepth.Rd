\name{quantile.localdepth}
\alias{quantile.localdepth}
\title{Quantiles of the distribution of the size of simplices or ellipsoids}
\description{
The function evaluates the size of the objects formed from a dataset according to a statistic and returns the corresponding quantiles. Depending on the method the objects are simplices or ellipsoids.
}
\usage{
quantile.localdepth(x, probs, use = c("volume", "diameter"), 
method = c("simplicial", "ellipsoid", "halfspace", "mahalanobis"), 
nsamp = "all", size = FALSE, dimension = NULL, covariance=NULL, \dots)
}
\arguments{
  \item{x}{numeric; vector, dataframe or matrix}
  \item{probs}{numeric; vector of probabilities with values in [0,1]. Quantile orders of the statistic}
  \item{use}{character; the statistic used to measure the size of the objects. For \code{method} equal to \code{"simplicial"} or \code{"ellipsoid"} allowed statistics are \code{"volume"} and \code{"diameter"}. For \code{method} equal to \code{"mahalanobis"} this parameter is not used and the only available statistic is pairwise Mahalanobis' distance}
  \item{method}{character; the type of (local) depth to be evaluated}
  \item{nsamp}{character or numeric; the number of objects that are considered. If \code{"all"}, the size of all \code{choose(NROW(x), NCOL(x)+1)} objects is evaluated. Otherwise, a simple random sample with replacement of size \code{nsamp} is performed from the set of all possible objects. See details below}
  \item{size}{logical; if \code{TRUE} the size of all possible objects is returned}
  \item{dimension}{numeric; only used with \code{method="ellipsoid"}. It is the squared length of the ellipsoid semimajor axis. If \code{dimension} is \code{NULL}, it is set to \code{NCOL(x)}}
  \item{covariance}{\code{NULL} or a numeric matrix; the \code{NCOL(x)*NCOL(x)} covariance matrix used in \code{method} equal to \code{"mahalanobis"}. If \code{NULL}, \code{cov(x)} is used}
  \item{\dots}{arguments passed to \code{\link{quantile.default}}}
}

\details{
When choosing between \code{nsamp='all'} or a numeric value for \code{nsamp} the value of \code{choose(NROW(x), NCOL(x)+1)} must be considered to avoid very long computing times. When the number of possible objects is greater than the maximum vector length manageable by the machine, an error is returned by the system. In this case it is necessary to use the approximated procedure through the parameter \code{nsamp}.
}

\value{
If \code{size} is \code{FALSE} then a vector with the quantiles is returned, otherwise a list with the following components:
  \item{quantile}{the requested quantiles}
  \item{stats}{the size of all objects}
  \item{call}{match.call()}
}
\author{Claudio Agostinelli and Mario Romanazzi}
\seealso{\code{\link{localdepth}}}
\examples{

set.seed(1234)
x <- matrix(rnorm(60, 0, 1), ncol=2)
volumesimplex <- quantile.localdepth(x, probs=c(0.1, 0.2), size=TRUE)
volumesimplex$quantile
diametersimplex <- quantile.localdepth(x, probs=c(0.1, 0.2), size=TRUE, use='diameter')
diametersimplex$quantile
par(mfrow=c(2,2))
plot(ecdf(volumesimplex$stats), xlab='volume', main="ECDF of the simplices's volume")
boxplot(volumesimplex$stats, xlab='volume', main="Boxplot of the simplices's volume")
plot(ecdf(diametersimplex$stats), xlab='volume', main="ECDF of the simplices's diameter")
boxplot(diametersimplex$stats, xlab='volume', main="Boxplot of the simplices's diameter")

}

\keyword{multivariate}
\keyword{nonparametric}
