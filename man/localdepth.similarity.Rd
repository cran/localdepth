\name{localdepth.similarity}
\alias{localdepth.similarity}

\title{Local depth similarity}
\description{
The function evaluates depth and local depth similarity for a set of points with respect to a dataset. 
}
\usage{
localdepth.similarity(x, y = NULL, tau, use = c("volume", "diameter"), 
  method = c("simplicial", "ellipsoid", "mahalanobis"), 
  type = c("exact", "approx"), nsamp = "all", nmax = 1, 
  tol = 10^(-9), dimension=NULL, location = NULL, covariance = NULL, 
  weight = NULL)
}

\arguments{
  \item{x}{numeric; vector, dataframe or matrix. If \code{x} is a \code{circular} vector, a circular version is used. Avoid ties by wiggling the data. The function only issues a warning for ties.}
  \item{y}{numeric; vector, dataframe or matrix with the same number of columns as \code{x}, or \code{NULL}. If \code{NULL}, \code{x} is used}
  \item{tau}{numeric; threshold value for the evaluation of the local depth. Use function \code{\link{quantile.localdepth}} to evaluate \code{tau} using a quantile of the size of the objects}
  \item{use}{character; the statistic used to measure the size of the objects. Currently, for \code{method} equal to \code{"simplicial"} or \code{"ellipsoid"} allowed statistics are \code{"volume"} and \code{"diameter"}. For \code{method} equal to \code{"mahalanobis"} this parameter is not used and the only available statistic is pairwise Mahalanobis' distance}
  \item{method}{character; the type of (local) depth similarity to be evaluated}
  \item{type}{character; how to evaluate membership. Only active for \code{method="simplicial"}. See details.}
  \item{nsamp}{character or numeric; the number of objects that are considered. If \code{"all"}, the size of all \code{choose(NROW(x), NCOL(x)+1)} objects is evaluated. Otherwise, a simple random sample with replacement of size \code{nsamp} is performed from the set of all possible objects.}
  \item{nmax}{numeric; maximum fraction (in the range (0,1]) of objects to be considered when \code{nsamp} is not equal to \code{all}. If \code{nmax=1} the number of searched objects can reach the number of possible objects (\code{choose(NROW(x), NCOL(x)+1)} for simplicial and ellipsoid depth)}
  \item{tol}{numeric; tolerance parameter to be fixed depending on the machine precision. Used to decide membership of points located near to the boundary of the objects}
  \item{dimension}{numeric; only used with \code{method="ellipsoid"}. It is the squared length of the ellipsoid semimajor axis. If \code{dimension} is \code{NULL}, it is set to \code{NCOL(x)}}
  \item{location}{\code{NULL} or a numeric vector; the \code{NCOL(x)} means vector used in \code{method} equal to \code{"mahalanobis"}. If \code{NULL}, \code{apply(x, 2, mean)} is used}
  \item{covariance}{\code{NULL} or a numeric matrix; the \code{NCOL(x)*NCOL(x)} covariance matrix used in \code{method} equal to \code{"mahalanobis"}. If \code{NULL}, \code{cov(x)} is used}
  \item{weight}{experimental parameter used to weight entries in the similarity matrix. Not implemented in each method, dimension.}
}
\details{
With \code{method="simplicial"} and \code{type="exact"}, membership of the points in simplices is evaluated; when \code{type="approx"}, an approximate membership function is used. See references below.
}

\value{
  The function returns an object of \code{\link{class}} \code{localdepth.similarity} with the following components:

  \item{localdepth}{matrix of the local depth similarities}
  \item{depth}{matrix of the depth similarities}
  \item{max.localdepth}{\code{max(localdepth)}}
  \item{max.depth}{\code{max(depth)}}
  \item{num}{vector with two components. \code{num[1]} gives the number of objects used for the evaluation of the depth similarity; \code{num[2]} is the number of objects used for the evaluation of the local depth similarity}
  \item{call}{match.call() result. Note that this is called from the internal function}
  \item{tau}{value of the corresponding input parameter}
  \item{use}{value of the corresponding input parameter}
  \item{tol}{value of the corresponding input parameter}
  \item{x}{value of the corresponding input parameter}
  \item{y}{value of the corresponding input parameter}
  \item{type}{value of the corresponding input parameter}
  \item{nsamp}{value of the corresponding input parameter}
  \item{method}{value of the corresponding input parameter}
}

\note{The function is not yet implemented for Ellipsoid (local) depth.}

\author{Claudio Agostinelli and Mario Romanazzi}
\references{
C. Agostinelli and M. Romanazzi (2007). Local depth of univariate distributions. Working paper n. 1/2007, Dipartimento di Statistica, Universita' Ca' Foscari, Venezia.

C. Agostinelli and M. Romanazzi (2008). Local depth of multidimensional data. Working paper n. 3/2008, Dipartimento di Statistica, Universita' Ca' Foscari, Venezia.

  R.Y. Liu, J.M. Parelius and K. Singh (1999) Multivariate analysis by data depth: descriptive statistics, graphics and inference. The Annals of Statistics, 27, 783-858.
}
\seealso{\code{\link{localdepth}}}
\examples{
  data(cork)
  tau <- quantile.localdepth(cork[,c(1,3)], probs=0.1, method='simplicial')
  sim <- localdepth.similarity(cork[,c(1,3)], tau=tau, method='simplicial')
  plot(hclust(d=as.dist(1-sim$localdepth/sim$max.localdepth)))
  plot(hclust(d=as.dist(1-sim$depth/sim$max.depth)))
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{cluster}
