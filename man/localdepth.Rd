\name{localdepth}
\alias{localdepth}
\title{Local depth}
\description{
  The function evaluates the depth and the local depth for a set of
  points with respect to a dataset. 
}
\usage{
localdepth(x, y = NULL, tau, use = c("volume", "diameter"), 
  method = c("simplicial", "ellipsoid",
  "halfspace", "mahalanobis", "hyperspheresimplicial"), 
  type = c("exact", "approx"), nsamp = "all", nmax = 1, 
  tol=10^(-9), dimension=NULL, location=NULL, covariance=NULL)
}

\arguments{
  \item{x}{numeric; vector, dataframe or matrix. If \code{x} is a
    \code{circular} vector, a circular version is used. Avoid ties
    by wiggling the data. The function only issues a warning for ties}
  \item{y}{numeric; vector, dataframe or matrix with the same number of
    columns as \code{x}, or \code{NULL}. If \code{NULL}, \code{x} is used}
  \item{tau}{numeric; threshold value for the evaluation of the local
    depth. Use function \code{\link{quantile.localdepth}} to evaluate
    \code{tau} using a quantile of the size of the objects}
  \item{use}{character; the statistic used to measure the size of the
    objects. Currently, for \code{method} equal to \code{"simplicial"}
    or \code{"ellipsoid"} allowed statistics are \code{"volume"} and
    \code{"diameter"}. For \code{method} equal to \code{"mahalanobis"}
    this parameter is not used and the only available statistic is
    pairwise Mahalanobis' distance}
  \item{method}{character; the type of (local) depth to be evaluated}
  \item{type}{character; how to evaluate membership. Only active for
    \code{method="simplicial"}. See details.}
  \item{nsamp}{character or numeric; the number of objects that are
    considered. If \code{"all"}, the size of all \code{choose(NROW(x),
      NCOL(x)+1)} objects is evaluated. Otherwise, a simple random
    sample with replacement of size \code{nsamp} is performed from
    the set of all possible objects.}
  \item{nmax}{numeric; maximum fraction (in the range (0,1]) of objects
    to be considered when \code{nsamp} is not equal to \code{all}.
    If \code{nmax=1} the number of searched objects can reach the number
    of possible objects (\code{choose(NROW(x), NCOL(x)+1)} for
    simplicial and ellipsoid depth)}
  \item{tol}{numeric; tolerance parameter to be fixed depending on the
    machine precision. Used to decide membership of points located near
    to the boundary of the objects}
  \item{dimension}{numeric; only used with \code{method="ellipsoid"}.
    It is the squared length of the ellipsoid semimajor axis.
    If \code{dimension} is \code{NULL}, it is set to \code{NCOL(x)}}
  \item{location}{\code{NULL} or a numeric vector; the \code{NCOL(x)}
    means vector used in \code{method} equal to \code{"mahalanobis"}.
    If \code{NULL}, \code{apply(x, 2, mean)} is used}
  \item{covariance}{\code{NULL} or a numeric matrix; the
    \code{NCOL(x)*NCOL(x)} covariance matrix used in \code{method} equal
    to \code{"mahalanobis"}. If \code{NULL}, \code{cov(x)} is used}
}

\details{
With \code{method="simplicial"} and \code{type="exact"}, membership of
the points in simplices is evaluated; when \code{type="approx"},
an approximate membership function is used. See references below.

With \code{method="hyperspheresimplicial"} the data matrix are normalize
to have vectors of length one, then the (local) angular simplicial depth
is evaluated. 
}

\value{
  The function returns an object of \code{\link{class}}
  \code{localdepth} with the following components:

  \item{localdepth}{vector of the local depth values for the given \code{tau}}
  \item{depth}{vector of the depth values}
  \item{max.localdepth}{\code{max(localdepth)}}
  \item{max.depth}{\code{max(depth)}}
  \item{num}{vector with two components. \code{num[1]} gives the number
    of objects used for the evaluation of the depth; \code{num[2]} is
    the number of objects used for the evaluation of the local depth}
  \item{call}{match.call() result. Note that this is called from the
    internal function}
  \item{tau}{value of the corresponding input parameter}
  \item{use}{value of the corresponding input parameter}
  \item{tol}{value of the corresponding input parameter}
  \item{x}{value of the corresponding input parameter}
  \item{y}{value of the corresponding input parameter}
  \item{type}{value of the corresponding input parameter}
  \item{nsamp}{value of the corresponding input parameter}
  \item{method}{value of the corresponding input parameter}
}

\references{
  C. Agostinelli and M. Romanazzi (2007). Local depth of univariate
  distributions. Working paper n. 1/2007, Dipartimento di Statistica,
  Universita' Ca' Foscari, Venezia.

  C. Agostinelli and M. Romanazzi (2008). Local depth of multidimensional
  data. Working paper n. 3/2008, Dipartimento di Statistica,
  Universita' Ca' Foscari, Venezia.

  C. Agostinelli and M. Romanazzi (2011). Local depth.
  Journal of Statistical Planning and Inference, 141:817-830.
  doi:10.1016/j.jspi.2010.08.001.

  C. Agostinelli e M. Romanazzi (2013). Nonparametric analysis of
  directional data based on data depth.
  Environmental and Ecological Statistics, 20(2):253-270.
  doi:10.1007/s10651-012-0218-z.
}

\author{Claudio Agostinelli and Mario Romanazzi}

\seealso{\code{\link{quantile.localdepth}}, \code{\link{plot.localdepth}}}
\examples{
# Iris Setosa dataset
data(iris)
setosa <- iris[iris$Species=="setosa",-5]
str(setosa)
n <- dim(setosa)[1];p <- dim(setosa)[2]
lab <- paste("U",1:n,sep="")
pairs(setosa)
# Wiggling data with gaussian error
set.seed(1234)
err <- matrix(c(rnorm(n,0,sd(setosa[,1])/1000),
  rnorm(n,0,sd(setosa[,2])/1000),
  rnorm(n,0,sd(setosa[,3])/1000),
  rnorm(n,0,sd(setosa[,4])/1000)),n,p,byrow=FALSE)
setosa1 <- setosa + as.data.frame(err)

\dontrun{
qd <- quantile.localdepth(setosa1, probs=c(seq(0.01,0.09,0.01),
  seq(0.1,0.9,0.1)), use="diameter", method="simplicial", nsamp="all", size=TRUE)
qv <- quantile.localdepth(setosa1, probs=c(seq(0.01,0.09,0.01),
  seq(0.1,0.9,0.1)), use="volume", method="simplicial", nsamp="all",size=TRUE)
}

## Here we explore 10% of simplices (0.1*choose(50,4+1))
qd <- quantile.localdepth(setosa1,
  probs=c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1)),
  use="diameter", method="simplicial", nsamp=211876, size=TRUE)
qv <- quantile.localdepth(setosa1,
  probs=c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1)),
  use="volume", method="simplicial", nsamp=211876,size=TRUE)

c(mean(qd$stats),median(qd$stats))
c(mean(qv$stats),median(qv$stats))

round(qd$quantile,2)
round(qv$quantile,8)

boxplot(qd$stats,names="Diameter",xlab="Simplex Size",
  horizontal=TRUE,main="Iris Setosa")
boxplot(qv$stats,names="Volume",xlab="Simplex Size",
  horizontal=TRUE,main="Iris Setosa")

plot(qd$stats,qv$stats,xlab="Simplex Diameter",
  ylab="Simplex Volume", main="Iris Setosa",pch=20)

# Simplicial global/local depth

sldv10 <- localdepth(setosa1,tau=qv$quantile[10],use="volume",
  method="simplicial", type="approx", nsamp="all", nmax=0.2)
str(sldv10)
sd10 <- sldv10$depth       # simplicial depth
sld10 <- sldv10$localdepth # simplicial local depth (volume, 10%)
lab[sd10 == (p+1)/n]        # convex hull vertices (25 points out of 50)
sdmax <- max(sd10);imax <- which.max(sd10);setosa[imax,] # deepest point

# DD-Plot

plot(sldv10, xlab="Simplicial Depth",
  ylab="Local Simplicial Depth", main="Iris Setosa",pch=20,ylim=c(0,1.02))
abline(a=0,b=1,lty="dashed")

\dontrun{
# Mahalanobis Depth of Iris Versicolor wrt Iris Setosa

iversicolor <- iris[iris$Species == "versicolor",-5]
str(iversicolor)
set.seed(2345)
err <- matrix(c(rnorm(n,0,sd(iversicolor[,1])/1000),
  rnorm(n,0,sd(iversicolor[,2])/1000),
  rnorm(n,0,sd(iversicolor[,3])/1000),
  rnorm(n,0,sd(iversicolor[,4])/1000)), n,p,byrow=FALSE)
iversicolor1 <- iversicolor + as.data.frame(err)
qm <- quantile.localdepth(setosa1,
  probs=c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1)),
  method="mahalanobis",nsamp="all",size=TRUE)
round(qm$quantile,2)
vers <- localdepth(x=setosa1, y=iversicolor1, tau=qm$quantile[10],
  method="mahalanobis", nsamp="all")
str(vers)
round(vers$depth,6)
round(vers$localdepth,6) # identically zero
}
}

\keyword{multivariate}
\keyword{nonparametric}
