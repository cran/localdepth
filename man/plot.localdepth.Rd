\name{plot.localdepth}
\title{A DD-plot for Local Depth versus Depth}
\usage{
\method{plot}{localdepth}(x, xlab="Depth", ylab="Local Depth", main="DD plot", mark=0.9, labels=NULL, \dots)
}

\alias{plot.localdepth}

\arguments{
 \item{x}{an object of class \code{localdepth}}
 \item{xlab}{a title for the x axis: see \code{\link{title}}}
 \item{ylab}{a title for the y axis: see \code{\link{title}}}
 \item{main}{an overall title for the plot: see \code{\link{title}}}
 \item{mark}{points with local depth or depth greater than \code{mark} are shown with their row indices in the dataset}
 \item{labels}{labels used to mark the points. It must have the same length as the number of rows of \code{x$y}}
 \item{\dots}{graphical parameters can be given as arguments}
}

\description{
  Provides DD-plot of normalized localdepth versus normalized depth.
}

\value{
  The function returns the \code{x} object as \code{invisible}.
}

\author{Claudio Agostinelli and Mario Romanazzi}

\seealso{\code{\link{localdepth}}}

\examples{
set.seed(1234)
x <- rnorm(20)
tau <- quantile.localdepth(x, probs=0.2, use="volume", method="simplicial")
res <- localdepth(x,tau=tau, use="volume", method="simplicial") 
plot(res, xlab="Simplicial Depth", ylab="Local Simplicial Depth", pch=20)
abline(a=0,b=1,lty="dashed")
}

\keyword{multivariate}
\keyword{nonparametric}


