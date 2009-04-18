\name{cork}
\alias{cork}
\docType{data}
\title{Cork dataset}
\description{
The cork dataset (Rao, 1948) gives the weights of cork borings of the trunk for 28 trees on the north (N), east (E), south (S) and west (W) directions.
}
\usage{data(cork)}
\format{
Data frame with 28 observations on the following 4 variables.
  \describe{
    \item{\code{N}}{North direction}
    \item{\code{E}}{East direction}
    \item{\code{S}}{South direction}
    \item{\code{W}}{West direction}
  }
}

\source{
  C.R. Rao (1948) Tests of significance in multivariate analysis. Biometrika, 35, 58-79. 
}

\references{
  K.V. Mardia, J.T. Kent and J.M. Bibby (1979) Multivariate Analysis, Academic Press.
}

\examples{
  data(cork)
  set.seed(1234)
  n <- nrow(cork)
  err1 <- rnorm(n,0,sd(cork[,1])/100)
  err2 <- rnorm(n,0,sd(cork[,2])/100)
  err3 <- rnorm(n,0,sd(cork[,3])/100)
  err4 <- rnorm(n,0,sd(cork[,4])/100)
  cork <- cork + data.frame(err1,err2,err3,err4)

  tau <- quantile.localdepth(cork, probs=0.1, size=TRUE, use='volume')
  boxplot(tau$stats)

  depth <- localdepth(x=cork, tau=tau$quantile, tol=0, method='simplicial', use='volume')
  plot(depth)
  abline(0, 1, lty=2)
}
\keyword{datasets}
