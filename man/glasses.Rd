\name{glasses}
\alias{glasses}
\docType{data}
\title{Roman-Venetian Glass dataset}
\description{
  The dataset gives the compositions of 73 glass specimens excavated in
  a series of archeological campaign in the island of Venice lagoon
  and nearby borders.
}
\usage{data(glasses)}
\format{
Data frame with 73 observations on the following 14 variables.
  \describe{
    \item{\code{ID}}{Specimens number}
    \item{\code{CLA}}{type of glass: 1, Venetian glass (XIVth-XVIth
      centuries); 2, Roman glass (Ist-Vth centuries A.D.);
      3, Middle Ages glass (VIth-XIIIth centuries)}
    \item{\code{SIO2}}{a numeric vector}
    \item{\code{AL2O3}}{a numeric vector}
    \item{\code{NA2O}}{a numeric vector}
    \item{\code{K2O}}{a numeric vector}
    \item{\code{CAO}}{a numeric vector}
    \item{\code{MGO}}{a numeric vector}
    \item{\code{FE2O3}}{a numeric vector}
    \item{\code{TIO2}}{a numeric vector}
    \item{\code{P2O5}}{a numeric vector}
    \item{\code{CL}}{a numeric vector}
    \item{\code{SO3}}{a numeric vector}
    \item{\code{MNO}}{a numeric vector}
  }
}
\details{
  Variables SIO2 to MNO are provide composition of glass specimens
  on the corresponding chemical compound.
}
\source{
Mario Romanazzi.
}
\references{
  M. Verita', A. Renier and S. Zecchin (2002) Chemical analyses of ancient
  glass findings excavated in the venetian lagoon. Journal of Cultural
  Heritage, 3, 261-271.
}
\examples{
glasses <- na.omit(glasses[,-1])
n <- nrow(glasses)
p <- ncol(glasses)
col <- rep("black",n)
col[glasses[,1] == 2] <- "red"
col[glasses[,1] == 3] <- "green"
lab <- rep("M",n)
lab[glasses[,1] == 1] <- "V"
lab[glasses[,1] == 2] <- "R"
lab1 <- paste(lab,1:n,sep="")
pc <- princomp(glasses[,-1],cor=TRUE)
summary(pc)
plot(scale(pc$scores[,1]), scale(pc$scores[,2]), pch=20, xlab="PC1 (36)",
  ylab="PC2 (28)", main="Roman Venetian Glasses Data")
text(scale(pc$scores[,1]), scale(pc$scores[,2]), labels=lab1, cex=0.6,
  pos=3, col=col)
abline(h=0, v=0, lty="dashed", col="grey")
plot(scale(pc$scores[,1]), scale(pc$scores[,3]), pch=20, xlab="PC1 (36)",
  ylab="PC3 (11)", main="Roman Venetian Glasses Data")
text(scale(pc$scores[,1]), scale(pc$scores[,3]), labels=lab1, cex=0.6,
  pos=3, col=col)
abline(h=0, v=0, lty="dashed", col="grey")
plot(scale(pc$scores[,2]), scale(pc$scores[,3]), pch=20, xlab="PC2 (28)",
  ylab="PC3 (11)", main="Roman Venetian Glasses Data")
text(scale(pc$scores[,2]), scale(pc$scores[,3]), labels=lab1, cex=0.6,
  pos=3, col=col)
abline(h=0, v=0, lty="dashed", col="grey")

set.seed(1234)
err1 <- rnorm(n,0,sd(glasses[,2])/100)
err2 <- rnorm(n,0,sd(glasses[,3])/100)
err3 <- rnorm(n,0,sd(glasses[,4])/100)
err4 <- rnorm(n,0,sd(glasses[,5])/100)
err5 <- rnorm(n,0,sd(glasses[,6])/100)
err6 <- rnorm(n,0,sd(glasses[,7])/100)
err7 <- rnorm(n,0,sd(glasses[,8])/100)
err8 <- rnorm(n,0,sd(glasses[,9])/100)
err9 <- rnorm(n,0,sd(glasses[,10])/100)
err10 <- rnorm(n,0,sd(glasses[,11])/100)
err11 <- rnorm(n,0,sd(glasses[,12])/100)
err12 <- rnorm(n,0,sd(glasses[,13])/100)
glasses <- glasses[,-1] + data.frame(err1,err2,err3,err4,err5,err6,err7,
  err8,err9, err10,err11,err12)

tau <- quantile.localdepth(glasses, probs=c(0.1, 0.9), method='mahalanobis')
gla10 <- localdepth(glasses, tau=tau[1], method='mahalanobis')
gla90 <- localdepth(glasses, tau=tau[2], method='mahalanobis')

plot(gla10)
abline(0, 1, lty=2)
plot(gla90)
abline(0, 1, lty=2)

}
\keyword{datasets}
