\name{chaetocnema}
\alias{chaetocnema}
\docType{data}
\title{Chaetocnema dataset}
\description{
The Chaetocnema dataset includes anatomical measures of 74 males of fleabeetles, of three different genera: Ch. Concinna (A, 21 cases), Ch. Heikertingeri (B, 31 cases) and Ch. Heptapotamica (C, 22 cases).
}
\usage{data(chaetocnema)}
\format{
Data frame with 74 observations on the following 8 variables.
  \describe{
    \item{\code{X10}}{width of the first joint of the first tarsus in microns (the sum of measurements for both tarsi)}
    \item{\code{X12}}{the same for the second joint}
    \item{\code{X14}}{the maximal width of the aedeagus in the fore-part in microns}
    \item{\code{X18}}{the front angle of the aedeagus (1 unit = 7.5 degrees)}
    \item{\code{X40}}{the maximal width of the head between teh external edges of the eyes in 0.01mm}
    \item{\code{X48}}{the aedeagus width from the side in microns}
    \item{\code{species}}{species of the fleabeetles. A factor with levels \code{A}: Ch. Concinna, \code{B}: Ch. Heikertingeri, \code{C}: Ch. Heptapotamica}
    \item{\code{area}}{colletion area of specimens. A factor with levels \code{A}: Environs of Uljianovsk, \code{B}: Khvalynsk, the Volga, \code{C}: Pern, \code{D}: Environs of Leningrad, \code{E}: the Ukraine, \code{F}: Ashkhabd, Turkmenistan, \code{G}: France, \code{H}: Ustj-Zilma, \code{I}: Gabra, Abkhazia, \code{L}: Ussuri district, \code{M}: Yakutsk district, \code{N}: Khabarovsk, \code{O}: Germany, \code{P}: Environs of Lake Issyk-Kul, Kirghizia, \code{Q}: Alma-ata, Kazachstan, \code{R}: Environs of Frunza, Kirghizia}
  }
}
\details{
The aedeagus is the median lobe of the male genital capsule that is surrounded by the phallobase and through which the sperm duct traverses. 
}
\source{
A.A. Lubischew (1962) On the use of discriminant functions in taxonomy. Biometrics, 18:4, 455--477.
}
\references{
Hand, D.J., et al. (1994) A Handbook of Small Data Sets, London: Chapman & Hall, 254-255.
}

\seealso{
The dataset is also present in \code{\link[amap]{lubisch}}.
}

\examples{
data(chaetocnema)
lab <- chaetocnema[,7]
chaetocnema <- chaetocnema[,4:5]
set.seed(1234)
chaetosd <- apply(chaetocnema, 2, sd)
chaetocnema <- chaetocnema+cbind(rnorm(nrow(chaetocnema), 0, chaetosd[1]/100), rnorm(nrow(chaetocnema), 0, chaetosd[2]/100))


gx <- seq(min(chaetocnema[,1]), max(chaetocnema[,1]), length=60)
gy <- seq(min(chaetocnema[,2]), max(chaetocnema[,2]), length=60)
y <- as.matrix(expand.grid(gx, gy))

tauvolume <- quantile.localdepth(chaetocnema, probs=0.20, use='volume')
system.time(depthcontours <- localdepth(x=chaetocnema, y=y, tau=tauvolume, use='volume', method='simplicial'))
system.time(depthcontoursapprox <- localdepth(x=chaetocnema, y=y, tau=tauvolume, use='volume', method='simplicial', type='approx'))

par(mfrow=c(2,2))
plot(chaetocnema, pch=as.character(lab), main='Exact Simplicial Depth')
contour(x=gx, y=gy, z=matrix(depthcontours$depth/depthcontours$max.depth, nrow=length(gx), ncol=length(gy)), add=TRUE)

plot(chaetocnema, pch=as.character(lab), main='Approx Simplicial Depth')
contour(x=gx, y=gy, z=matrix(depthcontoursapprox$depth/depthcontoursapprox$max.depth, nrow=length(gx), ncol=length(gy)), add=TRUE)

plot(chaetocnema, pch=as.character(lab), main='Exact Simplicial Local Depth, Tau=0.2')
contour(x=gx, y=gy, z=matrix(depthcontours$localdepth/depthcontours$max.localdepth, nrow=length(gx), ncol=length(gy)), add=TRUE)

plot(chaetocnema, pch=as.character(lab), main='Approx Simplicial Local Depth, Tau=0.2')
contour(x=gx, y=gy, z=matrix(depthcontoursapprox$localdepth/depthcontoursapprox$max.localdepth, nrow=length(gx), ncol=length(gy)), add=TRUE)
}
\keyword{datasets}
