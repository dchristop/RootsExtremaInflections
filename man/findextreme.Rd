\name{findextreme}
\alias{findextreme}
\title{
Implementation of Integration Extreme Finding Estimator (IEFE) for extreme points identification 
}
\description{
Given a noisy or not planar curve as a set of discrete  \eqn{\left\{(x_i,y_i),i=1,2,\ldots n\right\}}{{(x_i,y_i),i=1,2,...n}} points we use Integration Extreme Finding Estimator (IEFE) algorithm as it is described at [1] in ordfer to find the extreme point of it.
}
\usage{
findextreme(x, y, parallel = FALSE, silent = TRUE, tryfast = FALSE)
}
\arguments{
  \item{x}{
A numeric vector for the independent variable without missing values
}
  \item{y}{
A numeric vector for the dependent variable without missing values
}
  \item{parallel}{
Logical input, if TRUE then parallel processing will be used (default=FALSE)
}
  \item{silent}{
Logical input, if TRUE then no details will be printed out during code execution (default=TRUE)
}
  \item{tryfast}{
Logical input, if TRUE then instead 'BEDE' will be used from IEFE algorithm instead of BESE (default=FALSE)
}
}
\details{
The parallel=TRUE otpion must be used if length(x)>20000. The tryfast=TRUE can be used for big data sets, but BEDE is not so accuracy as BESE, so use it with caution.
}
\value{ A named vector with next components is returned:
\enumerate{
  \item{x1 }{the left endpoint of the final interval of BESE or BEDE iterations }
  \item{x2 }{the right endpoint of the final interval of BESE or BEDE iterations }
  \item{chi }{the estimation of extreme as x-abscissa}
  \item{chi }{the estimation of extreme as y-abscissa taken from the interpolation polynomial of 2nd degree for 
  the data points (x1,y1), (x2,y2), (chi,ychi)}}
}
\references{
[1]Demetris T. Christopoulos (2019). New methods for computing extremes and roots of a planar curve: introducing Noisy Numerical Analysis (2019). \emph{ResearchGate}. \url{http://dx.doi.org/10.13140/RG.2.2.17158.32324} \cr
}
\author{
Demetris T. Christopoulos
}
\note{
The 'yvalue' at output vector is an interpolation approxiamtion for the y-value of unknown function at its extreme point 'chi' and does not mean that it will be certainly accurate. Thta is the truth if underlying function can be well approximated by low order polynomials.
}
\seealso{
\code{\link{symextreme}}, \code{\link{findmaxtulip}}, \code{\link{findmaxbell}} 
}
\examples{
## Legendre polynomial 5th order
## True extreme point p=0.2852315165, y=0.3466277
f=function(x){(63/8)*x^5-(35/4)*x^3+(15/8)*x} 
x=seq(0,0.7,0.001);y=f(x)
plot(x,y,pch=19,cex=0.5)
a=findextreme(x,y)
a
##        x1        x2       chi    yvalue 
## 0.2840000 0.2860000 0.2850000 0.3466274 
sol=a['chi']
abline(h=0)
abline(v=sol)
abline(v=a[1:2],lty=2)
abline(h=f(sol),lty=2)
points(sol,f(sol),pch=17,cex=2)
#
## The same function with noise from U(-0.05,0.05)
set.seed(2019-07-26);r=0.05;y=f(x)+runif(length(x),-r,r)
plot(x,y,pch=19,cex=0.5)
a=findextreme(x,y)
a
##        x1        x2       chi    yvalue 
## 0.2890000 0.2910000 0.2900000 0.3895484 
sol=a['chi']
abline(h=0)
abline(v=sol)
abline(v=a[1:2],lty=2)
abline(h=f(sol),lty=2)
points(sol,f(sol),pch=17,cex=2)
#
}
\keyword{IEFE}
\keyword{BESE}
\keyword{BEDE}
\keyword{maximum}
\keyword{minimum}
