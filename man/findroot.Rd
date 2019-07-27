\name{findroot}
\alias{findroot}
\title{
Find the root for a  planar curve byn using Integration Root Finding Estimator (IEFE) algorithm
}
\description{
Given a noisy or not planar curve as a set of discrete \eqn{\left\{(x_i,y_i),i=1,2,\ldots n\right\}}{{(x_i,y_i),i=1,2,...n}} points we use Integration Root Finding Estimator (IRFE) algorithm as it is described at [1] in ordfer to find the root of it.
}
\usage{
findroot(x, y, parallel = FALSE, silent = TRUE, tryfast = FALSE)
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
[1]Demetris T. Christopoulos (2019). New methods for computing extremes and roots of a planar curve: introducing Noisy Numerical Analysis (2019). \emph{ResearchGate}. \url{http://dx.doi.org/10.13140/RG.2.2.17158.32324}
}
\author{
Demetris T. Christopoulos
}
\note{
The 'yvalue' at output vector is an interpolation approxiamtion for the y-value of unknown function at its extreme point 'chi' and does not mean that it will be certainly accurate. Thta is the truth if underlying function can be well approximated by low order polynomials.
}
\seealso{
\code{\link{scan_curve}}, \code{\link{scan_noisy_curve}}
}
\examples{
#
## Legendre polynomial 5th order
f=function(x){(63/8)*x^5-(35/4)*x^3+(15/8)*x} 
x=seq(0.2,0.8,0.001);y=f(x);ya=abs(y)
plot(x,y,pch=19,cex=0.5,ylim=c(min(y),max(ya)))
abline(h=0);
lines(x,ya,lwd=4,col='blue')
rt=findroot(x,y)
rt
##           x1            x2           chi        yvalue 
## 5.370000e-01  5.400000e-01  5.385000e-01 -7.442574e-05 
abline(v=rt['chi'])
abline(v=rt[1:2],lty=2);abline(h=rt['yvalue'],lty=2)
points(rt[3],rt[4],pch=17,col='blue',cex=2)
#
## Same curve but with noise from U(-0.5,0.5)
#
set.seed(2019-07-24);r=0.05;y=f(x)+runif(length(x),-r,r)
ya=abs(y)
plot(x,y,pch=19,cex=0.5,ylim=c(min(y),max(ya)))
abline(h=0)
points(x,ya,pch=19,cex=0.5,col='blue')
rt=findroot(x,y)
rt
##         x1          x2         chi      yvalue 
## 0.53400000  0.53700000  0.53550000 -0.01762159 
abline(v=rt['chi'])
abline(v=rt[1:2],lty=2);abline(h=rt['yvalue'],lty=2)
points(rt[3],rt[4],pch=17,col='blue',cex=2)
#
}
\keyword{root}
\keyword{IRFE}