\name{findmaxtulip}
\alias{findmaxtulip}
\title{
Implementation of Tulip Extreme Finding Estimator (TEFE) algorithm for a planar curve
}
\description{
For a curve that can be classified as 'tulip' a fast computation of its maximum is performed by applying Tulip Extreme Finding Estimator (TEFE) algorithm of [1].
}
\usage{
findmaxtulip(x, y, concave = TRUE)
}
\arguments{
  \item{x}{
A numeric vector for the independent variable without missing values
}
  \item{y}{
A numeric vector for the dependent variable without missing values
}
  \item{concave}{
Logical input, if TRUE then curve is supposed to have a maximum (default=TRUE)
}
}
\details{
If we want to compute minimum we just set concave=FALSE and proceed.
}
\value{
A named vector with next components is returned:
\enumerate{
  \item{j1 }{the index of x vextor that is the left endpoint of final symmetrization interval}
  \item{j1 }{the index of x vextor that is the right endpoint of final symmetrization interval}
  \item{chi }{the estimation of extreme as x-abscissa}
}
}
\references{
[1]Demetris T. Christopoulos (2019). New methods for computing extremes and roots of a planar curve: introducing Noisy Numerical Analysis (2019). \emph{ResearchGate}. \url{http://dx.doi.org/10.13140/RG.2.2.17158.32324}
}
\author{
Demetris T. Christopoulos
}
\note{
Please use function \code{\link{classify_curve}} if you have not visual inspection in order to find the extreme type. Do not use that function if curve shape is not 'tulip', use either \code{\link{symextreme}} or \code{\link{findextreme}}.
}
\seealso{
\code{\link{classify_curve}}, \code{\link{symextreme}}, \code{\link{findmaxbell}}, \code{\link{findextreme}}
}
\examples{
#
f=function(x){100-(x-5)^2}
x=seq(0,12,by=0.01);y=f(x)
plot(x,y,pch=19,cex=0.5)
cc=classify_curve(x,y)
cc$shapetype
## 1] "tulip"
a<-findmaxtulip(x,y)
a
## j1   j2  chi 
##  1 1001    5 
abline(v=a['chi'])
abline(v=x[a[1:2]],lty=2);abline(h=y[a[1:2]],lty=2)
points(x[a[1]:a[2]],y[a[1]:a[2]],pch=19,cex=0.5,col='blue')
#
## Same curve with noise from U(-1.5,1.5)
set.seed(2019-07-26);r=1.5;y=f(x)+runif(length(x),-r,r)
plot(x,y,pch=19,cex=0.5)
cc=classify_curve(x,y)
cc$shapetype
## 1] "tulip"
plot(x,y,pch=19,cex=0.5)
a<-findmaxtulip(x,y)
a
##    j1       j2      chi 
## 1.000 1002.000    5.005 
abline(v=a['chi'])
abline(v=x[a[1:2]],lty=2);abline(h=y[a[1:2]],lty=2)
points(x[a[1]:a[2]],y[a[1]:a[2]],pch=19,cex=0.5,col='blue')
#
}
\keyword{tulip}
\keyword{TEFE}
\keyword{maximum}
\keyword{minimum}