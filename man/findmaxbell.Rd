\name{findmaxbell}
\alias{findmaxbell}
\title{
Implementation of Bell Extreme Finding Estimator (BEFE) algorithm for a planar curve
}
\description{
For a curve that can be classified as 'bell' a fast computation of its maximum is performed by applying Bell Extreme Finding Estimator (BEFE) algorithm of [1].
}
\usage{
findmaxbell(x, y, concave = TRUE)
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
  \item{j1 }{the index of x-left}
  \item{j2 }{the index of x-right}
  \item{chi }{the estimation of extreme as x-abscissa}
}
}
\references{
[1]Demetris T. Christopoulos (2019). New methods for computing extremes and roots of a planar curve: introducing Noisy Numerical Analysis (2019). \emph{ResearchGate}. \url{http://dx.doi.org/10.13140/RG.2.2.17158.32324} \cr
}
\author{
Demetris T. Christopoulos
}
\note{
Please use function \code{\link{classify_curve}} if you have not visual inspection in order to find the extreme type. Do not use that function if curve shape is not 'bell', use either \code{\link{symextreme}} or \code{\link{findextreme}}.
}
\seealso{
\code{\link{classify_curve}}, \code{\link{symextreme}}, \code{\link{findmaxtulip}}, \code{\link{findextreme}}
}
\examples{
#
f=function(x){1/(1+x^2)}
x=seq(-2,2.0,by=0.01);y=f(x)
plot(x,y,pch=19,cex=0.5)
cc=classify_curve(x,y)
cc$shapetype
## 1] "bell"
a1<-findmaxbell(x,y)
a1
## j1           j2          chi 
## 1.770000e+02 2.250000e+02 1.110223e-16 
abline(v=a1['chi'])
abline(v=x[a1[1:2]],lty=2);abline(h=y[a1[1:2]],lty=2)
points(x[a1[1]:a1[2]],y[a1[1]:a1[2]],pch=19,cex=0.5,col='blue')
#
## Same curve with noise from U(-0.05,0.05)
set.seed(2019-07-26);r=0.05;y=f(x)+runif(length(x),-r,r)
plot(x,y,pch=19,cex=0.5)
cc=classify_curve(x,y)
cc$shapetype
## 1] "bell"
a1<-findmaxbell(x,y)
a1
##     j1     j2    chi 
## 169.00 229.00  -0.02 
abline(v=a1['chi'])
abline(v=x[a1[1:2]],lty=2);abline(h=y[a1[1:2]],lty=2)
points(x[a1[1]:a1[2]],y[a1[1]:a1[2]],pch=19,cex=0.5,col='blue')
#
}
\keyword{bell}
\keyword{BEFE}
\keyword{maximum}
\keyword{minimum}
