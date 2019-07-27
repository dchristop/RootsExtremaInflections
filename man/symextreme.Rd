\name{symextreme}
\alias{symextreme}
\title{
A wrapper for using either Tulip Extreme Finding Estimator (TEFE) or Bell Extreme Finding Estimator (BEFE) algorithm for a symmetric extreme point
}
\description{
Function uses \code{\link{classify_curve}} and then either \code{\link{findmaxbell}} or  \code{\link{findmaxtulip}} to proceed. See [1] for definitions and details of algorithms.
}
\usage{
symextreme(x, y, concave = NULL, type = NULL)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{type}{
A character string inpute denoting the shape of the curve, either 'bell' or 'tulip' (default=NULL)
}
}
\details{
This function is useful if we know that our curve has a symmetry around its extreme point but we cannot directly infer for the relevant shape.
}
\value{A list with next memebers is returned:
\enumerate{
 \item{maximum }{logical, if TRUE then curve has a maximum}
 \item{minimum }{logical, if TRUE then curve has a minimum}
 \item{results }{a named vector with components:
   \itemize{
      \item{j1 }{the index of x-left}
      \item{j2 }{the index of x-right}
      \item{chi }{the estimation of extreme as x-abscissa}
      }
      }
}
}
\references{
[1]Demetris T. Christopoulos (2019). New methods for computing extremes and roots of a planar curve: introducing Noisy Numerical Analysis (2019). \emph{ResearchGate}. \url{http://dx.doi.org/10.13140/RG.2.2.17158.32324} \cr
}
\author{
Demetris T. Christopoulos
}
\note{
You can Use the 'type' input if you are sure for the shape of the curve.
}
\seealso{
\code{\link{classify_curve}},  \code{\link{findmaxtulip}}, \code{\link{findmaxbell}}, \code{\link{findextreme}}
}
\examples{
#
## Bell curve
f=function(x){1/(1+x^2)}
x=seq(-2,2.0,by=0.01);y=f(x)
plot(x,y,pch=19,cex=0.5)
a=symextreme(x,y)
a
## $maximum
## [1] TRUE
## 
## $minimum
## [1] FALSE
## 
## $results
##           j1           j2          chi 
## 1.770000e+02 2.250000e+02 1.110223e-16 
abline(v=a$results['chi'])
#
## Tulip curve
f=function(x){100-(x-5)^2}
x=seq(0,12,by=0.01);y=f(x)
plot(x,y,pch=19,cex=0.5)
a=symextreme(x,y)
a
## $maximum
## [1] TRUE
## 
## $minimum
## [1] FALSE
## 
## $results
##     j1   j2  chi 
## 1 1001    5 
abline(v=a$results['chi'])
#
}
\keyword{bell}
\keyword{tulip}
\keyword{TEFE}
\keyword{BEFE}
\keyword{maximum}
\keyword{minimum}
\keyword{symmetry}