\name{classify_curve}
\alias{classify_curve}
\title{
Classifies a planar curve according to its convexity, symmetry and extreme type
}
\description{
Given a planar curve we want to know its convexity ranges, the index for use in ESE or EDE methods, the existence of a  maximum or minimum and the classification as tulip or bell.
}
\usage{
classify_curve(x, y)
}
\arguments{
  \item{x}{
A numeric vector for the independent variable without missing values
}
  \item{y}{
A numeric vector for the dependent variable without missing values
}
}
\details{
A first check for the existence of infiniteies is performed. If infinities exist, then all outputs are set to 'NA', otherwise main code works.
}
\value{
A list with next members:
\enumerate{
  \item{ctype }{the convexity type of curve: convex, concave, convex/concave or concave/convex}
  \item{index }{the relevant index for ESE or EDE usage: 0,1,0,1 respectively to 'ctype' previously presented values}
  \item{asymmetry }{a classification of asymmetry type, if exists}
  \item{totalconvexity }{the overall dominant convexity type of curve if we omit the fact probably has ranges of different types of convexity}
  \item{ismax }{logical value, TRUE if a maximum seems to exits, FALSE if a minimum seems to be the case, 'NA' otherwise}
  \item{shapetype }{the shape of extreme as tulip, bell or 'NA'}
}
}
\references{
[1]Demetris T. Christopoulos (2019). New methods for computing extremes and roots of a planar curve: introducing Noisy Numerical Analysis (2019). \emph{ResearchGate}. \url{http://dx.doi.org/10.13140/RG.2.2.17158.32324} \cr

[2]Demetris T. Christopoulos (2014). Developing methods for identifying the inflection point of a convex/concave curve. arXiv:1206.5478v2 [math.NA]. \url{https://arxiv.org/pdf/1206.5478v2.pdf} \cr

[3]Demetris T. Christopoulos (2016). On the efficient identification of an inflection point.International Journal of Mathematics and Scientific Computing, (ISSN: 2231-5330), vol. 6(1). \url{https://veltech.edu.in/wp-content/uploads/2016/04/Paper-04-2016.pdf} 
}
\author{
Demetris T. Christopoulos
}
\note{
Results of current function have an approximation type, since not all kind of curves can be classified by a given procedure.
Caution has been taken in order to be able to infer for the very basic attributes.
}
\seealso{
\code{\link{symextreme}}, \code{\link{scan_curve}}, \code{\link{scan_noisy_curve}}
}
\examples{
## Lets create a convex/concave curve and classify it:
f=function(x){5+5*tanh(x-5)}
x=seq(0,8,0.01)
y=f(x)
plot(x,y,pch=19,cex=0.2)
cc=classify_curve(x,y)
cc$ctype
## [1] "convex_concave"
cc$index
## [1] 0
## Use 'index':
ede(x,y,cc$index)
##      j1  j2 chi
## EDE 369 633   5
## Lets create an 'almost convex' curve and see it:
f=function(x){-x^3+5*x-6}
x=seq(-2.5,1.5,0.01)
y=f(x)
plot(x,y,pch=19,cex=0.2)
cc=classify_curve(x,y)
cc$totalconvexity
## [1] "convex"
## Check for existence of a maximum for a noisy curve:
f=function(x){100-(x-5)^2}
x=seq(0,12,by=0.01);
r=1.5;y=f(x)+runif(length(x),-r,r)
plot(x,y,pch=19,cex=0.2)
cc=classify_curve(x,y)
cc$ismax
## [1] TRUE
}
\keyword{curve}
\keyword{convexity}
\keyword{asymmetry}
\keyword{tulip}
\keyword{bell}
\keyword{extreme}
