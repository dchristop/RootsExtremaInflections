\name{scan_noisy_curve}
\alias{scan_noisy_curve}
\title{
Finds roots, extrema and inflections for a planar noisy curve
}
\description{
Given a noisy planar curve as a set of discrete \eqn{\left\{(x_i,y_i),i=1,2,\ldots n\right\}}{{(x_i,y_i),i=1,2,...n}} points we use Integration Root Finding Estimator (IRFE), Integration Extreme Finding Estimator (IEFE) and BESE in order to find all roots and the enclosed between them extrema and inflction points. Estimators are defined and described at [1] and [2], [3].
}
\usage{
scan_noisy_curve(x, y, noise = NULL, rootsoptim = TRUE, findextremes = TRUE, 
  findinflections = TRUE, silent = FALSE, plots = TRUE)
}
\arguments{
  \item{x}{
A numeric vector for the independent variable without missing values
}
  \item{y}{
A numeric vector for the dependent variable without missing values
}
  \item{noise}{
Either NULL (default value) or TRUE/FALSE input for suggesting that curve is indeed a noisy one
}
  \item{rootsoptim}{
Logical input, if TRUE then IRFE algorithm will also be used for roots (default=TRUE)
}
  \item{findextremes}{
Logical input, if TRUE then find all existing extrema between the found roots (default=TRUE)
}
  \item{findinflections}{
Logical input, if TRUE then find all existing inflection points between the found roots (default=TRUE)
}
  \item{silent}{
Logical input, if TRUE then no details will be printed out during code execution (default=TRUE)
}
  \item{plots}{
Logical input, if TRUE then a plot with all roots, extrema and inflctions will be created (default=TRUE)
}
}
\details{
The function first uses \code{\link{scan_curve}} in norder to perform a first study. Then it uses an extesnsive if-else environment in order to cover all possible cases without error breaks.
}
\value{A list with next memebers is returned:
\enumerate{
 \item{study }{a data frame with all intervals that curve can be divided and next columns:
    \itemize{
      \item{j }{the index of x vector for the left endpoint of current interval}
      \item{dj }{the number of x points of current interval}
      \item{interval }{logical, if TRUE then it is a unique convexity interval}
      \item{i1 }{the index of x vector for the left endpoint of current interval}
      \item{i2 }{the index of x vector for the right endpoint of current interval}
      \item{interval }{logical, if TRUE then a root exists inside current interval}
    }
  }
  \item{roots_average }{a data frame with all averaged roots and next columns:
    \itemize{
      \item{x1 }{the left endpoint of the averaged interval }
      \item{x2 }{the right endpoint of the averaged interval }
      \item{chi }{the estimation of root as x-abscissa}
      \item{chi }{the estimation of root as y-abscissa taken from the interpolation polynomial of 2nd degree for 
      a properly collected triad of (x[i],y[i]) points}
    }
  }
  \item{roots_optim }{a data frame with all roots from IRFE and next columns:
    \itemize{
      \item{x1 }{the left endpoint of the final interval of BESE iterations }
      \item{x2 }{the right endpoint of the final interval of BESE iterations }
      \item{chi }{the estimation of root as x-abscissa}
      \item{chi }{the estimation of root as y-abscissa taken from the interpolation polynomial of 2nd degree for 
      the data points (x1,y1), (x2,y2), (chi,ychi)}
    }
  }
   \item{extremes }{a data frame with all extremes between roots and next columns:
    \itemize{
      \item{x1 }{the left endpoint of the final interval of BESE iterations }
      \item{x2 }{the right endpoint of the final interval of BESE iterations }
      \item{chi }{the estimation of extreme as x-abscissa}
      \item{chi }{the estimation of extreme as y-abscissa taken from the interpolation polynomial of 2nd degree for 
      the data points (x1,y1), (x2,y2), (chi,ychi)}
      }
    }
  \item{inflections }{a data frame with all inflections between roots and next columns:
    \itemize{
      \item{x1 }{the left endpoint of the final interval of BESE iterations }
      \item{x2 }{the right endpoint of the final interval of BESE iterations }
      \item{chi }{the estimation of inflection as x-abscissa}
      \item{chi }{the estimation of inflection as y-abscissa taken from the interpolation polynomial of 2nd degree for 
      the data points (x1,y1), (x2,y2), (chi,ychi)}
      }
    }
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
If curve is not a noisy one, then \code{\link{scan_curve}} is used.
}
\seealso{
\code{\link{scan_curve}}, \code{\link{classify_curve}}
}
\examples{
#
## Legendre polynomial 5th order
f=function(x){(63/8)*x^5-(35/4)*x^3+(15/8)*x} 
x=seq(-1,1,0.001)
set.seed(2019-07-26);r=0.05;y=f(x)+runif(length(x),-r,r)
plot(x,y,pch=19,cex=0.5)
rn=scan_noisy_curve(x,y)
rn
## $study
##       j  dj interval   i1   i2  root
## 3    97 351     TRUE   97  448 FALSE
## 18  477 502     TRUE  477  979 FALSE
## 39 1021 505     TRUE 1021 1526 FALSE
## 54 1558 343     TRUE 1558 1901 FALSE
## 
## $roots_average
##       x1     x2     chi       yvalue
## 1 -0.906 -0.904 -0.9050 -0.002342389
## 2 -0.553 -0.524 -0.5385  0.005003069
## 3 -0.022  0.020 -0.0010  0.003260937
## 4  0.525  0.557  0.5410 -0.007956680
## 5  0.900  0.911  0.9055 -0.008015683
## 
## $roots_optim
##       x1     x2     chi       yvalue
## 1 -0.909 -0.901 -0.9050 -0.023334404
## 2 -0.531 -0.527 -0.5290  0.029256059
## 3  0.001  0.003  0.0020  0.001990572
## 4  0.530  0.565  0.5475  0.019616283
## 5  0.909  0.912  0.9105  0.009288338
## 
## $extremes
##          x1     x2     chi     yvalue
## [1,] -0.773 -0.766 -0.7695  0.4102010
## [2,] -0.280 -0.274 -0.2770 -0.3804006
## [3,]  0.308  0.316  0.3120  0.3372764
## [4,]  0.741  0.744  0.7425 -0.4414494
## 
## $inflections
##          x1     x2     chi       yvalue
## [1,] -0.772 -0.275 -0.5235 -0.076483193
## [2,] -0.275  0.281  0.0030 -0.007558037
## [3,]  0.301  0.776  0.5385  0.018958334
#
}
\keyword{noise}
\keyword{IEFE}
\keyword{IRFE}
\keyword{BEDE}
\keyword{root}
\keyword{maximum}
\keyword{minimum}
\keyword{inflection}
