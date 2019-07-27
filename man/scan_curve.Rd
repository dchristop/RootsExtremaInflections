\name{scan_curve}
\alias{scan_curve}
\title{
Finds roots, extrema and inflections for a planar not noisy curve
}
\description{
Given a not noisy planar curve as a set of discrete \eqn{\left\{(x_i,y_i),i=1,2,\ldots n\right\}}{{(x_i,y_i),i=1,2,...n}} points we use Integration Root Finding Estimator (IRFE), Integration Extreme Finding Estimator (IEFE) and BESE in order to find all roots and the enclosed between them extrema and inflction points. Estimators are defined and described at [1] and [2], [3].
}
\usage{
scan_curve(x, y, findroots = TRUE, findextremes = TRUE,
   findinflections = TRUE, silent = FALSE, plots = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector for the independent variable without missing values
}
  \item{y}{
A numeric vector for the dependent variable without missing values
}
  \item{findroots}{
Logical input, if TRUE then find all existing roots (default=TRUE)
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
The function first makes a study for the curve based on the interval classification done by absolute y-values. Then it uses an extensive if-else environment in order to cover all possible cases without error breaks.
}
\value{A list with next memebers is returned:
\enumerate{
 \item{study }{a data frame with all intervals that curve can be divided and next columns:
    \itemize{
      \item{x }{the root estimation taken from absolute value }
      \item{y }{the y-value for root estimation x}
      \item{y }{the absolute y-value for root estimation x}
      \item{zero }{if TRUE, then all values of transformed yi's are zero for current interval}
      \item{dif }{monotonicty index, if it is +1 (-1) then curve is increasing (decreasing)at root x }
      \item{ja }{the index of x vector for the left endpoint of current interval}
      \item{jb }{the index of x vector for the right endpoint of current interval}
      \item{root_monotonicity }{monotonicty characterization for current root x}
      \item{extremity }{extremity characterization for current interval}
      \item{sigmoidicity }{sigmoidicity characterization for current interval}
    }
  }
  \item{roots }{a data frame with all roots and next columns:
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
A heuristic is first used for determination of noise: if curve is classified as a noisy one, then only 'study' and a rough plot are available.
}
\seealso{
\code{\link{scan_noisy_curve}}, \code{\link{classify_curve}}
}
\examples{
#
## Legendre polynomial 5th order
f=function(x){(63/8)*x^5-(35/4)*x^3+(15/8)*x} 
x=seq(-1,1,0.001);y=f(x)
plot(x,y,pch=19,cex=0.5)
abline(h=0)
rall=scan_curve(x,y)
rall$study
rall$roots
##          x1     x2           chi        yvalue
## [1,] -0.907 -0.905 -9.060000e-01  1.234476e-03
## [2,] -0.540 -0.537 -5.385000e-01  7.447856e-05
## [3,] -0.001  0.001  5.551115e-17  1.040844e-16
## [4,]  0.537  0.540  5.385000e-01 -7.444324e-05
## [5,]  0.905  0.907  9.060000e-01 -1.234476e-03
rall$extremes
##          x1     x2    chi     yvalue
## [1,] -0.766 -0.764 -0.765  0.4196969
## [2,] -0.286 -0.284 -0.285 -0.3466274
## [3,]  0.284  0.286  0.285  0.3466274
## [4,]  0.764  0.766  0.765 -0.4196969
rall$inflections
##          x1     x2           chi        yvalue
## [1,] -0.579 -0.576 -5.775000e-01  9.659939e-02
## [2,] -0.001  0.001  5.551115e-17  1.040829e-16
## [3,]  0.576  0.579  5.775000e-01 -9.659935e-02
#
}
\keyword{IEFE}
\keyword{IRFE}
\keyword{BEDE}
\keyword{root}
\keyword{maximum}
\keyword{minimum}
\keyword{inflection}

