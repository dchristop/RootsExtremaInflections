\name{inflexi}
\alias{inflexi}
\title{
Function to Find the Inflection Point of a Planar Curve
}
\description{
It takes as input the x, y numeric vectors, the indices for the range to be searched plus some other 
options and finds the inflection point for that interval, while it plots data, Taylor polynomial and and the computed \eqn{|a_2|}{|a_2|} coefficients.
}
\usage{
inflexi(x, y, i1, i2, nt, alpha = 5, xlb = "x", ylb = "y", xnd = 3,  ynd = 3,
plots = TRUE, plotpdf = FALSE, doparallel=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector for the independent variable 
}
  \item{y}{
A numeric vector for the dependent variable 
}
  \item{i1}{
The first index for choosing a specific interval \eqn{[a,b]=[x_{i1},x_{i2}]}{[a,b]=[x_{i1},x_{i2}]}
}
  \item{i2}{
The second index for choosing a specific interval \eqn{[a,b]=[x_{i1},x_{i2}]}{[a,b]=[x_{i1},x_{i2}]}
}
  \item{nt}{
The degree of the Taylor polynomial that will be fitted to the data 
}
  \item{alpha}{
The level of statistical significance for the confidence intervals of coefficients \eqn{a_0, a_1,..., a_{nt-1}}{a_0, a_1,..., a_{nt-1}} (default value = 5)
}
  \item{xlb}{
A label for the x-variable (default value = "x")
}
  \item{ylb}{
A label for the y-variable (default value = "y")
}
  \item{xnd}{
The number of digits for plotting the x-axis (default value = 3)
}
  \item{ynd}{
The number of digits for plotting the y-axis (default value = 3)
}
  \item{plots}{
If plots=TRUE then a plot is created on default monitor (default value = TRUE) 
}
  \item{plotpdf}{
If plotpdf=TRUE then a pdf plot is created and stored on working directory (default value = FALSE) 
}
  \item{doparallel}{
If doparallel=TRUE then parallel computing is applied, based on the available workers of current machine (default value = FALSE)    
}
}
\details{
The point \eqn{x_i}{x_i} which makes the relevant \eqn{|a_2|}{|a_2|} minimum is the estimation for the function's inflection point at the interval \eqn{[x_{i1},x_{i2}]}{[x_{i1},x_{i2}]}.
}
\value{
It returns an environment with two components:
  \item{an}{a matrix with 3 columns: lower, upper bound of confidence interval and middle value for each coefficient an}
  \item{fextr}{a list with 2 members: the position i and the value of the estimated inflection point \eqn{\rho=x_{i}}{\rho=x_i}}
}
\references{
Demetris T. Christopoulos (2014). Roots, extrema and inflection points by using a proper Taylor regression procedure. \emph{SSRN}. \url{https://dx.doi.org/10.2139/ssrn.2521403}
}
\author{
Demetris T. Christopoulos
}
\section{Warnings }{When you are using RStudio it is necessary to leave enough space for the plot window in order for the plots to appear normally.
The data should come from a function at least \eqn{C^{(2)}}{C^(2)} in order to be able to find an inflection point, if exists.
} 
\examples{
#Load data:
#
data(xydat)
#
#Extract x and y variables:
#
x=xydat$x;y=xydat$y
#
#Find inflection point, plot results, print Taylor coefficients and rho estimation:
#
d<-inflexi(x,y,1,length(x),5,5,plots=TRUE);d$an;d$finfl;
#
#Find multiple inflection points.
#Let's create some data:
#
f=function(x){3*cos(x-5)};xa=0.;xb=9;
set.seed(12345);x=sort(runif(101,xa,xb));r=0.1;y=f(x)+2*r*(runif(length(x))-0.5);plot(x,y)
#
#The first inflection point is
d1<-inflexi(x,y,20,50,5,5,plots=TRUE);d1$an;d1$finfl;
#         2.5 %     97.5 %           an
# a0  0.1483905  0.2377617  0.193076089
# a1  2.9024852  3.0936024  2.998043835
# a2 -0.2053120  0.2220390  0.008363525
# a3 -0.5845597 -0.3426017 -0.463580702
# a4 -0.2431038  0.1136244 -0.064739689
# a5 -0.0893246  0.0687848 -0.010269897
# [1] 19.000000  3.493296
#Compare it with the actual rho_1=3.429203673
#
#The second inflection point is
# d2<-inflexi(x,y,50,length(x),5,5,plots=TRUE);d2$an;d2$finfl;
#           2.5 %       97.5 %            an
# a0 -0.000875677  0.057156356  0.0281403394
# a1 -3.058363342 -2.942026810 -3.0001950762
# a2 -0.056224101  0.044135857 -0.0060441222
# a3  0.433135897  0.528446241  0.4807910691
# a4 -0.011774733  0.012002414  0.0001138404
# a5 -0.026899286 -0.009520899 -0.0182100925
# [1] 23.000000  6.567948
#You have to compare it with the actual value of rho_2=6.570796327
}
\keyword{ extremexi }