\name{rootxi}
\alias{rootxi}
\title{
Function to Find the Root of a Planar Curve
}
\description{
It takes as input the x, y numeric vectors, the indices for the range to be searched plus some other 
options and finds the root for that interval, while it plots data, Taylor polynomial and and the computed \eqn{|a_0|}{|a_0|} coefficients.
}
\usage{
rootxi(x, y, i1, i2, nt, alpha = 5, xlb = "x", ylb = "y", xnd = 3,  ynd = 3,
plots = TRUE, plotpdf = FALSE, doparallel=FALSE)
}
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
The point \eqn{x_i}{x_i} which makes the relevant \eqn{|a_0|}{|a_0|} minimum is the estimation for the function's root at the interval \eqn{[x_{i1},x_{i2}]}{[x_{i1},x_{i2}]}.
}
\value{
It returns an environment with two components:
  \item{an}{a matrix with 3 columns: lower, upper bound of confidence interval and middle value for each coefficient a_n}
  \item{froot}{a list with 2 members: the position i and the value of the estimated root \eqn{\rho=x_{i}}{\rho=x_i}}
}
\references{
Demetris T. Christopoulos (2014). Roots, extrema and inflection points by using a proper Taylor regression procedure. \emph{SSRN}. \url{https://dx.doi.org/10.2139/ssrn.2521403}
}
\author{
Demetris T. Christopoulos
}
\section{Warnings }{When you are using RStudio it is necessary to leave enough space for the plot window in order for the plots to appear normally.
The data should come from a function at least \eqn{C^{(0)}}{C^(0)} in order to find the root, provided that such a root exists.
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
#Find root, plot results, print Taylor coefficients and rho estimation:
#
b<-rootxi(x,y,1,length(x),5,5,plots=TRUE);b$an;b$froot;
#
#Find multiple roots.
#Let's create some data:
#
f=function(x){3*cos(x-5)};xa=0.;xb=9;
set.seed(12345);x=sort(runif(101,xa,xb));r=0.1;y=f(x)+2*r*(runif(length(x))-0.5);plot(x,y)
#
#The first root is
#
b1<-rootxi(x,y,1,20,5,5,plots=TRUE);b1$an;b1$froot;
#          2.5 %      97.5 %          an
# a0 -0.09380972  0.03295954 -0.03042509
# a1 -3.63025679 -2.89908741 -3.26467210
# a2 -0.90435090  0.80658742 -0.04888174
# a3 -1.27911360  6.88168053  2.80128346
# a4 -8.77763032  2.51983279 -3.12889877
# a5 -1.10798564  3.38419904  1.13810670
# [1] 5.0000000 0.3108189
#Compare it with the actual rho_1=0.2876110196
#
#The second root is
#
b2<-rootxi(x,y,20,50,5,5,plots=TRUE);b2$an;b2$froot;
#         2.5 %     97.5 %           an
# a0  0.1483905  0.2377617  0.193076089
# a1  2.9024852  3.0936024  2.998043835
# a2 -0.2053120  0.2220390  0.008363525
# a3 -0.5845597 -0.3426017 -0.463580702
# a4 -0.2431038  0.1136244 -0.064739689
# a5 -0.0893246  0.0687848 -0.010269897
# [1] 19.000000  3.493296
#You have to compare it with the actual value of rho_2=3.429203673
#
#Finally the third root is
b3<-rootxi(x,y,50,90,5,5,plots=TRUE);b3$an;b3$froot;
#           2.5 %       97.5 %            an
# a0 -0.002269152  0.058784414  0.0282576308
# a1 -3.090980046 -2.938875341 -3.0149276930
# a2 -0.089893659  0.075094637 -0.0073995112
# a3  0.403040978  0.591836654  0.4974388159
# a4 -0.035442477  0.037165754  0.0008616385
# a5 -0.048414145  0.005815106 -0.0212995192
# [1] 23.000000  6.567948
#You have to compare it with the actual value of rho_3=6.570796327
}
\keyword{rootxi }