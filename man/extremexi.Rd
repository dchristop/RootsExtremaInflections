\name{extremexi}
\alias{extremexi}
\title{
Function to Find the Extreme Point of a Planar Curve
}
\description{
It takes as input the x, y numeric vectors, the indices for the range to be searched plus some other 
options and finds the extreme point for that interval, while it plots data, Taylor polynomial and the computed \eqn{|a_1|}{|a_1|} coefficients.
}
\usage{
extremexi(x, y, i1, i2, nt, alpha = 5, xlb = "x", ylb = "y", xnd = 3,  ynd = 3,
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
The point \eqn{x_i}{x_i} which makes the relevant \eqn{|a_1|}{|a_1|} minimum is the estimation for the function's extreme point at the interval \eqn{[x_{i1},x_{i2}]}{[x_{i1},x_{i2}]}.
}
\value{
It returns an environment with two components:
  \item{an}{a matrix with 3 columns: lower, upper bound of confidence interval and middle value for each coefficient an}
  \item{fextr}{a list with 2 members: the position i and the value of the estimated extreme point \eqn{\rho=x_{i}}{\rho=x_i}}
}
\references{
Demetris T. Christopoulos (2014). Roots, extrema and inflection points by using a proper Taylor regression procedure. \emph{SSRN}. \url{https://dx.doi.org/10.2139/ssrn.2521403}
}
\author{
Demetris T. Christopoulos
}
\section{Warnings }{When you are using RStudio it is necessary to leave enough space for the plot window in order for the plots to appear normally.
The data should come from a function at least \eqn{C^{(1)}}{C^(1)} in order to be able to find an extreme point, if exists.
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
#Find extreme point, plot results, print Taylor coefficients and rho estimation:
#
c<-extremexi(x,y,1,length(x),5,5,plots=TRUE);c$an;c$fextr;
#
#Find multiple extrema.
#Let's create some data:
#
f=function(x){3*cos(x-5)};xa=0.;xb=9;
set.seed(12345);x=sort(runif(101,xa,xb));r=0.1;y=f(x)+2*r*(runif(length(x))-0.5);plot(x,y)
#
#The first extreme point is
c1<-extremexi(x,y,1,40,5,5,plots=TRUE);c1$an;c1$fextr;
#          2.5 %      97.5 %           an
# a0 -3.02708631 -2.94592364 -2.986504975
# a1  0.07660314  0.24706531  0.161834227
# a2  1.42127770  1.58580632  1.503542012
# a3 -0.09037154  0.10377241  0.006700434
# a4 -0.14788899 -0.08719428 -0.117541632
# a5 -0.03822416  0.01425066 -0.011986748
# [1] 22.000000  1.917229
#Compare it with the actual rho_1=1.858407346
#
#The second extreme point is
c2<-extremexi(x,y,50,80,5,5,plots=TRUE);c2$an;c2$fextr;
#          2.5 %       97.5 %         an
# a0  2.89779980  3.064703163  2.9812515
# a1  0.27288720  0.541496278  0.4071917
# a2 -1.81454401 -0.677932480 -1.2462382
# a3 -1.76290384  0.216201349 -0.7733512
# a4  0.02548354  1.269671304  0.6475774
# a5 -0.25156866  0.007565154 -0.1220018
# [1] 7.000000 4.896521
#You have to compare it with the actual value of rho_2=5.0
#
#Finally the third extreme point is
c3<-extremexi(x,y,80,length(x),5,5,plots=TRUE);c3$an;c3$fextr;
#         2.5 %     97.5 %         an
# a0 -3.0637461 -2.9218614 -2.9928037
# a1 -0.2381605  0.2615635  0.0117015
# a2  0.7860259  2.0105383  1.3982821
# a3 -1.4187417  0.7472155 -0.3357631
# a4 -0.7943208  1.0876143  0.1466468
# a5 -0.6677733  1.7628833  0.5475550
# [1] 11.000000  8.137392
#You have to compare it with the actual value of rho_3=8.141592654
}
\keyword{ extremexi }