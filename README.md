Overview
--------

RootsExtremaInflection is a package that finds roots, extrema and
inflection points of a planar curve which is given as a data frame of
discrete (xi,yi) points.

Basic functions are:

-   `rootxi()` Taylor Regression Estimator for roots
-   `extremexi()` Taylor Regression Estimator for extrema
-   `inflexi()` Taylor Regression Estimator for inflections
-   `classify_curve()` curve classification by convexity and shape type
-   `findmaxbell()` Bell Extreme Finding Estimator for symmetric extrema
-   `findmaxtulip()` Tulip Extreme Finding Estimator for symmetric
    extrema
-   `findextreme()` Integration Extreme Finding Estimator for all types
    of extrema
-   `findroot()` Integration Root Finding Estimator for roots
-   `scan_curve()` scans a not noisy curve and finds all roots and
    extrema, inflections between them
-   `scan_noisy_curve()` scans a noisy curve and finds all roots and
    extrema, inflections between them

Installation
------------

``` r
# Install with dependencies:
install.packages('RootsExtremaInflection',dependencies=TRUE)
```

Usage
-----

``` r
library(RootsExtremaInflection)

Load data:
data(xydat)
#
#Extract x and y variables:
x=xydat$x;y=xydat$y
#
#Find root, plot results, print Taylor coefficients and rho estimation:
b<-rootxi(x,y,1,length(x),5,5,plots=TRUE);b$an;b$froot;
#
#Find extreme, plot results, print Taylor coefficients and rho estimation:
c<-extremexi(x,y,1,length(x),5,5,plots=TRUE);c$an;c$fextr;
#
#Find inflection point, plot results, print Taylor coefficients and rho estimation:
d<-inflexi(x,y,1,length(x),5,5,plots=TRUE);d$an;d$finfl;
# Create a relative big data set...
f=function(x){3*cos(x-5)};xa=0.;xb=9;
set.seed(12345);x=sort(runif(5001,xa,xb));r=0.1;y=f(x)+2*r*(runif(length(x))-0.5);
#
#Find root, plot results, print Taylor coefficients and rho estimation in parallel:
#b1<-rootxi(x,y,1,round(length(x)/2),5,5,plots=TRUE,doparallel = TRUE);b1$an;b1$froot;
# Available workers are 12 
# Time difference of 5.838743 secs
#           2.5 %       97.5 %           an
# a0 -0.006960052  0.004414505 -0.001272774
# a1 -2.982715739 -2.933308292 -2.958012016
# a2 -0.308844145 -0.213011162 -0.260927654
# a3  0.806555336  0.874000586  0.840277961
# a4 -0.180720951 -0.161344935 -0.171032943
# a5  0.007140500  0.009083859  0.008112180
# [1] 177.0000000   0.2924279
# Compare with exact root = 0.2876110196
#Find extreme, plot results, print Taylor coefficients and rho estimation in parallel:
#c1<-extremexi(x,y,1,round(length(x)/2),5,5,plots=TRUE,doparallel = TRUE);c1$an;c1$fextr;
# Available workers are 12 
# Time difference of 5.822514 secs
#            2.5 %       97.5 %           an
# a0 -3.0032740050 -2.994123850 -2.998698927
# a1 -0.0006883998  0.012218393  0.005764997
# a2  1.4745326519  1.489836668  1.482184660
# a3 -0.0340626683 -0.025094859 -0.029578763
# a4 -0.1100798736 -0.105430525 -0.107755199
# a5  0.0071405003  0.009083859  0.008112180
# [1] 1022.000000    1.852496
# Compare with exact extreme = 1.858407346
#Find inflection point, plot results, print Taylor coefficients and rho estimation in parallel:
#d1<-inflexi(x,y,1090,2785,5,5,plots=TRUE,doparallel = TRUE);d1$an;d1$finfl;
# Available workers are 12 
# Time difference of 4.343851 secs
#           2.5 %       97.5 %            an
# a0 -0.008238016  0.002091071 -0.0030734725
# a1  2.995813560  3.023198534  3.0095060468
# a2 -0.014591465  0.015326175  0.0003673549
# a3 -0.531029710 -0.484131902 -0.5075808056
# a4 -0.008253975  0.007556465 -0.0003487551
# a5  0.016126428  0.034688019  0.0254072236
# [1] 800.000000   3.427705
# Compare with exact inflection = 3.429203673
# Or execute rootexinf() and find a set of them at once and in same time:
#a<-rootexinf(x,y,100,round(length(x)/2),5,plots = TRUE,doparallel = TRUE);
#a$an0;a$an1;a$an2;a$frexinf;
# Available workers are 12 
# Time difference of 5.565372 secs
#           2.5 %      97.5 %           an0
# a0 -0.008244278  0.00836885  6.228596e-05
# a1 -2.927764078 -2.84035634 -2.884060e+00
# a2 -0.447136449 -0.30473094 -3.759337e-01
# a3  0.857290490  0.94794071  9.026156e-01
# a4 -0.198104383 -0.17360676 -1.858556e-01
# a5  0.008239609  0.01059792  9.418764e-03
#           2.5 %      97.5 %          an1
# a0 -3.005668018 -2.99623116 -3.000949590
# a1 -0.003173501  0.00991921  0.003372854
# a2  1.482600580  1.50077450  1.491687542
# a3 -0.034503271 -0.02551597 -0.030009618
# a4 -0.115396537 -0.10894117 -0.112168855
# a5  0.008239609  0.01059792  0.009418764
#           2.5 %       97.5 %          an2
# a0  0.083429390  0.092578772  0.088004081
# a1  3.007115452  3.027343849  3.017229650
# a2 -0.009867779  0.006590042 -0.001638868
# a3 -0.517993955 -0.497886933 -0.507940444
# a4 -0.043096158 -0.029788902 -0.036442530
# a5  0.008239609  0.010597918  0.009418764
#            index     value
# root          74 0.2878164
# extreme      923 1.8524956
# inflection  1803 3.4604842
#
## Next examples are for the
## Legendre polynomial of 5th order:
#
f=function(x){(63/8)*x^5-(35/4)*x^3+(15/8)*x} 
#
### findextreme()
#
## True extreme point p=0.2852315165, y=0.3466277
x=seq(0,0.7,0.001);y=f(x)
plot(x,y,pch=19,cex=0.5)
a=findextreme(x,y)
a
##        x1        x2       chi    yvalue 
## 0.2840000 0.2860000 0.2850000 0.3466274 
sol=a['chi']
abline(h=0)
abline(v=sol)
abline(v=a[1:2],lty=2)
abline(h=f(sol),lty=2)
points(sol,f(sol),pch=17,cex=2)
#
## The same function with noise from U(-0.05,0.05)
set.seed(2019-07-26);r=0.05;y=f(x)+runif(length(x),-r,r)
plot(x,y,pch=19,cex=0.5)
a=findextreme(x,y)
a
##        x1        x2       chi    yvalue 
## 0.2890000 0.2910000 0.2900000 0.3895484 
sol=a['chi']
abline(h=0)
abline(v=sol)
abline(v=a[1:2],lty=2)
abline(h=f(sol),lty=2)
points(sol,f(sol),pch=17,cex=2)
#
### findroot()
#
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
### scan_curve()
#
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
### scan_noisy_curve()
#
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

```

Why should I use RootsExtremaInflection package in R?
-----------------------------------------------------

-   Because it can give you a reliable estimation by using Taylor
    Regression, if you prefer non parametric statistical methods
-   Because it can find the desired point(s) by using geometric and
    analytic methods, without any kind of model adoption
-   Because it is the core of a new field called ‘Noisy Numerical
    Analysis’ which combines Geometry, Calculus, Numerical Analysis and
    Statistics in order to solve known problems of Numerical Analysis

Contact
-------

Please send comments, suggestions or bug breports to
dchristop$econ.uoa.gr
