findmaxtulip=function(x,y,concave=TRUE)
{
  # Function to compute tulip type symmetrical extremes
  if(!concave){y=-y}
  n=length(x)
  xm=c(x[1],x[n])
  ym=c(y[1],y[n])
  ii=which.max(ym);#print(ii);
  if(ii==1){
    i1=1
    y1=y[1]
    # It is the area above the horizontal line y=y1 and the curve
    # s=sapply(2:(n-1), function(j,x,y){findareatpl(x,y,j)},x,y)
    # sli=vapply(1:(n-1),function(ii,x,y,y1){0.5*(x[ii + 1] - x[ii])*(y[ii] + y[ii + 1] - 2*y1)},FUN.VALUE=numeric(1),x,y,y1)
    dxl = diff(x[1:n], 1, 1)
    fl = y[1:n] - y1
    sl = 0.5 * (fl[1:(n - 1)] + fl[2:n])*dxl
    s=cumsum(sl)
    jm=which.max(s)
    im=jm+1
    i2=im 
  }
  else
  {
    i2=n
    yn=y[n]
    # It is the area above the horizontal line y=yn and the curve
    # s=sapply((n-1):2, function(j,x,y){findareatpr(x,y,j)},x,y)
    # sri=vapply((n-1):1,function(ii,x,y,yn){0.5*(x[ii + 1] - x[ii])*(y[ii] + y[ii + 1] - 2*yn)},FUN.VALUE=numeric(1),x,y,yn)
    dxr = diff(x[1:n], 1, 1)
    fr = y[1:n] - y[n]
    sr = rev(0.5 * (fr[1:(n - 1)] + fr[2:n])*dxr)
    s=cumsum(sr)
    jm=which.max(s)
    im=n+1-jm #--> n-jm
    i1=im
  }
  out=c(i1,i2,(x[i1]+x[i2])/2)
  names(out)=c("j1","j2","chi")
  return(out)
}  
