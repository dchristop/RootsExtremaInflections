findmaxbell=function(x,y,concave=TRUE)
{
  # Function to compute bell type symmetrical extremes
  if(!concave){y=-y}
  sll=rep(0,length(x)-1);
  srr=sll;
  for (i in (2:length(x))){
    a=findipl(x,y,i)
    sll[i-1]=a[3];   
    srr[i-1]=a[4];
  }
  jl=which.min(sll)+1;
  jr=which.min(srr)+1;
  xs=0.5*(x[jl]+x[jr])
  out=c(jl,jr,xs)
  names(out)=c("j1","j2","chi")
  return(out)
}