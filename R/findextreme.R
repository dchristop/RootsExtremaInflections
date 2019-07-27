findextreme=function(x,y,parallel=FALSE,silent=TRUE,tryfast=FALSE){
  # Function to compute general extremes
  # Function to compute integrals
  sumint=function (x, y, j)
  {
    dx = diff(x[1:j], 1, 1)
    fx = y[1:j] 
    strap = 0.5 * (fx[1:(j - 1)] + fx[2:j])
    sj = sum(dx * strap)
    c(j, x[j],sj)
  }
  # Polynomial interpolation of 2nd degree
  polint2=function(x,x1,x2,x3,y1,y2,y3){
    y1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3))+y2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3))+y3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2))
  }
  n=length(y)
  if(n<7){
    stop('Number of data points must be at least 7. Please increase the number of points.')
  }
  # Main algorithm
  # Compute integrals
  if(parallel){
    #
    tp1=Sys.time()
    runfind=function(j,x,y){
      dx = diff(x[1:j], 1, 1)
      fx = y[1:j] 
      strap = 0.5 * (fx[1:(j - 1)] + fx[2:j])
      return(sum(dx * strap))
    }
    environment(runfind) <- .GlobalEnv
    cl <- makeCluster(detectCores());registerDoParallel(cl);
    clusterEvalQ(cl=cl,list(library("inflection")))
    fs=parallel::parSapply(cl=cl,2:(n-1),runfind,x,y)
    stopCluster(cl)
    tp2=Sys.time();print(tp2-tp1) 
    t12=paste0("Time for computing surfaces was ",format(tp2-tp1,units="sec"))
    if(!silent){message(t12)}
    #
  }else{
    t1=Sys.time() 
    fs=c()
    fs=sapply(2:(n-1), function(j,x,y){sumint(x,y,j)[3]},x,y)
    t2=Sys.time()
    t12=paste0("Time for computing surfaces was ",format(t2-t1,units="sec"))
    if(!silent){message(t12)}
  }
  # Check sigmoidicity type:
  cc=check_curve(1:length(fs),fs)
  cc
  # Try BEDE for fast run, although not so accurate as BESE:
  if(tryfast){
    t1=Sys.time() 
    dd=bede(1:length(fs),fs,cc$index)
    t2=Sys.time()
    t12=paste0("Time for computing BEDE was ",format(t2-t1,units="sec"))
    if(!silent){message(t12)}
    # Check convergence for BEDE
    if(is.nan(dd$iters$EDE[dim(dd$iters)[1]])){
      message(paste0('BEDE was not able to find an inflection at [',x[1],' , ',x[n],'] \n Continue by using BESE...'))
    }else{
      solint=x[t(dd$iters[dim(dd$iters)[1],c("a","b")]+1)]  #Plus 1!
      sol=mean(solint)
      jj=c(unlist(dd$iters[dim(dd$iters)[1],c("a","b","EDE")]));jj
      j1=jj[1];j2=round(jj[3]);j3=jj[2];c(j1,j2,j3)
      yex=polint2(sol,x[j1],x[j2],x[j3],y[j1],y[j2],y[j3]);yex
      out=c(solint,sol,yex)
      names(out)=c("x1","x2","chi","yvalue")
      return(out)
    }
  }
  # Try BESE if BEDE failed
  t1=Sys.time() 
  bb=bese(1:length(fs),fs,cc$index,doparallel = parallel)
  t2=Sys.time()
  t12=paste0("Time for computing BESE was ",format(t2-t1,units="sec"))
  if(!silent){message(t12)}
  # Check existence of local extremes
  if(is.nan(bb$iters$ESE[dim(bb$iters)[1]])){
    message(paste0('It seems that no local extremes exist at [',x[1],' , ',x[n],']'))
    out=c("x1"=x[1],"x2"=x[n],"chi"=NA,"yvalue"=NA)
    return(out)
  }else{
    solint=x[t(bb$iters[dim(bb$iters)[1],c("a","b")]+1)]  #Plus 1!
    sol=mean(solint)
    jj=c(unlist(bb$iters[dim(bb$iters)[1],c("a","b","ESE")]));jj
    j1=jj[1];j2=round(jj[3]);j3=jj[2];c(j1,j2,j3)
    yex=polint2(sol,x[j1],x[j2],x[j3],y[j1],y[j2],y[j3]);yex
    out=c(solint,sol,yex)
    names(out)=c("x1","x2","chi","yvalue")
    return(out)
  }
}
