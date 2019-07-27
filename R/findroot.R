findroot=function(x,y,parallel=FALSE,silent=TRUE,tryfast=FALSE){
  # Function to compute a root inside a given interval
  # Function to compute integrals
  sumint=function (x, y, j){
    dx = diff(x[1:j], 1, 1)
    fx = y[1:j] 
    strap = 0.5 * (fx[1:(j - 1)] + fx[2:j])
    sj = sum(dx * strap)
    c(j, x[j],sj)
  }
  # Lagrange polynomial of the 2nd order
  polint2=function(x,x1,x2,x3,y1,y2,y3){
    y1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3))+y2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3))+y3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2))
  }
  n=length(y)
  # Check for insufficient number of points
  if(n<7){
    stop('Number of data points must be at least 7. Please increase the number of points.')
  }
  # Keep y values
  yi=y
  # Check if all y-values are nonzero
  if(sum(y>0)==n | sum(y<0)==n){
    if(!silent){message(paste0('It seems that no roots exist at [',x[1],' , ',x[n],']'))}
    out=c("x1"=x[1],"x2"=x[n],"chi"=NA,"yvalue"=NA)
    return(out)
  }
  # Take absolute values and convert root to extreme
  y2=abs(y)-y
  jnozero=which(y2!=0)
  if(length(jnozero)!=0){y[jnozero]=-y[jnozero]}
  # Main run
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
    fs=parallel::parSapply(cl=cl,2:(n-1),runfind,x,y)
    stopCluster(cl)
    tp2=Sys.time();print(tp2-tp1) 
    t12=paste0("Time for computing surfaces was ",format(tp2-tp1,units="sec"))
    if(!silent){message(t12)}
  }else{
    t1=Sys.time() 
    fs=sapply(2:(n-1), function(j,x,y){sumint(x,y,j)[3]},x,y)
    t2=Sys.time()
    t12=paste0("Time for computing surfaces was ",format(t2-t1,units="sec"))
    if(!silent){message(t12)}
  }
  # Find inflection point now if exists...
  # Try BEDE for fast run, although not so accurate as BESE:
    if(tryfast){
      t1=Sys.time() 
      # dd=bede(1:length(fs),fs,cc2$index)
      dd=bede(1:length(fs),fs,1)
      t2=Sys.time()
      t12=paste0("Time for computing BEDE was ",format(t2-t1,units="sec"))
      if(!silent){message(t12)}
      # Check convergence for BEDE
      if(is.nan(dd$iters$EDE[dim(dd$iters)[1]])){
        if(!silent){message(paste0('BEDE was not able to find an inflection at [',x[1],' , ',x[n],'] \n Continue by using BESE...'))}
      }else{
        solint=x[t(dd$iters[dim(dd$iters)[1],c("a","b")]+1)]  #Plus 1!
        sol=mean(solint)
        jj=c(unlist(dd$iters[dim(dd$iters)[1],c("a","b","EDE")]));jj
        j1=jj[1];j3=jj[2];j2=ifelse(j3-j1!=1,round(jj[3]),j3+1);c(j1,j2,j3)
        yval=polint2(sol,x[j1],x[j2],x[j3],yi[j1],yi[j2],yi[j3]);yval
        out=c(solint,sol,yval)
        names(out)=c("x1","x2","chi","yvalue")
        return(out)
      }
    }
    # Try BESE if BEDE has failed
    t1=Sys.time() 
    # bb=bese(1:length(fs),fs,cc2$index,doparallel = parallel)
    bb=bese(1:length(fs),fs,1,doparallel = parallel)
    t2=Sys.time()
    t12=paste0("Time for computing BESE was ",format(t2-t1,units="sec"))
    if(!silent){message(t12)}
    # Check existence of local extremes and for roots at interval edges...
    if(is.nan(bb$iplast)){
      if(yi[1]==0){
        out=c(x[1],NA,x[1],yi[1])
        names(out)=c("x1","x2","chi","yvalue")
        return(out)
      }else if(yi[n]==0){
        out=c(NA,x[n],x[n],yi[n])
        names(out)=c("x1","x2","chi","yvalue")
        return(out)
      }else{
        if(!silent){message(paste0('It seems that no roots exist at [',x[1],' , ',x[n],'] or multiple roots exist'))}
        out=c("x1"=x[1],"x2"=x[n],"chi"=NA,"yvalue"=NA)
        return(out)
      }
    }else{
      solint=x[t(bb$iters[dim(bb$iters)[1],c("a","b")]+1)]  #Plus 1!
      sol=mean(solint)
      jj=c(unlist(bb$iters[dim(bb$iters)[1],c("a","b","ESE")]));jj
      j1=jj[1];j3=jj[2];j2=ifelse(j3-j1!=1,round(jj[3]),j3+1);c(j1,j2,j3)
      yval=polint2(sol,x[j1],x[j2],x[j3],yi[j1],yi[j2],yi[j3]);yval
      out=c(solint,sol,yval)
      names(out)=c("x1","x2","chi","yvalue")
      return(out)
    }
}