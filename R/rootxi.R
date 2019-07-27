rootxi <-
function(x,y,i1,i2,nt,alpha=5,xlb="x",ylb="y",xnd=3,ynd=3,plots=TRUE,plotpdf=FALSE,doparallel=FALSE){
  #Find root of discrete data (xi,yi) with Taylor regression and degree of polynomial nt
  #Choose the desired range i1,i2 of data in order to search at interval [x_i1,x_i2] 
  #Choose the level of statistical significance as a number, ie 5 means 5%
  #Plot results with xlb=label for x-axis, ylb=label for y-axis
  #Set number of digits for x-axis (xnd) and for y-axis (ynd)
  #Check for proper polynomial degree nt:
  #Check for nt negative or less than one:
  if(nt<0 | nt<1){stop('The argument "nt" must be nt >= 1')}
  #Check for nt non integer:
  if(nt-round(nt)!=0.0){warning('The argument "nt" must be integer, now it will be rounded to closest integer');nt=round(nt)}
  #Check for sufficient number of xy-points:
  if(length(i1:i2)<(nt+2)){stop('Number of xy-points is not sufficient for regression. Try decreasing nt...')}
  #Initialize:
  lev=(100-alpha)/100;
  x1<-x[i1:i2];y1<-y[i1:i2];
  #Function to Search for all available root-points p=x_i
  i=NULL;
  #
  fs=function(i,x1,y1){
    pr1<-x1[i];
    xm<-cbind();df<-NULL;
    for (j in (1:nt)){xm<-cbind(xm,cbind((x1-pr1)^j))}
    df<-as.data.frame(xm,row.names = NULL, optional = FALSE)
    xnam <- paste0("V", 1:(nt+0))
    fmla <- as.formula(paste("y1 ~ ", paste(xnam, collapse= "+")))
    c1<-lm(fmla,data=df);
    yp1<-predict(c1);
    ci1<-confint(c1,level = lev);
    am<-matrix(ci1,ncol=2,dimnames=list(c(paste0("a",0:nt)),c(colnames(ci1))));
    v0n=as.double(abs(c1$coeff[1]))
    points=cbind(v0n,pr1)
    names(points)=c("a0","x0")
    out=list(c(v0n,pr1),yp1,am)
    names(out)=c("points","predicted","confint")
    return(out)
  }
  #
  #Do parallel computing on request only:
  #
  if(doparallel){
    ncores=detectCores();
    cat(paste0('Available workers are ',ncores),'\n')
    #
    t1=Sys.time();
    cl <- makeCluster(ncores);
    registerDoParallel(cl)
    m3=foreach(i=1:length(i1:i2)) %dopar% {fs(i,x1,y1)};
    stopCluster(cl)
    t2=Sys.time();print(as.POSIXlt(t2, "GMT")-as.POSIXlt(t1, "GMT"),quote=F);#Time difference of 13.20788 secs #It worked!
    #
  }else{
    m3=lapply(1:length(i1:i2),fs,x1=x1,y1=y1)
  }
  #
  #Process results...
  mans=as.data.frame(do.call(rbind,sapply(m3,function(x){x[1]})))
  colnames(mans)=c("a0","x0")
  rownames(mans)=1:dim(mans)[1]
  #Find minimum value |a_0| and corresponding root x_i=p
  n0<-which.min(mans$a0);
  ipc<-mans[n0,"x0"];
  ys=as.data.frame(do.call(cbind,sapply(m3,function(x){x[2]})))
  yf1<-ys[,n0];
  v0n=mans$a0
  vpr=mans$x0
  ################
  #Plot results if plots
  if(plots){
  if(plotpdf){pdf('rootplot.pdf')}
  par(mfrow=c(1,2))
  ymin<-min(y1);
  ymax<-max(y1);
  dy<-(ymax-ymin)/5;
  plot(x1,y1,ylim=range(y1),col='blue',pch=19,cex=.5,axes=FALSE,ylab=ylb ,xlab=xlb);
  lines(x1,yf1,col='black',lty=1,lwd=2)
  dx<-(x1[length(x1)]-x1[1])/5;
  xticks<-round(c(x1[1],x1[n0],x1[length(x1)]),digits=xnd)
  axis(1,at=xticks,cex.axis=0.7,las=2)
  yticks<-round((seq(ymin,ymax,by=dy)),digits=ynd);
  axis(2,at=yticks,cex.axis=0.7)
  abline(h=yticks, v=xticks, col="gray", lty=3)
  abline(v=x1[n0],lty=2,col="red",lwd=3)
  legend('top',col=c('blue'),pch=c(19),legend=c('data'),bty='n',cex=0.6);
  legend('left',col=c('black','red'),lty=c(1,2),lwd=c(2,3),legend=c(paste0('Taylor fit (',nt,' )'),'root'),bty='n',cex=0.7)
  stit<-paste0('Data for [',toString(round(x1[1],digits=2)),',',toString(round(x1[length(x1)],digits=xnd)),'] \n');
  title(paste0(stit,'Taylor Regression n = ',nt,' , a = ',alpha,' % ','\n Find root'),cex.main=0.7)
  box()
  plot(vpr,v0n,col='blue',pch=19,cex=.5,axes=FALSE,ylab=expression(paste('|',alpha[0],'|')) ,xlab=expression(rho));
  abline(v=ipc,lty=3,col="blue",lwd=2)
  dp1<-(ipc-vpr[1])/2;
  dp2<-(vpr[length(vpr)]-ipc)/3;
  xticks<-round(c(vpr[1],ipc,vpr[length(vpr)]),digits=xnd)
  axis(1,at=xticks,cex.axis=0.7,las=2)
  dvn0<-(max(v0n)-min(v0n))/5;
  yticks<-c(seq(min(v0n),max(v0n),by=dvn0))
  axis(2,at=yticks,cex.axis=0.7,labels=round(yticks,digits=ynd))
  abline(h=yticks, v=xticks, col="gray", lty=3)
  legend('top',col=c('blue'),pch=c(19),legend=c(expression(paste('|',alpha[0],'|'))),bty='n',cex=0.6);
  title(main=expression(paste('Plot of all available |',alpha[0],'|')),cex.main=0.7)  
  box() 
  if(plotpdf){dev.off()}
  par(mfrow=c(1,1))
  }
  if(plotpdf){cat("File 'rootplot.pdf' has been created","\n");dev.off()}
  #Return output...
  amf=m3[[n0]]$confint;amf
  ans<-new.env();
  ans$an<-matrix(cbind(amf,0.5*rowSums(amf)),nrow=nt+1,ncol=3,byrow=F,list(c(paste0("a",c(0:nt))),c(colnames(amf),"an")));
  ans$froot<-c(n0,vpr[n0]);
  return(ans)
}