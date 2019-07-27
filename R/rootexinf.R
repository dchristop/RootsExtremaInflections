rootexinf <-
function(x,y,i1,i2,nt,alpha=5,xlb="x",ylb="y",xnd=3,ynd=3,plots=TRUE,plotpdf=FALSE,doparallel=FALSE){
  #Find root, extreme and inflection point for discrete data (xi,yi) with Taylor regression and degree of polynomial nt
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
  #Function to Search for all available root, extreme and inflection points 
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
    #root
    v0n0=as.double(abs(c1$coeff[1]))
    #extreme
    v0n1=as.double(abs(c1$coeff[2]))
    #inflection
    v0n2=as.double(abs(c1$coeff[3]))
    #all
    points=c(pr1,v0n0,v0n1,v0n2)
    out=list(points,yp1,am)
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
    t2=Sys.time();print(as.POSIXlt(t2, "GMT")-as.POSIXlt(t1, "GMT"),quote=F);#Time difference of 13.64839 secs
    #
  }else{
    m3=lapply(1:length(i1:i2),fs,x1=x1,y1=y1)
  }
  #
  #Process results...
  mans=as.data.frame(do.call(rbind,sapply(m3,function(x){x[1]})))
  colnames(mans)=c("rho","a0","a1","a2")
  rownames(mans)=1:dim(mans)[1]
  #Find minimum values for |a_0|, |a_1|, |a_2| and corresponding data x-points 
  n03=apply(mans[,c("a0","a1","a2")],2,which.min)
  # n0<-which.min(mans$a0);
  ipc3<-mans[n03,"rho"];
  ys=as.data.frame(do.call(cbind,sapply(m3,function(x){x[2]})))
  yf13<-ys[,n03];
  v0n3=mans[,c("a0","a1","a2")]
  vpr=mans$rho
  ################
  #Plot results if plots
  if(plots){
  if(plotpdf){pdf('root_extreme_inflection_plot.pdf')}
  par(mfrow=c(3,2))
  xmin<-min(x1);xmax<-max(x1);dx<-(xmax-xmin)/5;xticks<-round((seq(xmin,xmax,by=dx)),digits=xnd);
  ymin<-min(y1);ymax<-max(y1);dy<-(ymax-ymin)/5;yticks<-round((seq(ymin,ymax,by=dy)),digits=ynd);
  # dx<-(x1[length(x1)]-x1[1])/5;xticks<-round(c(x1[1],x1[n03[1]],x1[length(x1)]),digits=xnd)
  #root
  plot(x1,y1,ylim=range(y1),col='blue',pch=19,cex=.5,axes=FALSE,ylab=ylb ,xlab=xlb);
  lines(x1,yf13[,1],col='black',lty=1,lwd=2)
  axis(1,at=xticks,cex.axis=0.7,las=2);axis(2,at=yticks,cex.axis=0.7);abline(h=yticks, v=xticks, col="gray", lty=3)
  abline(v=x1[n03[1]],lty=2,col="red",lwd=3)
  legend('top',col=c('blue'),pch=c(19),legend=c('data'),bty='n',cex=0.6);
  legend('left',col=c('black','red'),lty=c(1,2),lwd=c(2,3),legend=c(paste0('Taylor fit (',nt,' )'),'root'),bty='n',cex=0.7)
  stit<-paste0('Data for [',toString(round(x1[1],digits=2)),',',toString(round(x1[length(x1)],digits=xnd)),'] \n');
  title(paste0(stit,'Taylor Regression n = ',nt,' , a = ',alpha,' % ','\n Find root'),cex.main=0.7)
  box()
  plot(vpr,v0n3[,"a0"],col='blue',pch=19,cex=.5,axes=FALSE,ylab=expression(paste('|',alpha[0],'|')) ,xlab=expression(rho));
  abline(v=ipc3[1],lty=3,col="blue",lwd=2)
  dp1<-(ipc3[1]-vpr[1])/2;
  dp2<-(vpr[length(vpr)]-ipc3[1])/3;
  xticks1<-round(c(vpr[1],ipc3[1],vpr[length(vpr)]),digits=xnd)
  axis(1,at=xticks1,cex.axis=0.7,las=2)
  dvn0<-(max(v0n3[,"a0"])-min(v0n3[,"a0"]))/5;
  yticks1<-c(seq(min(v0n3[,"a0"]),max(v0n3[,"a0"]),by=dvn0))
  axis(2,at=yticks1,cex.axis=0.7,labels=round(yticks,digits=ynd))
  abline(h=yticks1, v=xticks1, col="gray", lty=3)
  legend('top',col=c('blue'),pch=c(19),legend=c(expression(paste('|',alpha[0],'|'))),bty='n',cex=0.6);
  title(main=expression(paste('Plot of all available |',alpha[0],'|')),cex.main=0.7)  
  box() 
  #extreme
  plot(x1,y1,ylim=range(y1),col='blue',pch=19,cex=.5,axes=FALSE,ylab=ylb ,xlab=xlb);
  axis(1,at=xticks,cex.axis=0.7,las=2);axis(2,at=yticks,cex.axis=0.7);abline(h=yticks, v=xticks, col="gray", lty=3)
  lines(x1,yf13[,2],col='black',lty=1,lwd=2)
  abline(v=x1[n03[2]],lty=2,col="red",lwd=3)
  legend('top',col=c('blue'),pch=c(19),legend=c('data'),bty='n',cex=0.6);
  legend('left',col=c('black','red'),lty=c(1,2),lwd=c(2,3),legend=c(paste0('Taylor fit (',nt,' )'),'extreme'),bty='n',cex=0.7)
  stit<-paste0('Data for [',toString(round(x1[1],digits=2)),',',toString(round(x1[length(x1)],digits=xnd)),'] \n');
  title(paste0(stit,'Taylor Regression n = ',nt,' , a = ',alpha,' % ','\n Find extreme'),cex.main=0.7)
  box()
  plot(vpr,v0n3[,"a1"],col='blue',pch=19,cex=.5,axes=FALSE,ylab=expression(paste('|',alpha[1],'|')) ,xlab=expression(rho));
  abline(v=ipc3[2],lty=3,col="blue",lwd=2)
  dp1<-(ipc3[2]-vpr[1])/2;
  dp2<-(vpr[length(vpr)]-ipc3[2])/3;
  xticks2<-round(c(vpr[1],ipc3[2],vpr[length(vpr)]),digits=xnd)
  axis(1,at=xticks2,cex.axis=0.7,las=2)
  dvn0<-(max(v0n3[,"a1"])-min(v0n3[,"a1"]))/5;
  yticks2<-c(seq(min(v0n3[,"a1"]),max(v0n3[,"a1"]),by=dvn0))
  axis(2,at=yticks2,cex.axis=0.7,labels=round(yticks,digits=ynd))
  abline(h=yticks2, v=xticks2, col="gray", lty=3)
  legend('top',col=c('blue'),pch=c(19),legend=c(expression(paste('|',alpha[1],'|'))),bty='n',cex=0.6);
  title(main=expression(paste('Plot of all available |',alpha[1],'|')),cex.main=0.7)  
  box() 
  #inflection
  plot(x1,y1,ylim=range(y1),col='blue',pch=19,cex=.5,axes=FALSE,ylab=ylb ,xlab=xlb);
  axis(1,at=xticks,cex.axis=0.7,las=2);axis(2,at=yticks,cex.axis=0.7);abline(h=yticks, v=xticks, col="gray", lty=3)
  lines(x1,yf13[,3],col='black',lty=1,lwd=2)
  abline(v=x1[n03[3]],lty=2,col="red",lwd=3)
  legend('top',col=c('blue'),pch=c(19),legend=c('data'),bty='n',cex=0.6);
  legend('left',col=c('black','red'),lty=c(1,2),lwd=c(2,3),legend=c(paste0('Taylor fit (',nt,' )'),'inflection'),bty='n',cex=0.7)
  stit<-paste0('Data for [',toString(round(x1[1],digits=2)),',',toString(round(x1[length(x1)],digits=xnd)),'] \n');
  title(paste0(stit,'Taylor Regression n = ',nt,' , a = ',alpha,' % ','\n Find inflection'),cex.main=0.7)
  box()
  plot(vpr,v0n3[,"a2"],col='blue',pch=19,cex=.5,axes=FALSE,ylab=expression(paste('|',alpha[2],'|')) ,xlab=expression(rho));
  abline(v=ipc3[3],lty=3,col="blue",lwd=2)
  dp1<-(ipc3[3]-vpr[1])/2;
  dp2<-(vpr[length(vpr)]-ipc3[3])/3;
  xticks3<-round(c(vpr[1],ipc3[3],vpr[length(vpr)]),digits=xnd)
  axis(1,at=xticks3,cex.axis=0.7,las=2)
  dvn0<-(max(v0n3[,"a2"])-min(v0n3[,"a2"]))/5;
  yticks3<-c(seq(min(v0n3[,"a2"]),max(v0n3[,"a2"]),by=dvn0))
  axis(2,at=yticks3,cex.axis=0.7,labels=round(yticks,digits=ynd))
  abline(h=yticks3, v=xticks3, col="gray", lty=3)
  legend('top',col=c('blue'),pch=c(19),legend=c(expression(paste('|',alpha[2],'|'))),bty='n',cex=0.6);
  title(main=expression(paste('Plot of all available |',alpha[2],'|')),cex.main=0.7)  
  box() 
  if(plotpdf){dev.off()}
  par(mfrow=c(1,1))
  }
  if(plotpdf){cat("File 'root_extreme_inflection_plot.pdf' has been created","\n");dev.off()}
  #Return output...
  amf1=m3[[n03[1]]]$confint;amf1
  amf2=m3[[n03[2]]]$confint;amf2
  amf3=m3[[n03[3]]]$confint;amf3
  ans<-new.env();
  #ci
  ans$an0<-matrix(cbind(amf1,0.5*rowSums(amf1)),nrow=nt+1,ncol=3,byrow=F,list(c(paste0("a",c(0:nt))),c(colnames(amf1),"an0")));
  ans$an1<-matrix(cbind(amf2,0.5*rowSums(amf2)),nrow=nt+1,ncol=3,byrow=F,list(c(paste0("a",c(0:nt))),c(colnames(amf2),"an1")));
  ans$an2<-matrix(cbind(amf3,0.5*rowSums(amf3)),nrow=nt+1,ncol=3,byrow=F,list(c(paste0("a",c(0:nt))),c(colnames(amf2),"an2")));
  #root, extreme, inflection
  frexinf=cbind(n03,vpr[n03]);rownames(frexinf)=c("root","extreme","inflection");colnames(frexinf)=c("index","value")
  ans$frexinf<-frexinf
  return(ans)
}