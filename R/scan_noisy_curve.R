scan_noisy_curve=function(x,y,noise=NULL,rootsoptim=TRUE,findextremes=TRUE,findinflections=TRUE,silent=FALSE,plots=TRUE){
  # Function for computing roots, extrema and inflections for a noisy curve
  # Useful functions
  round2 <- function(x) { trunc(x + sign(x) * 0.5)}
  # Check if curve is without error
  check_noise=function(y,nlag){
    sds=sapply(1:nlag, function(i,y){round(sd(diff(y,i)),6)},y=y)
    if(sds[1]==0){
      return(FALSE)
    }else{
      sumsds=sum(round(sds/sds[1],1))
      nseq=nlag*(nlag+1)/2
      ifelse(abs(sumsds-nseq)>(nlag-2),{return(TRUE)},{return(FALSE)})
    }
  }
  # Interpolation
  polint2=function(x,x1,x2,x3,y1,y2,y3){
    y1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3))+y2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3))+y3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2))
  }
  # Check for noise if not given true
  if(is.null(noise)){
    ifelse(length(x)>10,{nlag=5},{nlag=3})
    noise=check_noise(y,nlag)
  }
  #
  if(noise){
    scan1=scan_curve(x,y,findroots = FALSE,findextremes = FALSE,findinflections = FALSE,silent=TRUE,plots=plots)
    sc=scan1$study
    sc
    if(length(sc)==1){
      return(scan1)
    }else{
      if(plots){points(sc$x,sc$y,col='gold',pch=19)}
      jj=as.integer(rownames(sc))      
      dj=as.integer(diff(jj))
      djlogs=round(log10(dj))
      djlogsu=unique(djlogs)
      if(sum(djlogsu)==1){one=TRUE}else{one=FALSE}
      #
      ###
      if(length(dj)!=0){
        #
        dd=data.frame("j"=jj[1:(length(jj)-1)],"dj"=dj)
        dd=dd[order(dd$dj,decreasing = T),]
        dd2=dd[dd$dj>1,]        
        if(dim(dd2)[1]>=5){
          max2=max(dd2$dj);max2
          dj2=(dd2$dj-max2)^2;dj2          
          nint=ede(1:length(dj2),dj2,0)[2]-1
        }else{
          nint=1
        }
        dint=dd[1:nint,]
        dint
        if(nint==1 | one){
          onlyone=TRUE
          jint=dint$j
          djint=dint$dj
          dr=data.frame("j"=jint,"dj"=djint)
          dr[,"interval"]=FALSE
          dr[,"i1"]=dr[,"j"]
          dr[,"i2"]=dr[,"j"]+dr[,"dj"]
          dr$root=rep(TRUE,dim(dr)[1])
          dr
          sol1=x[jint];yval1=y[jint];c(sol1,yval1)
          if(dim(dr)[1]==1){
            sol2=mean(x[dr$i1:dr$i2]);yval2=mean(y[dr$i1:dr$i2]);c(sol2,yval2)
            sol3=x[jint+djint];yval3=y[jint+djint];c(sol3,yval3)
          }else{
            sol2=mean(x[jj]);yval2=mean(y[jj]);c(sol2,yval2)
            sol3=mean(x[jint+djint]);yval3=mean(y[jint+djint]);c(sol3,yval3)
          }          
          vlogs=abs(round2(log10(abs(c(sol1,sol2,sol3)))));vlogs
          if(length(unique(vlogs))==1 | one){
            sol=mean(x[jj]);yval=mean(y[jj]);c(sol,yval)
            xyroots=data.frame("x1"=x[jj[1]],"x2"=x[jj[length(jj)]],"chi"=sol,"yvalue"=yval)
          }else{
            xyroots=data.frame("x1"=x[jj[1]],"x2"=x[jj[length(jj)]],"chi"=c(sol1,sol2,sol3),"yvalue"=c(yval1,yval2,yval3))
            xyroots=xyroots[abs(xyroots$yvalue)<1,]
          }
          xyroots
          #
          if(plots){
            abline(v=xyroots$chi,lty=2)
            points(xyroots$chi,xyroots$yvalue,col='green',pch=20,cex=2)
          }
        }else{
          onlyone=FALSE
          dd=dd[order(dd$j),]
          dd$interval=rep(FALSE,dim(dd)[1])
          dd[rownames(dd)%in%rownames(dint),"interval"]=TRUE
          dd$i1=rep(NA,dim(dd)[1])
          dd$i2=rep(NA,dim(dd)[1])
          dd[dd$interval,"i1"]=dd[dd$interval,"j"]
          dd[dd$interval,"i2"]=dd[dd$interval,"j"]+dd[dd$interval,"dj"]
          dd$root=rep(FALSE,dim(dd)[1])
          dd[!dd$interval,"root"]=TRUE
          dd
          dr=dd[dd$interval,]
          dr
          # Find roots, attempt 1
          if(dim(dr)[1]>=1 & !onlyone){
            # List of "roots" around a real root...
            #######################################
            lroots=list()
            lroots[[1]]=dd$j[1]:dr[1,"i1"]
            if(dim(dr)[1]>1){
              for(k in 2:dim(dr)[1]){
                lroots[[k]]=dr[k-1,"i2"]:dr[k,"i1"]
              }
            }
            lroots[[length(lroots)+1]]=dr[dim(dr)[1],"i2"]:dd[dim(dd)[1],'j']
            # Find average roots:
            xyroots=data.frame(t(sapply(lroots, function(lr,x,y){
              out=c(x[lr[1]],x[lr[length(lr)]],mean(x[lr]),mean(y[lr]))
              names(out)=c("x1","x2","chi","yvalue")
              return(out)
            },x=x,y=y)))
            #
            if(plots){points(xyroots$chi,xyroots$yvalue,col='green',pch=20,cex=2)}
            #
          }else{
            if(!silent){message(paste0('It seems that no roots exist at [',x[1],' , ',x[length(x)],']'))}
            xyroots=NA
          }
          #ok roots 1st attempt done.
        }
        #
        if(rootsoptim & !onlyone){
          # roots 2nd attempt now...use optimization
          # 1st root
          irangefirst=round(mean(1:dr[1,"i1"])):round(mean(dr[1,"i1"]:dr[1,"i2"]))
          x1first=x[irangefirst];y1first=y[irangefirst]
          ifelse(length(irangefirst)>=7,{rfirst=findroot(x1first,y1first)},
                 {rfirst=c("x1"=x[irangefirst[1]],x2=x[irangefirst[2]],"chi"=NA,"yvalue"=NA) })
          # indermediate roots
          rootsinter=t(sapply(1:(dim(dr)[1]-1), function(i,dextrs){
            irange=round(mean(dr[i,"i1"]:dr[i,"i2"])):round(mean(dr[i+1,"i1"]:dr[i+1,"i2"]))
            if(length(irange)>=7){
              x1=x[irange];y1=y[irange]
              cc=check_curve(x1,y1)
              cc
              out=findroot(x1,y1)
              return(out)
            }else{
              out=c("x1"=x[irange[1]],x2=x[irange[2]],"chi"=NA,"yvalue"=NA)
              return(out)
            }
          },dr))
          rootsinter
          # last root
          irangelast=round(mean(dr[dim(dr)[1],"i1"]:dr[dim(dr)[1],"i2"])):round(mean(dr[dim(dr)[1],"i2"]:length(x)))
          x1last=x[irangelast];y1last=y[irangelast]
          ifelse(length(irangelast)>=7,{rlast=findroot(x1last,y1last)},
                 {rlast=c("x1"=x[irangelast[1]],x2=x[irangelast[2]],"chi"=NA,"yvalue"=NA) })
          #
          xyroots2=as.data.frame(rbind(rfirst,rootsinter,rlast))
          rownames(xyroots2)=1:dim(xyroots2)[1]
          colnames(xyroots2)=c("x1","x2","chi","yvalue")
          xyroots2
          #
          if(plots){points(xyroots2$chi,xyroots2$yvalue,col='green',pch="|",cex=2)}
        }else if(rootsoptim & onlyone){
          xyroots2=findroot(x,y)
        }
        #
        if(findextremes){
            # Find extreme points
            if(dim(dr)[1]>=1 &!onlyone){
              xt=t(sapply(1:dim(dr)[1], function(i,dr){
                irange=dr[i,"i1"]:dr[i,"i2"];irange
                if(length(irange)>=7){
                  out=findextreme(x[irange],y[irange])
                  return(out)
                }else{
                  out=c("x1"=x[irange[1]],x2=x[irange[2]],"chi"=NA,"yvalue"=NA)
                  return(out)
                }
              },dr))
              #
              if(plots & sum(sapply(xt[1,],is.na))==0){
                points(xt[,"chi"],xt[,"yvalue"],col='black',pch=15,cex=2)
                xtx=xt[,'chi'];xty=xt[,'yvalue']
                segments(x0=xtx,y0=rep(0,length(xtx)),x1=xtx,y1=xty,lty=2,col='black')
              }
              #
            }else{
              if(!silent){message(paste0('It seems that no extreme between roots exist at [',x[1],' , ',x[length(x)],']'))}
              xt=NA
            }
        }else{
          xt=NA
        }
        # ok extremes found.
        #
        if(findinflections){  
            # Find inflection points
            if(dim(dr)[1]>1 &!onlyone){
              pinf=t(sapply(1:(dim(dr)[1]-1), function(i,dextrs){
                irange=dr[i,"i1"]:dr[i+1,"i2"]
                if(length(irange)>=4){
                  x1=x[irange];y1=y[irange]
                  cc=check_curve(x1,y1)
                  cc
                  xd=ede(x1,y1,cc$index)
                  solint=x1[xd[1:2]]
                  sol=xd[3]
                  j1=xd[1,1]
                  j3=xd[1,2]
                  j2=round(mean(c(j1,j3)))
                  c(j1,j2,j3)
                  yval=polint2(sol,x1[j1],x1[j2],x1[j3],y1[j1],y1[j2],y1[j3]);yval
                  out=c(solint,sol,yval)
                  names(out)=c("x1","x2","chi","yvalue")
                  return(out)
                }else{
                  out=c("x1"=x[irange[1]],x2=x[irange[2]],"chi"=NA,"yvalue"=NA)
                  return(out)
                }
              },dr))
              pinf
              inflects=pinf
              #
              if(plots){
                points(inflects[,"chi"],inflects[,"yvalue"],col='black',pch=18,cex=2)
                pinfx=pinf[,'chi'];pinfy=pinf[,'yvalue']
                segments(x0=pinfx,y0=rep(0,length(pinfx)),x1=pinfx,y1=pinfy,lty=2,col='black')
                
              }
              #
            }else{
              if(!silent){message(paste0('It seems that no inflections between roots exist at [',x[1],' , ',x[length(x)],']'))}
              inflects=NA
            }
        }else{
          inflects=NA
        }
        # ok inflections found and roots 2nd attempt done.
        ###
        # Return all results:
        out=list("study"=dr,"roots_average"=xyroots,"roots_optim"=xyroots2,"extremes"=xt,"inflections"=inflects)
        #
        return(out)
      }else{
        # return(sc)
        dr=data.frame("j"=sc[1,7],"dj"=0,"i1"=1,"i2"=sc[1,7],"root"=TRUE)
        rownames(dr)=sc[1,7]
        xyroots=data.frame("chi"=sc[1,1],"yvalue"=sc[1,2])
        xr2=findroot(x,y);xr2=matrix(xr2,ncol=4,nrow=1)
        xyroots2=as.data.frame(xr2)
        colnames(xyroots2)=c("x1","x2","chi","yvalue")
        if(plots){
          points(xyroots$chi,xyroots$yvalue,col='green',pch=20,cex=2)
          points(xyroots2$chi,xyroots2$yvalue,col='green',pch="|",cex=2)
          abline(v=c(xyroots$chi,xyroots2$chi),lty=2)
        }
        out=list("study"=dr,"roots_average"=xyroots,"roots_optim"=xyroots2,"extremes"=NA,"inflections"=NA)
        #
        return(out)
      }
    }
  }else{
    if(!silent){message("Curve seems to be not a noisy one, use of 'scan_curve()' now...")}
    out=scan_curve(x,y,findroots = TRUE,findextremes = findextremes,findinflections =findinflections,silent=silent,plots=plots)
    #
    return(out)
  }
}
