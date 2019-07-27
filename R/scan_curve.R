scan_curve=function(x,y,findroots=TRUE,findextremes=TRUE,findinflections=TRUE,silent=FALSE,plots=TRUE){
  # Functions
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
  # Transformation
  fabs=function(y){ifelse(y<0,-y,0)}
  # Monotonicity
  fmonoton=function(index){
    if(index==-1){
      mon="decreasing"
    }else if(index==1){
      mon="increasing"
    }else{
      mon=NA
    }
    return(mon)
  }
  # Extremity
  fextremity=function(d1,d2){
    if(d1==1 & d2==-1){
      out="max"
    }else if(d1==-1 & d2==1){
      out="min"
    }else{
      out=NA
    }
    return(out)
  }
  # Sigmoidicity
  fsigmoid=function(d1,d2){
    if(d1=="max" & d2=="min"){
      ctype="concave_convex"
      index=1
    }else if(d1=="min" & d2=="max"){
      ctype="convex_concave"
      index=0
    }else{
      ctype=NA
      index=NA
    }
    out=list("type"=ctype,"index"=index)
    return(out)
  }
  # Interpolation
  polint2=function(x,x1,x2,x3,y1,y2,y3){
    y1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3))+y2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3))+y3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2))
  }
  # Plot
  if(plots){plot(x,y,pch=19,cex=0.2);abline(h=0)}
  # Check for noise:
  ifelse(length(x)>10,{nlag=5},{nlag=3})
  noise=check_noise(y,n=nlag)
  if(noise){
    findroots=FALSE
    findextremes=FALSE
    findinflections=FALSE
    if(!silent){message('The curve seems to be a noisy one, so only a visualization is possible...')}
  }
  # Check for no roots:
  # Check if all y-values are nonzero
  n=length(x)
  if(sum(y>0)==n | sum(y<0)==n){
    if(!silent){message(paste0('It seems that no roots exist at [',x[1],' , ',x[n],']'))}
    roots=rep(NA,4);names(roots)=c("x1","x2","chi","yvalue")
  }
    # data frame
    df=data.frame("x"=x,"y"=y)
    rownames(df)=1:length(x)
    # Transform and store
    y2=fabs(y)
    df$y2=y2
    df$zero=rep(FALSE,dim(df)[1])
    df[df$y2==0,"zero"]=TRUE
    df$dif=rep(0,dim(df)[1])
    df$dif[1:(dim(df)[1]-1)]=diff(df$zero)
    # Flat regions
    izero=which(df$dif==0)
    # Cases:
    if(length(izero)==length(x)){
      jmin=which.min(y)
      xmin=x[jmin]
      ymin=y[jmin]
      xt=findextreme(x,y)
      ifelse(ymin==0,
             {
               roots=c(NA,NA,xmin,ymin)
               names(roots)=c("x1","x2","chi","yvalue")
             },
             {
               if(!silent){message(paste0('It seems that no roots exist at [',x[1],' , ',x[length(x)],']'))}
               roots=rep(NA,4)
               names(roots)=c("x1","x2","chi","yvalue")
               }
               )
      out=list("study"=NA,"roots"=roots,"extremes"=xt,"inflections"=NA)
      if(plots){points(xmin,ymin,pch=19,cex=1.5)}
      return(out)
    }else{
      # Roots:
      iroot=which(df$dif!=0)
      # Cases:
      if(length(iroot)==0){
        if(!silent){message(paste0('It seems that no roots exist at [',x[1],' , ',x[length(x)],']'))}
        out=list("study"=NA,"roots"=NA,"extremes"=NA,"inflections"=NA)
        return(out)
      }else{
        # Root intervals
        df$ja=rep(NA,dim(df)[1])
        df$jb=rep(NA,dim(df)[1])
        df[iroot[1],"ja"]=1
        df[iroot[1],"jb"]=iroot[1]
        #
        if(length(iroot)>1){
          for(j in 2:length(iroot)){
            df[iroot[j],"ja"]=iroot[j-1]+1
            df[iroot[j],"jb"]=iroot[j]
          }
        }
        #
        dfroot=df[c(iroot,dim(df)[1]),]
        dfroot[length(iroot)+1,"ja"]=iroot[length(iroot)]+1
        dfroot[length(iroot)+1,"jb"]=dim(df)[1]
        dfroot[length(iroot)+1,"dif"]=0
        # Monotonicity
        dfroot$root_monotonicity=sapply(dfroot$dif,fmonoton)
        dfroot
        # Extremity
        dfroot$extremity=rep(NA,dim(dfroot)[1])
        #
        dd=dfroot$dif[1:(dim(dfroot)[1]-1)]
        dd
        #
        if(length(dd)>1){
          extrs=sapply(2:length(dd), function(i,dd){fextremity(dd[i-1],dd[i])},dd)
          dfroot[2:length(dd),"extremity"]=extrs
        }else{
          extrs=NULL
        }
        #
        # Sigmoidicity
        dfroot$sigmoidicity=rep(NA,dim(dfroot)[1])
        # Only roots:
        droot=dfroot[1:dim(dfroot)[1]-1,]        
        #
        # Find all roots if asked so:
        if(findroots){
          j1=c(1,iroot,dim(df)[1]);j1
          if(length(j1)>1){
            j2=sapply(2:length(j1), function(i,j1){round(mean(j1[(i-1):i]))},j1);j2
            rt=list()
            for(k in 1:(length(j2)-1))   
            {
              x1=x[j2[k]:j2[k+1]]
              y1=y[j2[k]:j2[k+1]]
              if(length(x1)>=7){
                rt[[k]]=findroot(x1,y1)
              }else{
                if(!silent){message(paste0('Insufficient number of points for computing root at [',x1[1],' , ',x1[length(x1)],']'))}
                rt[[k]]=c("x1"=x1[1],x2=x1[length(x1)],"chi"=NA,"yvalue"=NA)
              }
            }
            roots=do.call(rbind,rt)
          }else{
            k1=round(mean(j1[1:2]))
            k2=round(mean(j1[2:3]))
            x1=x[k1:k2]
            y1=y[k1:k2]
            roots=findroot(x1,y1)
          }
        }else{
          roots=NA
        }
        #
        # Find extremes if asked so:
        if(findextremes){
          dextrs=droot[!is.na(droot$extremity),]
          if(dim(dextrs)[1]!=0){
            xt=t(sapply(1:dim(dextrs)[1], function(i,dextrs){
              irange=dextrs[i,"ja"]:dextrs[i,"jb"]
              if(length(irange)>=7){
                out=findextreme(x[irange],y[irange])
                return(out)
              }else{
                out=c("x1"=x[irange[1]],x2=x[irange[2]],"chi"=NA,"yvalue"=NA)
                return(out)
              }
            },dextrs))
            extremes=xt
          }else{
            if(!silent){message(paste0('It seems that no extremes between roots exist at [',x[1],' , ',x[length(x)],']'))}
            extremes=NA
          }
        }else{
          extremes=NA
        }
        #
        # Find inflections if asked so:
        if(findinflections){
          dextrs=droot[!is.na(droot$extremity),]
          dextrs
          if(dim(dextrs)[1]>1){
            pinf=t(sapply(1:(dim(dextrs)[1]-1), function(i,dextrs){
              cs=fsigmoid(dextrs[i,"extremity"],dextrs[i+1,"extremity"])
              cs
              irange=dextrs[i,"ja"]:dextrs[i+1,"jb"]
              irange
              if(length(irange)>=4){
                x1=x[irange];y1=y[irange]
                bb=bese(x1,y1,cs$index)
                bb
                #
                solint=t(bb$iters[dim(bb$iters)[1],c("a","b")])
                sol=mean(solint)
                nint=bb$iters[dim(bb$iters)[1],"n"]
                xab=bb$iters[dim(bb$iters)[1],c("a","b")]
                xa=as.double(xab[1]);j1=which(x1==xa)
                j2=j1+1
                xb=as.double(xab[2]);j3=which(x1==xb)
                c(j1,j2,j3)
                yval=polint2(sol,x1[j1],x1[j2],x1[j3],y1[j1],y1[j2],y1[j3]);yval
                out=c(solint,sol,yval)
                names(out)=c("x1","x2","chi","yvalue")
                return(out)
              }else{
                out=c("x1"=x[irange[1]],x2=x[irange[2]],"chi"=NA,"yvalue"=NA)
                return(out)
              }
            },dextrs))
            pinf
            inflects=pinf
            infs=data.frame(t(sapply(2:length(dextrs$extremity), function(i,extrs){fsigmoid(extrs[i-1],extrs[i])},extrs=dextrs$extremity)))
            dextrs[2:dim(dextrs)[1],"sigmoidicity" ]=do.call(c, infs$type)
            droot[2:dim(droot)[1],]=dextrs
          }else{
            if(!silent){message(paste0('It seems that no inflections between roots exist at [',x[1],' , ',x[length(x)],']'))}
            inflects=NA
          }
        }else{
          inflects=NA
        }
        #
        # Plots
        if(plots){
          #
          for(i in 1:dim(droot)[1]){
            # Root points
            if(dfroot[i,"dif"]==1){
              ptpch=24
            }else if(dfroot[i,"dif"]==-1){
              ptpch=25
            }else{
              ptpch=23
            }
            # Extreme and inflectionpoints
            if(is.na(dfroot[i,"extremity"])){
              colr="black"
            }else if(dfroot[i,"extremity"]=="max"){
              colr="red"
            }else if(dfroot[i,"extremity"]=="min"){
              colr="blue"
            }else{
              colr="black"
            }
            #
            lines(x[dfroot$ja[i]:dfroot$jb[i]],y[dfroot$ja[i]:dfroot$jb[i]],col=colr,lwd=3)
            points(dfroot$x[i],dfroot$y[i],pch=ptpch,cex=1.0)
            #
            if(findextremes){
              if(i>1){
                j=i-1
                segments(x0=xt[j,'chi'],y0=0,x1=xt[j,'chi'],y1=xt[j,'yvalue'],lty=2,col=colr)
                points(xt[j,'chi'],xt[j,'yvalue'],col='black',cex=1,pch=19)
              }
            }
            #
            if(findinflections){
              if(i>2){
                j=i-2
                segments(x0=pinf[j,'chi'],y0=0,x1=pinf[j,'chi'],y1=pinf[j,'yvalue'],lty=2,col=colr)
                points(pinf[j,'chi'],pinf[j,'yvalue'],col='black',cex=1,pch=18)
              }
            }
          }
          #
        }
        #
        out=list("study"=droot,"roots"=roots,"extremes"=extremes,"inflections"=inflects)
        return(out)
      }
    }
}
