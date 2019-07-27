classify_curve=function(x,y){
  # Check for infinities
  infp=which(y==Inf)
  if(length(infp)!=0){
    mssg=paste0('Infinities were detected for x-points ',paste0(infp,collapse=","))
    message(mssg)
    message(paste0(x[infp],collapse=","))
    out=list("ctype"=NA,"index"=NA,"asymmetry"=NA,
             "totalconvexity"=NA,"ismax"=NA,"shapetype"=NA)
    return(out)
  }else{
    # Function to decide for the convexity type of the curve
    N=length(x)
    #
    if(N>21){
      quants=seq(0,1,0.05)
    }else{
      quants=seq(0,1,0.1)
    }
    #
    jj=sapply(quants, function(q,N){as.integer(quantile(1:N,q))},N)
    #
    LR=t(sapply(jj, function(j,x,y){findipl(x,y,j)[3:4]},x,y))
    LR
    #
    sr0=LR[1,2]
    sln=LR[dim(LR)[1],1]
    #
    sleft=LR[2:dim(LR)[1],1]
    sumleft=sum(sleft,na.rm = TRUE)
    sumleft
    sright=LR[1:(dim(LR)[1]-1),2]
    sumright=sum(sright,na.rm = TRUE)
    sumright
    #
    leftsigns=sign(sleft)
    uleft=unique(leftsigns)
    rightsigns=sign(sright)
    uright=unique(rightsigns)
    #
    leftsigns
    rightsigns
    # Check convexity or at least convexity at the beginning
    # Left
    if(length(uleft)==1){
      ifelse(uleft>0,{cleft="concave"},{cleft="convex"})
    }else{
      ifelse(sum(head(leftsigns,3))>0,{cleft="concave"},{cleft="convex"})
    }
    # Right
    if(length(uright)==1){
      ifelse(uright>0,{cright="concave"},{cright="convex"})
    }else{
      ifelse(sum(tail(rightsigns,3))>0,{cright="concave"},{cright="convex"})
    }
    # Total
    if(cleft=="convex" & cright=="concave"){
      ctype="convex_concave"
      index=0
    }else if(cleft=="concave" & cright=="convex"){
      ctype="concave_convex"
      index=1
    }else if(cleft=="convex" & cright=="convex"){
      ctype="convex"
      index=0
    }else if(cleft=="concave" & cright=="concave"){
      ctype="concave"
      index=1
    }
    # Check for data asymmetry
    nposl=sum(leftsigns>0)
    nnegl=sum(leftsigns<0)
    nposr=sum(rightsigns>0)
    nnegr=sum(rightsigns<0)
    if(ctype=="convex" | ctype=="convex_concave"){
      if(nnegl==nposr){
        asym="data_symmetry"
      } else if(nnegl<nposr){
        asym="data_left_asymmetry"
      } else{
        asym="data_right_asymmetry"
      }
    } else {
      if(nposl==nnegr){
        asym="data_symmetry"
      } else if(nposl<nnegr){
        asym="data_left_asymmetry"
      } else{
        asym="data_right_asymmetry"
      }
    }
    # Check if curve has a maximum
    if(sr0>0 & sln>0){
      ismax=TRUE
    }else if(sr0<0 & sln<0){
      ismax=FALSE
    }else{
      ismax=NA
    }
    # Check if it is overall concave (for maximum computation) or overall convex (for minimum computation)
    sint=findipl(x,y,N)[3]
    ifelse(sint>0,{totalconvexity="concave"},{totalconvexity="convex"})
    # Check shape type
    ifelse(sum(-1==leftsigns)>=3 | sum(-1==rightsigns)>=3,{stype="bell"},{stype="tulip"})
    # Return all
    out=list("ctype"=ctype,"index"=index,"asymmetry"=asym,
             "totalconvexity"=totalconvexity,"ismax"=ismax,"shapetype"=stype)
    return(out)
  }
}