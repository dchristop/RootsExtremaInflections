symextreme<-function(x,y,concave=NULL,type=NULL){
  # Function to compute symmetrical extremes
  # Find concavity of curve and determine the type of extreme
  if(is.null(concave) & is.null(type)){
    cc=classify_curve(x,y)
    concave1=cc$totalconvexity
    # Find shape of curve
    type=cc$shapetype
    # Set concavity of curve
    ifelse(concave1=="concave",{concave=TRUE},{concave=FALSE}) 
  }else if(is.null(concave)){
    cc=classify_curve(x,y)
    concave1=cc$totalconvexity
    # Set concavity of curve
    ifelse(concave1=="concave",{concave=TRUE},{concave=FALSE}) 
  }else if(is.null(type)) {
    # Find shape of curve
    cc=classify_curve(x,y)
    type=cc$shapetype
  }
  # Set extreme type
  ifelse(concave,{ismax=TRUE;ismin=FALSE},{ismax=FALSE;ismin=TRUE})
  # Compute extreme
  if(type=="bell")
  {
    a<-findmaxbell(x,y,concave)    
  }
  else if(type=="tulip")
  {
    a<-findmaxtulip(x,y,concave)    
  }else{
    stop("You must provide type='bell' or type='tulip'")
  }
  out=list("maximum"=ismax,"minimum"=ismin,"results"=a)  
  return(out)
}
