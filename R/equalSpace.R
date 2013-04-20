equalSpace <- function(
    y, 
    x=NULL, 
    numcls=1, 
    propForSigma=c(0,1), 
    shrinkTarget=NULL, 
    shrinkQuantile=0.975, 
    vb=TRUE
){
    if(numcls==1 | is.null(x))
        return(getAdjustDisp(y, propForSigma= propForSigma, shrinkTarget, 
            shrinkQuantile, verbose=vb)$adj)
    if(!is.null(shrinkTarget) & length(shrinkTarget)!=numcls){
        print(paste("Warning: the number of shrink targes is unequal to the",
            "number of pre-decied groups. Only the first target is used."))
        shrinkTarget=shrinkTarget[1]
        numcls=1
    }
    
    if(sum(is.na(x))>0 ) 
        print("The NA values in the dependent variable were ignored.")
    if( length(y) != length(x) ) 
        stop(paste("Error: check the input of equalSpace. y and x have", 
            "unequal lengths in equalSpace function."))
    rgx=range(x[x>-Inf]) 
    cut=seq(from=rgx[1], to=rgx[2], length=numcls+1)
    cls=rep(1, length(y))
    cls[x<=cut[2]]=1
    cls[x>cut[numcls]]=numcls
    for(i in 2:(numcls-1)){
        cls[x>cut[i] & x<=cut[i+1]]=i
    }
    sizes=tapply(rep(1, length(cls)), cls, sum)
    js=y
    mean.y=mean(y)
    for(i in 1:length(sizes)){
        if(sizes[i]>2){
            x.sub=x[cls==i]
            if(!is.null(shrinkTarget)){
                mixr=getAdjustDisp(y[cls==i], propForSigma=propForSigma, 
                    shrinkTarget[i], shrinkQuantile, verbose=vb)
            } else {
                mixr=getAdjustDisp(y[cls==i], propForSigma= propForSigma, 
                    shrinkTarget=NULL, shrinkQuantile=shrinkQuantile, 
                    verbose=vb)
            }
            js[cls==i]=mixr$adj 
        } else {
            js[cls==i]=mean.y
        }
    }
    return(js)
}
