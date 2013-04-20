getQ <- function(
    countsTable, 
    sizeFactors=NULL, 
    q.vec=NULL, 
    plotASD=FALSE, 
    numPart=1, 
    propForSigma=c(0,1), 
    verbose=TRUE, 
    shrinkTarget=NULL, 
    shrinkQuantile=NULL
){
    counts=as.matrix(countsTable)
    if(is.null(sizeFactors)) { 
        sizeFactors=getNormFactor(countsTable)         
    }
        
    if (is.null(q.vec)){
        q.vec=seq(0.05, 0.995, 0.005)
    }
    names(q.vec)=1:length(q.vec)
    allAdjDisp=list()
    
    normc=as.matrix(t(t(counts)/sizeFactors))
    normc.m=rowMeans(normc)
    normc.v=rowVars(normc)
    s_st=mean(1/sizeFactors)
    s.m=1/s_st
  
    disp=s.m*(normc.v - normc.m*s_st )/(normc.m)^2
    disp[is.na(disp) ]=0
    disp[disp<=0]=min(disp[disp>0])
    disp.m= mean(disp)
    asd.mle=round(mean((disp-disp.m)^2, na.rm=T), 4)
    xx=quantile(disp, prob=q.vec) #define the X-axis for ASD plot.    

    #calculate average squared difference
    normc.m[normc.m<=0]=1
    asd=rep(0, length(q.vec))
    for(i in 1:length(q.vec)){
        allAdjDisp[[i]]=equalSpace(disp, log(normc.m), numPart, 
            propForSigma= propForSigma,  shrinkQuantile=q.vec[i], vb=FALSE) 
        allAdjDisp[[i]]=pmax(1e-8, allAdjDisp[[i]])  
        names(allAdjDisp[[i]])=rownames(counts)
        asd[i]=mean((allAdjDisp[[i]]- disp)^2, na.rm=T)
    }      
   
    diff.q=diff.asd=rep(0, length(asd))
    maxASD=max(asd, na.rm=T)
    maxASD.pnt=which( asd == maxASD)
    for ( i in 1:length(asd) ){ #
        diff.asd[i]=maxASD - asd[i]
        diff.q[i]=q.vec[maxASD.pnt] - q.vec[i] 
    }

    numAdjPoints=4
    len.asd=length(asd) - numAdjPoints + 1
    slope1=rep(1, len.asd)
    for (i in 1:len.asd){ 
        slope1.xx=q.vec[ i:(i+numAdjPoints-1) ]
        slope1.yy=asd[ i:(i+numAdjPoints-1) ]
        slope1[i]=cov(slope1.xx, slope1.yy) / var(slope1.xx)
    }
    maxSlope1=max(abs(slope1)) 
    maxSlope1.pnt=which(abs(slope1) == maxSlope1) 
    sub.slope1=abs(slope1)[maxSlope1.pnt : len.asd] 
    sub.diff.asd=diff.asd[maxSlope1.pnt:length(diff.asd)] 
    pred.diff=matrix(NA, nrow=length( sub.diff.asd), ncol=numAdjPoints)
    for(i in 1:length(sub.diff.asd)){
        for(j in 1:numAdjPoints ){ 
            if( i-j >= 0){
                pred.diff[i,j]=sub.diff.asd[i] / sub.slope1[i-j+1]
            }
        }
    }
    max.pred= max(pred.diff, na.rm=T) 
    max.rowInd=which( apply(pred.diff, 1, max, na.rm=T) == max.pred ) 
   
    max.pnt=max.rowInd + maxSlope1.pnt - 1 
   
    if(!is.null(shrinkQuantile)){
        max.pnt=1 
        q.vec=shrinkQuantile
        adjDisp1=equalSpace(disp, log(normc.m), numPart, 
            propForSigma=propForSigma, shrinkQuantile=q.vec[1], vb=FALSE) 
        adjDisp1=pmax(1e-8, adjDisp1)  
        names(adjDisp1)=rownames(counts)
        asd.target=mean((adjDisp1- disp)^2, na.rm=T)         
    }
    if(!is.null(shrinkTarget)){
        max.pnt=1 
        disp.tm=c(disp[!is.na(disp)], shrinkTarget[1])
        q.vec=round(rank(disp.tm)[disp.tm==shrinkTarget[1]]/
            length(disp[!is.na(disp)]), 3)
        adjDisp1=equalSpace(disp, log(normc.m), numPart, 
            propForSigma=propForSigma, shrinkQuantile=q.vec[1], vb=FALSE) 
        adjDisp1=pmax(1e-8, adjDisp1)  
        names(adjDisp1)=rownames(counts)
        asd.target=mean((adjDisp1- disp)^2, na.rm=T)
    }

    if(!is.null(shrinkTarget)){
        target=shrinkTarget[1]
    } else {   
        target=round(quantile(disp, prob=q.vec[max.pnt[1]]),3) 
    }
   
    if (plotASD) {
        tt1="ASD vs shrinkage target"
        plot(asd~xx, cex=0.5, pch=16, main=tt1, 
            ylab="Average Sqared Difference", 
            xlab="the common information in dispersion estimates")    
        if(!is.null(shrinkTarget) | !is.null(shrinkQuantile)){
            abline(h=asd.target, lty=2, col="blue")
            abline(v=target, lty=2, col="blue")         
        } else {
            asd.target=asd[max.pnt[1]]
            abline(h=asd.target,  lty=2, col="blue")
            abline(v=target, lty=2, col="blue")                
        }
        abline(h=asd.mle, lty=1, col="red")
        mtext(text=asd.mle, side=4, at=asd.mle, las=2)
        mtext(text=round(asd.target, 3), side=4, at=asd.target, las=2)
        mtext(text=round(target,3), side=3, at=target, las=0)
    }
   
    if(verbose){
        if(length(q.vec)>1) {
            print(paste("The selected quantile target is ", 
                q.vec[max.pnt[1]], sep=''))   
        }else{
            if(!is.null(shrinkQuantile) | !is.null(shrinkTarget) ){
                print(paste("The selected quantile target is", q.vec))
            } else { 
                print(paste("The quantile target is", q.vec)) 
            }
        }
    }
    return(list(q=q.vec[max.pnt[1]], target=target))
}
