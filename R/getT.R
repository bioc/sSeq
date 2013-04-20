getT <- function(
    countsTable, 
    sizeFactors=NULL, 
    q.vec=NULL, 
    plotASD=FALSE, 
    numPart=1, 
    propForSigma=c(0,1), 
    verbose=TRUE, 
    shrinkTarget=NULL, 
    shrinkQuantile=NULL, 
    shrinkVar=FALSE, 
    eSlope=0.05, 
    disp=NULL, 
    dispXX=NULL, 
    normalize=FALSE, 
    lwd1=4.5, 
    cexlab1=1.2
){
    if(!is.null(countsTable)) {counts=as.matrix(countsTable)}
    if(is.null(countsTable) & is.null(disp)){
        stop("at least provide the initial dispersion estimates.")
    }
    if(is.null(sizeFactors) & !is.null(countsTable)) { 
        sizeFactors=getNormFactor(countsTable)    	
    }
    
    if(is.null(eSlope)){
        eSlope=0.002
    }else{
        if(length(eSlope)>1 & verbose) 
            print("Note: only the first value in eSlope is used for tests.")
    }
		
    allAdjDisp=list()	
    if(is.null(disp) & !is.null(countsTable)){
        normc=as.matrix(t(t(counts)/sizeFactors))
        normc.m=rowMeans(normc)
        normc.v=rowVars(as.matrix(t(t(counts)/sqrt(sizeFactors))))
        s_st=mean(1/sizeFactors)
        disp=(normc.v - normc.m )/(normc.m)^2 
        normc.m[normc.m<=0]=1
        log.normc.m=log(normc.m)
    } else if (!is.null(dispXX)){
        normc.m=dispXX 
        normc.m[normc.m<=0]=1 
        log.normc.m=log(normc.m)
    } else {    
        normc.m=NULL 
        log.normc.m=NULL 
    }
    
    if(shrinkVar&is.null(disp)){
        disp=normc.v 
        if(verbose) 
            print("Shrinkage estimates on variance are used.")
    } else {  
        if(verbose)	
            print("Shrinkage estimates on dispersion are used for the tests.")
    }
    disp[is.na(disp) ]=0
    disp[disp<=0]=0  
    
    if(numPart==1){  
        disp.m= mean(disp)
        asd.mle=round(mean((disp-disp.m)^2, na.rm=T), 4)
        rg.xx=quantile(disp[is.finite(disp)], prob=c(0.05,0.995))
        xx=seq(rg.xx[1], rg.xx[2], length.out=200)
        asd=rep(0, length(xx))    
        for(i in 1:length(xx)){
            allAdjDisp[[i]]=equalSpace(disp, log.normc.m, 1, 
                propForSigma=propForSigma,  shrinkTarget=xx[i], vb=FALSE) 
            allAdjDisp[[i]]=pmax(1e-8, allAdjDisp[[i]])  
            names(allAdjDisp[[i]])=1:length(disp)
            asd[i]=mean((allAdjDisp[[i]]- disp)^2, na.rm=T)
        }
        diff.q=diff.asd=rep(0, length(asd))
        maxASD=max(asd, na.rm=T)
        maxASD.pnt=which( asd == maxASD)
        for ( i in 1:length(asd) ){ 
            diff.asd[i]=maxASD - asd[i]
            diff.q[i]=xx[maxASD.pnt] - xx[i] 
        }
        numAdjPoints=6
        len.asd=length(asd) - numAdjPoints + 1
        slope1=rep(1, len.asd)
        if(normalize){
            xx1=xx/sd(xx)
            yy1=asd/sd(asd)   	
            eSlope=eSlope*5
        } else {
            xx1=xx
            yy1=asd
        }
        for (i in 1:len.asd){ 
            slope1.xx=xx1[ i:(i+numAdjPoints-1) ]
            slope1.yy=yy1[ i:(i+numAdjPoints-1) ]
            slope1[i]=cov(slope1.xx, slope1.yy)/var(slope1.xx)
        }
        maxSlope1=max(abs(slope1))
        maxSlope1.pnt=which(abs(slope1) == maxSlope1)
        sub.slope1=abs(slope1)[maxSlope1.pnt:len.asd]
        sub.diff.asd=diff.asd[maxSlope1.pnt:length(diff.asd)] 
        pred.diff=matrix(NA, nrow=length(sub.diff.asd), ncol=numAdjPoints)
        for(i in 1:length(sub.diff.asd)){
            for(j in 1:numAdjPoints ){ 
                if( i-j >= 0){
                    pred.diff[i,j]=sub.diff.asd[i]/sub.slope1[i-j+1]
                }
            }
        }
        max.pred= max(pred.diff, na.rm=T)
        max.rowInd=which(apply(pred.diff, 1, max, na.rm=T)==max.pred) 
        temp.max.pnt=max.rowInd+maxSlope1.pnt-1-ceiling(numAdjPoints/2)
        max.pnt=rep(0, length(eSlope))
        for(k in 1:length(eSlope)){
            max.pnt[k]=temp.max.pnt[1] 
            while(-slope1[max.pnt[k]-ceiling(numAdjPoints/2)]<eSlope[k]&
                slope1[max.pnt[k]-ceiling(numAdjPoints/2)]<0
            ){
                max.pnt[k]=max.pnt[k] - 1
            }
        }
        if(!is.null(shrinkQuantile)){
            max.pnt=1 
            q.vec=shrinkQuantile
            adjDisp1=equalSpace(disp, log.normc.m, numPart, 
                propForSigma=propForSigma, shrinkQuantile=q.vec[1], vb=FALSE) 
            adjDisp1=pmax(1e-8, adjDisp1)  
            names(adjDisp1)=1:length(disp)
            asd.target=mean((adjDisp1- disp)^2, na.rm=T)   
            target=round(quantile(disp, prob=q.vec[1]),3) 
        }
   
        if(!is.null(shrinkTarget)){
            max.pnt=1 
            disp.tm=c(disp[!is.na(disp)], shrinkTarget[1])
            q.vec=round(rank(disp.tm)[disp.tm==shrinkTarget[1]]/
                length(disp[!is.na(disp)]), 3)
            adjDisp1=equalSpace(disp, log.normc.m, numPart, 
                propForSigma=propForSigma, shrinkQuantile=q.vec[1], vb=FALSE) 
            adjDisp1=pmax(1e-8, adjDisp1)  
            names(adjDisp1)=1:length(disp)
            asd.target=mean((adjDisp1- disp)^2, na.rm=T)
            target=shrinkTarget[1]
        }

        if(is.null(shrinkQuantile) & is.null(shrinkTarget)){
            target=asd.target=q.vec=rep(0,length(eSlope))
            for(k in 1:length(eSlope)){
                target[k]=xx[max.pnt[k]][1]
                asd.target[k]=asd[max.pnt[k]]
                disp.tm=c(disp[!is.na(disp)], target[k])
                q.vec[k]=round(rank(disp.tm)[disp.tm==target[k]]/
                    length(disp[!is.na(disp)]), 3)
            }
        }
   
        if (plotASD) {
            tt1="ASD vs shrinkage target"
            plot(asd~xx, cex=0.5, type="l", lwd=lwd1, main="", ylab='', 
                xlab='', axes=FALSE)
            mtext(text=expression(paste("ASD(",xi, ")", sep='')), side=2, 
                padj=-0.9, cex=cexlab1)
            mtext(text=expression(xi), side=1, padj=1.3, cex=cexlab1)
            axis(1, padj=-1, cex.axis=1.2) 
            axis(2, padj=1, cex.axis=1.2)
            if(!is.null(shrinkTarget) | !is.null(shrinkQuantile)){
                abline(h=asd.target, lty=2, col="blue")
                abline(v=target, lty=2, col="blue")     	
                mtext(text=round(asd.target, 3), side=4, at=asd.target, las=2)
                mtext(text=round(target,3), side=3, at=target, las=0)
            } else {
                prev.max.pnt=0
                for(k in 1:length(eSlope)){
                    if(prev.max.pnt==max.pnt[k]) 
                        next       	
                    abline(h=asd.target[k], lty=2, col="gray20")
                    abline(v=target[k], lty=2, col="gray20")
                    mtext(text=round(asd.target[k], 3), side=4, 
                        at=asd.target[k], las=2)
                    if(k==1){
                        leg.k=as.expression( 
                            bquote(hat(italic(xi))~.(paste("=", 
                            round(target[k],3),sep='')) ))
                    } else {
                        leg.k=target[k]
                    }
                    mtext(text=leg.k, side=3, at=round(target[k],3), las=0)
                    prev.max.pnt=max.pnt[k]
                }
            }
            abline(h=asd.mle, lty=1, col="gray40")
            mtext(text=asd.mle, side=4, at=asd.mle, las=2)
        }
   
        if(verbose){
            if(!is.null(shrinkQuantile) | !is.null(shrinkTarget) ){
                print(paste("The selected shrink target is", target[1]))
                print(paste("The selected shrink quantile is", q.vec[1]))
            } else { 
                print(paste("The shrink target is", target[1])) 
                print(paste("The shrink quantile is", q.vec[1])) 
            }
        }
        return(list(q=q.vec[1], target=target[1]))
   
    } else if (numPart>1) {
        if(is.null(log.normc.m)){
            stop("Error in getT: log.normc.m can not be NULL.")
        }
        out=getTgroup(y=disp, x=log.normc.m, numPart=numPart, plotASD=plotASD,
            verbose=verbose, eSlope=eSlope, lwd1=lwd1, cexlab1=cexlab1)
        return(out)	
    } else {
        stop("Error: numPart must be a non-negative integer.")
    }
}
