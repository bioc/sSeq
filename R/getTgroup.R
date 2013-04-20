getTgroup <- function(
    y, 
    x, 
    numPart=10, 
    plotASD=FALSE, 
    verbose=FALSE, 
    eSlope=0.05, 
    lwd1=4.5, 
    cexlab1=1.2
){
    rgx=range(x[is.finite(x)]) 
    cut=seq(from=rgx[1], to=rgx[2], length=numPart+1)
    cls=rep(1, length(y)) 
    cls[x<=cut[2]]=1
    cls[x>cut[numPart]]=numPart
    for(i in 2:(numPart-1)){
        cls[x>cut[i]&x<=cut[i+1]]=i
    }
    sizes=tapply(rep(1, length(cls)), cls, sum)
    qall.vec=targetall=rep(1.0, numPart)
    for(gp in 1:numPart){
        allAdjy=list()
        y1=y[cls==gp]
        x1=x[cls==gp]
        y.m=mean(y1)
        asd.mle=round(mean((y1-y.m)^2, na.rm=T), 4)
        rg.xx=quantile(y[is.finite(y1)], prob=c(0.05,0.995))
        xx=seq(rg.xx[1], rg.xx[2], length.out=200)
        asd=rep(0, length(xx))    
        for(i in 1:length(xx)){
            allAdjy[[i]]=equalSpace(y=y1, x=x1, numcls=1, 
                shrinkTarget=xx[i], vb=FALSE) 
            allAdjy[[i]]=pmax(1e-8, allAdjy[[i]])  
            names(allAdjy[[i]])=1:length(y1)
            asd[i]=mean((allAdjy[[i]]- y1)^2, na.rm=T)
        }
        diff.q=diff.asd=rep(0, length(asd))
        maxASD=max(asd, na.rm=T)
        maxASD.pnt=which( asd == maxASD)
        maxASD.pnt=max(maxASD.pnt)
        for ( i in 1:length(asd) ){ 
            diff.asd[i]=maxASD - asd[i]
            diff.q[i]=xx[maxASD.pnt] - xx[i] 
        }
        numAdjPoints=6
        len.asd=length(asd) - numAdjPoints + 1
        slope1=rep(1, len.asd)
        xx1=xx
        y11=asd
        for (i in 1:len.asd){ 
            slope1.xx=xx1[ i:(i+numAdjPoints-1) ]
            slope1.y1=y11[ i:(i+numAdjPoints-1) ]
            slope1[i]=cov(slope1.xx, slope1.y1)/var(slope1.xx)
        }
        maxSlope1=max(abs(slope1))
        maxSlope1.pnt=which(abs(slope1)==maxSlope1)
        sub.slope1=abs(slope1)[maxSlope1.pnt : len.asd]
        sub.diff.asd=diff.asd[maxSlope1.pnt : length(diff.asd)]
        pred.diff=matrix(NA, nrow=length( sub.diff.asd), ncol=numAdjPoints)
        for(i in 1:length(sub.diff.asd)){
            for(j in 1:numAdjPoints ){ 
                if( i-j >= 0){
                    pred.diff[i,j]=sub.diff.asd[i]/sub.slope1[i-j+1]
                }
            }
        }
        max.pred=max(pred.diff, na.rm=T) 
        max.rowInd=which( apply(pred.diff, 1, max, na.rm=T) == max.pred ) 
      
        temp.max.pnt=max.rowInd + maxSlope1.pnt - 1 -ceiling(numAdjPoints/2) 
        max.pnt=rep(0,length(eSlope)) 
        for(k in 1:length(eSlope)){
            max.pnt[k]=temp.max.pnt[1] 
            tm1=-slope1[max.pnt[k]-ceiling(numAdjPoints/2)]
            tm2= -tm1
            while( !is.na(tm1) & tm1[1] < eSlope[k]  & tm2[1]<0){
                max.pnt[k]=max.pnt[k] - 1
                tm1=-slope1[max.pnt[k]-ceiling(numAdjPoints/2)]
                tm2= -tm1     	
            if(length(tm1)==0)       
                break
            }
        }    

        target=asd.target=q.vec=rep(0,length(eSlope))
        for(k in 1:length(eSlope)){
            target[k]=xx[max.pnt[k]][1]
            asd.target[k]=asd[max.pnt[k]][1]
            y.tm=c(y1[!is.na(y1)], target[k])
            q.vec[k]=round(rank(y.tm)[y.tm==target[k]]/
                length(y1[!is.na(y1)]), 3)
        }

        if (plotASD) {
            tt1=paste("ASD vs shrinkage target in group", gp)
            plot(asd~xx, cex=0.5, type="l", lwd=lwd1, main=tt1,
                ylab='', xlab='', axes=FALSE)
            mtext(text=expression(paste("ASD(",xi, ")", sep='')), side=2, 
                padj=-0.9, cex=cexlab1)
            mtext(text=expression(xi), side=1, padj=1.3, cex=cexlab1)
            box(lwd=1.8)
            axis(1, padj=-1, cex.axis=1.2) 
            axis(2, padj=1, cex.axis=1.2)
            prev.max.pnt=0
            for(k in 1:length(eSlope)){
                if(prev.max.pnt==max.pnt[k]) 
                    next 	
                abline(h=asd.target[k],  lty=2, col="gray20")
                abline(v=target[k], lty=2, col="gray20")
                mtext(text=round(asd.target[k], 3), side=4, 
                    at=asd.target[k], las=2)
                if(k==1){
                    leg.k=as.expression( bquote(hat(italic(xi))~.(paste("=", 
                        round(target[k],3),sep='')) ))
                } else {
                    leg.k=target[k]
                }
                mtext(text=leg.k, side=3, at=round(target[k],3), las=0)
                prev.max.pnt=max.pnt[k]
            }
            abline(h=asd.mle, lty=1, col="gray40")
            mtext(text=asd.mle, side=4, at=asd.mle, las=2)
        }
   
        if(verbose){
            print(paste("In group",gp, "the average of the values on", 
                "X-axis for", sizes[gp], "genes is", mean(x1, na.rm=T)))
            print(paste("shrinkTarget ", target[1], " and shrinkQuantile ",
                q.vec[1], ".", sep='')) 
        }

        qall.vec[gp]=q.vec[1]
        targetall[gp]=target[1]
    }

    return(list(q=qall.vec, target=targetall))
}	