ecdfAUC <- function(dd, 
    col.line=NULL,
    main="ECDF",
    cex.leg=1,
    drawRef=FALSE,
    rm1=FALSE, 
    lineType=NULL,
    addLeg=TRUE,
    xlab="p-value",
    ylab="ECDF",
    cex.axis=1.5,
    cex.main=1.8, 
    cex.lab=1.2,
    axis.padj=c(-1, 1),
    lab.padj=c(-1.5, 1),
    lwd=1,
    box.lwd=1.2
){
    require(caTools)
    dd[is.na(dd)]=1
    nr=nrow(dd)
    nc=ncol(dd)
    colNm=colnames(dd)
    no1<-function(x){
        x=x[x<1]
    }
    
    if(is.null(col.line)) 
        col.line=1:nc
    if(is.null(lineType)) { 
        lineType=2:(nc+1) 
    } else {
        lineType=lineType[1:nc]
    }
    auc=rep(0,nc)
    j=1
    if (rm1) {
        yy=no1(dd[,j])
    } else {
        yy=dd[,j]
    }
    od=order(yy)
    yy=yy[od]    
    p1=ecdf(yy)
    auc[j]=1-trapz(p1(yy), yy)

    if(is.null(lwd)){
        lwd1=c(1, rep(2, ncol(dd)-1))
    } else  {
        lwd1=lwd
    }
    
    plot(p1(yy)~yy, main=main, lty=lineType[1], lwd=lwd1-0.5, col=col.line[j], 
        type="l", ylab="", xlab="", ylim=c(-0.05,1.05), xlim=c(-0.05,1.05), 
        cex.main=cex.main, axes=FALSE)
    mtext(text=ylab, side=2, padj=lab.padj[1], cex=cex.lab)
    mtext(text=xlab, side=1, padj=lab.padj[2], cex=cex.lab)
    axis(1, padj=axis.padj[1], cex.axis=cex.axis, lwd=box.lwd-0.01) 
    axis(2, padj=axis.padj[2], cex.axis=cex.axis, lwd=box.lwd-0.01)
    box(lwd=box.lwd)
    
    if(nc>1){
        for(j in 2:nc){
            if (rm1)
                yy=no1(dd[,j])
            else 
                yy=dd[,j]
            od=order(yy)
            yy=yy[od]    
            p1=ecdf(yy)
            auc[j]=1-trapz(p1(yy), yy)
            lines(p1(yy)~yy,  lty=lineType[j], col=col.line[j], lwd=lwd1)      
        }
    }  
    
    if(drawRef)
        abline(a=0, b=1, col="gray")
    
    leg1=colNm    
    dd1=dd
    colnames(dd1)=leg1
    if(addLeg){
        legend("bottomright",legend=leg1, col=col.line, lty=lineType, lwd=lwd1, 
            cex=cex.leg)
    }
    
    names(auc)=colNm
    return(auc)
}
