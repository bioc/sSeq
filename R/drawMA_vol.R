drawMA_vol <- function(                                                         
   y,                                   
   groups2,                           
   pv,                                
   cutoff=NULL,
   xlab1="(log2(A)+log2(B))/2",       
   ylab1="log2(A)-log2(B)",  
   tt1="MA plot",            
   tt2="volcano plot",       
   log2FoldChange =NULL,
   col1=c("black","red")
 ){
    AB=unique(groups2)
    A=AB[1]
    B=AB[2]
    if(length(pv)!=dim(y)[[1]])
        stop("check if pv is matched to y for drawMAfunction.")
    if(ncol(y)!=length(groups2) & length(unique(groups2))!=2) 
        stop("check y and groups2 for drawMA function.")
    muA=rowMeans(as.matrix(y[,groups2%in%A]))
    muA[muA==0]=1
    muB=rowMeans(as.matrix(y[,groups2%in%B]))
    muB[muB==0]=1
    if(is.null(log2FoldChange))
        log2FoldChange=log2(muA/muB)
    baseMean=log2(muA*muB)/2
   
    op<-par(mfrow=c(1,2))
    if( is.null(xlab1) |  is.null(ylab1) ){
   	    cond1 = unique(groups2)
   	    leg.A = cond1[1]
   	    leg.B = cond1[2]
   	    xlab1=paste("(log2(", leg.A, ")+log2(", leg.B, "))/2", sep='')
   	    ylab1=paste("(log2(", leg.A, ")-log2(", leg.B, "))/2", sep='')
    }
    plot(y=log2FoldChange, x=baseMean, pch=16, cex=0.3, col=col1[1],
        xlab=xlab1, ylab=ylab1, main=tt1)

    if(is.null(cutoff)){
        cutoff = min(pv[order(pv)][floor(nrow(y)*0.05)], 0.05)
        print(paste('cutoff', cutoff))
    }
    sel = pv<cutoff
    points(x=baseMean[sel], y=log2FoldChange[sel], col=col1[2], cex=0.6, pch=16)
    abline(h=0, col="blue")
   
    log2pv = -log2(pv)
    plot(y=log2pv, x=log2FoldChange,  pch=16, cex=0.3, col=col1[1],
         ylab="-log2(pvalue)", xlab="log2(fold change)", main=tt2 )
    points(y=log2pv[sel], x=log2FoldChange[sel], col=col1[2], cex=0.3, pch=16)          
    abline(v=0, col="blue")
 }
