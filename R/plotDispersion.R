plotDispersion <- function(
    DispSH,
    extraOutput=NULL, 
    plotMethod="logDisp", 
    ylim1=NULL,      
    legPos="topleft", 
    myCol=brewer.pal(9, "Set1"), 
    tt=NULL 
){
    sumDisps.raw=DispSH$raw
    sumDisps.mixSH=DispSH$SH
    mus=DispSH$mus
    if(is.list(extraOutput)){
        extra.names=names(extraOutput)
        extra.num=length(extraOutput)
    } else {
        extra.num=0 
        extra.names=""
    }

    if(plotMethod=="logDisp"){
        log.sum=log(sumDisps.raw) 
        log.sum[is.na(log.sum)]=0
        if(is.null(ylim1)) {
          com.disp=c(log.sum, log(c(sumDisps.mixSH, sumDisps.raw)) ) 
          com.disp=com.disp[is.finite(com.disp)]
          rgy=range( com.disp, na.rm=T)
        }  else {   rgy=ylim1 }
        if(is.null(tt)){tt=""}
        smoothScatter( log.sum~log(mus),  cex=0.9,  ylim=rgy,  
            ylab="log(dispersions)", xlab="log(means)",  main=tt,  
            col= blues9[5], cex.axis=1.5, cex.lab=1.5 )    
        if(extra.num>0){
            extra.col=rep("", extra.num)
            for(i in 1:extra.num){
                if(length(extraOutput[[i]]) != length(mus)){
                    print("input error for", extra.names[i]) 
                    next
                }
                points(log(extraOutput[[i]])~log(mus), 
                    col=ifelse(i==3, 1, myCol[i]), pch=1, cex=1)
                extra.col[i]=ifelse(i==3, 1, myCol[i])
            }
            leg1=c("sSeq", names(extraOutput))
            col1=c(c(myCol[3]), extra.col)
        } else {
            leg1=c(expression(paste("log(", hat(phi)[g]^"sSeq", ")", 
                sep='')), expression(paste("log(", hat(phi)[g]^"MM", ")",
                sep='')))
            col1=c("black", blues9[5])
            legPos="bottomright"
        }
        points(log(sumDisps.mixSH)~log(mus), col=col1[1], pch=16)
        legend(legPos, legend=leg1, col=col1, pch=16, bty="n", cex=1.5)
    }
    
    if(plotMethod=="disp"){
        if(is.null(ylim1)) 
            rgy=range(c(sumDisps.raw, sumDisps.mixSH), na.rm=T)
        else 
            rgy=ylim1
        plot(sumDisps.raw~log(mus), cex=0.9, ylim=rgy, 
            ylab="dispersions", xlab="log(means)", 
            main=ifelse(is.null(tt), "dispersions"))
        if(extra.num>0){
            extra.col=rep("", extra.num)
            for(i in 1:extra.num){
                if(length(extraOutput[[i]]) != length(mus)){
                    print("input error for", extra.names[i]) 
                    next
                }
                points( (extraOutput[[i]])~log(mus), 
                    col=ifelse(i==3, 1, myCol[i]), pch=1, cex=1)
                extra.col[i]=ifelse(i==3, 1, myCol[i])
            }
            leg1=c("sSeq", names(extraOutput))
            col1=c(c(myCol[3]), extra.col)
        } else {
            leg1="sSeq"
            col1=myCol[3]
        }
        points( (sumDisps.mixSH)~log(mus), col=myCol[3], pch=16, cex=0.2)
        legend(legPos, legend=leg1, col=col1, pch=c(1), bg="white")                     
    }
}
