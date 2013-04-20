sim <- function(
    ngenes, 
    true_mean1,
    conds, 
    alpha=function(m){rep(0.1, length(m))}, 
    mean_DE=0, 
    sd_DE=2,  
    s0=NULL, 
    s0_mean=2, 
    s0_sd=3, 
    true_isDE_proportion =0.3 
){
    if(is.null(s0)){
        s0=round(exp(rnorm(length(conds), mean=s0_mean, sd=s0_sd)),3) 	
    }
    true_isDE=(1:ngenes) %in% 
        sample(1:ngenes, size=ceiling(true_isDE_proportion*ngenes))
    true_randChg =(rnorm( ngenes, mean=mean_DE, sd=sd_DE ))*true_isDE
    nms =paste(1:ngenes, true_isDE, sep="_")
    true_mean2 =data.frame(A=true_mean1/exp(true_randChg), 
        B=true_mean1*exp(true_randChg)) 
    true_mean2[true_mean2<=0] =1
    rownames(true_mean2) =nms 
    true_dispA=alpha(true_mean2$A)
    true_dispB=alpha(true_mean2$B)
    true_disp=data.frame(A=true_dispA, B=true_dispB)
    true_v =data.frame(A=true_mean2$A+true_mean2$A^2*true_dispA, 
          B=true_mean2$B+true_mean2$B^2 * true_dispB ) 
    rownames(true_v) =nms
    countsTable =matrix(0, nrow=ngenes, ncol=length(conds))
    for (i in 1:length(conds)){
        q=true_mean2[[ conds[i] ]] * s0[i]
        d=true_disp[[conds[i] ]]/s0[i]
        for(j in 1:nrow(countsTable)){
            countsTable[j,i] =rnbinomMV( 1, mu=q[j], v=q[j]+q[j]*d[j] )
        }
    }
    countsTable[is.na(countsTable)]=0
    sum( rowSums(countsTable)==0 ) 
    nonzero =rowSums(countsTable)>0    
    rownames(countsTable)=nms
    colnames(countsTable)=paste(conds, 1:length(conds), sep="")
    return(list(countsTable=countsTable[nonzero,], 
        true_isDE=true_isDE[nonzero], true_mean=true_mean2[nonzero,], 
        true_var=true_v[nonzero,], true_disp=true_disp[nonzero,], s=s0))
}
