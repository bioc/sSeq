exactNBtest1 <- function(
    kA,
    kB,
    mu,
    disp,
    sA=1,
    sB=1,
    rA=0.5,
    rB=0.5
){
    ks=kA + kB
    allps1=dnbinom(0:ks, mu=sA*mu*rA, size=sA/disp)*
        dnbinom(ks:0, mu=sB*mu*rB, size=sB/disp)    
    pobs1=dnbinom(kA, mu=sA*mu*rA, size=sA/disp)*
        dnbinom(kB, mu=sB*mu*rB, size=sB/disp)                 
    sumsel=sum(allps1[allps1<=pobs1], na.rm=T)
    sumall=sum(allps1, na.rm=T)
    pval=sumsel/sumall  
    pval
}
