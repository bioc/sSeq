nbinomTestForMatricesSH <- function(
    countsA, 
    countsB, 
    sizeFactorsA, 
    sizeFactorsB, 
    numPart=1, 
    SHonly=FALSE, 
    propForSigma=c(0,1), 
    shrinkTarget=NULL, 
    shrinkQuantile=NULL, 
    cLA, 
    cLB, 
    contrast=NULL, 
    keepLevelsConsistant=TRUE,
    useMMdisp=FALSE,
    shrinkVariance=FALSE,
    pairedDesign=FALSE, 
    pairedDesign.dispMethod="per-pair", 
    useFisher=FALSE,
    Dispersions=NULL,
    eSlope=NULL,
    plotASD=FALSE,
    lwd_ASD=4.5, 
    cex_ASD=1.2
){
    cl.nm=sort(unique(c(cLA, cLB))) 
    cntA=as.matrix(countsA) 
    cntB=as.matrix(countsB) 
    sfA=sizeFactorsA 
    sfB=sizeFactorsB 
    s.m=1/mean(1/c(sfA, sfB))  
    s_st=mean(1/c(sfA, sfB))    
    s.tilde1=sum(c(sfA, sfB)) 
    normA=t(t(cntA)/sfA) 
    normB=t(t(cntB)/sfB)
    norm=cbind(normA, normB)
    norm.m=rowMeans(norm)
    norm.v=rowVars(norm)
    norm.mA=rowMeans(normA)
    norm.vA=rowVars(normA)
    norm.mB=rowMeans(normB) 
    norm.vB=rowVars(normB)    
    orig.shrinkQuantile=shrinkQuantile
    orig.shrinkTarget=shrinkTarget
    adj.getT=FALSE
    ll=list()    
    for(L in 1:length(cl.nm)) {
        countsA=as.matrix(cntA[,cLA==cl.nm[L]]) 
        nA=ncol(countsA)
        countsB=as.matrix(cntB[,cLB==cl.nm[L]]) 
        nB=ncol(countsB)
        if( is.null(orig.shrinkQuantile) & is.null(shrinkTarget) & 
            is.null(Dispersions)
        ){
            if(length(cl.nm)>1) 
                print(paste("Get shrinkage target at level", cl.nm[L]))
            countsTable1=cbind(countsA, countsB)
            sf1=getNormFactor(countsTable1)
            out.getT=getT(countsTable1, sf1, numPart=numPart, plotASD=FALSE, 
                propForSigma=propForSigma, verbose=FALSE, 
                shrinkVar=shrinkVariance, eSlope=eSlope, disp=NULL, 
                lwd1=lwd_ASD, cexlab1=cex_ASD) 
            shrinkTarget=out.getT$target
            if(min(out.getT$target)>=0.8 & min(out.getT$q)>=0.965){
                adj.getT=TRUE
            }
        }
        N=nA + nB 
        rA=nA/N  
        rB=nB/N
        cnt=cbind(countsA, countsB) 
        ng=nrow(cnt)
        kAs=rowSums(cbind(countsA)) 
        kBs=rowSums(cbind(countsB))
        conds0=c(rep('A',nA),rep('B',nB))
        countsTable0=cbind(countsA, countsB)
        sizeFactors=getNormFactor(countsTable0)    
        sizeFactorsA=sizeFactors[conds0=="A"]
        sizeFactorsB=sizeFactors[conds0=="B"]
        normcA=t(t(countsA)/sizeFactorsA)
        normcB=t(t(countsB)/sizeFactorsB)   
        normc=as.matrix(cbind(normcA, normcB), ncol=nA+nB)
        normc.m=rowMeans(normc)
        normc.v=rowVars(normc)
        s=sum(sizeFactors)
        s.mL=1/mean(1/sizeFactors)
        s.tilde=sum(sizeFactors)
        dispc=s.mL*(normc.v - normc.m*s_st )/(normc.m)^2
        dispc[is.na(dispc) ]=0
        dispc[dispc<=0]=0 
        xx=normc.m
        xx[xx<=0]=1
        if(shrinkVariance){
            normc.v[is.na(normc.v)]=0
            normc.v=equalSpace(normc.v, log(xx), numPart, 
                propForSigma=propForSigma, shrinkTarget=shrinkTarget, 
                shrinkQuantile=shrinkQuantile, vb=FALSE)
            adjdispc=s.m*(normc.v-normc.m*s_st )/(normc.m)^2
            adjdispc[is.na(adjdispc) ]=0 
            adjdispc[adjdispc<=0]=0
        } else {
            if(is.null(Dispersions)){
                adjdispc=equalSpace(dispc, log(xx), numPart, 
                    propForSigma=propForSigma, shrinkTarget=shrinkTarget, 
                    shrinkQuantile=shrinkQuantile, vb=FALSE) 
            } else {
                adjdispc=Dispersions
            }
        }
        sumDispsc=adjdispc 
        names(sumDispsc)=rownames(countsA)
        ll[[L]]=list(countsA=countsA, countsB=countsB, 
                    normcA=normcA, normcB=normcB,
                    sumDisps=adjdispc, disp=dispc, sfA=sizeFactorsA, 
                    sfB=sizeFactorsB, mus=normc.m, s=s, 
                    sA=sum(sizeFactorsA), sB=sum(sizeFactorsB),
                    s.tilde=s.tilde, nA=nA, nB=nB, s.m=s.mL, 
                    counts=cnt, kAs=kAs, kBs=kBs,  rA=rA, rB=rB)
    }
    s.m1=1
    for(L in 1:length(cl.nm)){ 
        s.m1=max(ll[[L]]$s.m, s.m1) 
    }
    shrinkQuantile=orig.shrinkQuantile
    shrinkTarget=orig.shrinkTarget    
    verbose=TRUE
    if(pairedDesign){
        if(pairedDesign.dispMethod=="per-pair"){
            print(paste("For paired design, the per-pair dispersion", 
                "estimates are used and shrunk separately."))
            verbose=FALSE
            disp=sumDisps=NULL
            for(L in 1:length(cl.nm)){
                disp=cbind(disp, ll[[L]]$disp)
                sumDisps=cbind(sumDisps, ll[[L]]$sumDisps)
            }
        } else if (pairedDesign.dispMethod=="pooled"){
            print(paste("For paired design, the aveaged dispersion", 
                "estimates across paires are used."))
            disp.pair=NULL
            for(L in 1:length(cl.nm)){
                disp.pair=cbind(disp.pair, ll[[L]]$disp)
            }
            disp=rowMeans(disp.pair, na.rm=T)
        } else {
            stop(paste("Error in defining the method of dispersion", 
                "estimates for paired design. \nPlease input 'per-pair'", 
                "or 'pooled'."))
        }
    } else {
        disp=s.m*(norm.v-norm.m*s_st)/(norm.m)^2
        disp[ is.na(disp) ]=0  
        disp[disp<=0]=0 
    }
    if(shrinkVariance){
        xx=norm.m 
        xx[xx<=0]=1    
        norm.v[is.na(norm.v)]=0  
        if(is.null(shrinkQuantile) & is.null(shrinkTarget)){
            shrinkTarget=getT(countsTable=NULL, sizeFactors=NULL, 
                numPart=numPart, propForSigma=propForSigma, verbose=FALSE,
                plotASD=plotASD, shrinkVar=shrinkVariance, eSlope=eSlope, 
                disp=norm.v, dispXX=norm.m, lwd1=lwd_ASD, 
                cexlab1=cex_ASD)$target
        }
        norm.v=equalSpace(norm.v, log(xx), numPart, 
            propForSigma=propForSigma, shrinkTarget=shrinkTarget, 
            shrinkQuantile=shrinkQuantile, vb=FALSE)    
        adjdisp=s.m*( norm.v - norm.m*s_st )/(norm.m)^2
        adjdisp[ is.na(adjdisp) ]=0  
        adjdisp[adjdisp<=0]=0
        sumDisps=adjdisp
    } else {
        if(useMMdisp){
            sumDisps=disp
            adjdisp=rep(NA, length(sumDisps))
        } else {  
            if(is.null(shrinkQuantile) & is.null(shrinkTarget) & 
                is.null(Dispersions)
            ){
                out.getT=getT(countsTable=NULL, sizeFactors=NULL, 
                    numPart=numPart, plotASD=plotASD, 
                    propForSigma=propForSigma, verbose=verbose, 
                    shrinkVar=shrinkVariance, eSlope=eSlope, disp=disp, 
                    dispXX=norm.m, lwd1=lwd_ASD, cexlab1=cex_ASD)
                shrinkTarget=out.getT$target
                if(min(out.getT$target)>=0.8 & min(out.getT$q)>=0.965){
                    adj.getT=TRUE
                }
            }
            if(!pairedDesign | 
                (pairedDesign & pairedDesign.dispMethod=="pooled")
            ){ 
                xx=norm.m
                xx[xx<=0]=1      
                if(is.null(Dispersions)){
                    adjdisp=equalSpace(disp, log(xx), numPart, 
                        propForSigma=propForSigma, vb=FALSE, 
                        shrinkTarget=shrinkTarget, 
                        shrinkQuantile=shrinkQuantile)     
                } else {
                    adjdisp=Dispersions
                }
                sumDisps=adjdisp
            }
        }
    }
    if(length(ll)>1){ 
        temCount=ll[[1]]$mus
        for(L in 2:length(ll)){
            temCount=cbind(temCount, ll[[L]]$mus)
        }
        lm=rowMeans(log(temCount))
        int.func=function(x1,lm){ 
            exp(median((log(x1)-lm)[is.finite(lm)]))
        }
        s_L=apply(temCount, 2, int.func, lm)
        s_L1 =sum(s_L)
        s_L.m=s_L1/length(s_L)
    } else {
        s_L=1 
        s_L.m=1
    }
    if(SHonly) {
        return(data.frame(SH=sumDisps, raw=disp, mus=norm.m))
    }
    if(!is.null(Dispersions)){
        if(length(Dispersions)!=length(norm.m)){
            stop(paste("Please let the length of input Dispersion", 
                "equal to the number of rows in countsTable, or", 
                "set it as NULL."))
        }
        print("The known dispersion values are used.")
        sumDisps=Dispersions
    }    
    perc=c(1,3,5,7,9, 10)
    progress=round(perc/10 * ng, 0) 
    progress.id=1
    pval=rep(1, ng)
    if(pairedDesign){
        for(i in 1:ng){
            if(progress.id<length(progress)&i==progress[progress.id]){
                progress.id=progress.id + 1
                print(paste(perc[progress.id]*10, "% processed.", sep=''))
            }
            log.p=0
            for(L in 1:length(cl.nm)){
                kA1=ll[[L]]$kAs[i] 
                kB1=ll[[L]]$kBs[i] 
                norm.mA1=s_L.m*norm.m[i]*ll[[L]]$sA 
                norm.mB1=s_L.m*norm.m[i]*ll[[L]]$sB 
                sA1=s_L.m*ll[[L]]$s.tilde
                sB1=s_L.m*ll[[L]]$s.tilde
                ks=kA1 + kB1       
                if(pairedDesign.dispMethod=="per-pair"){
                    allps1=dnbinom( 0:ks, mu= norm.mA1,  
                        size=sA1/ll[[L]]$sumDisps[i])*dnbinom(ks:0, 
                        mu=norm.mB1, size=sB1/ll[[L]]$sumDisps[i] ) 
                    pobs1=dnbinom(kA1, mu=norm.mA1, 
                        size=sA1/ll[[L]]$sumDisps[i])*dnbinom(kB1,
                        mu=norm.mB1, size=sB1/ll[[L]]$sumDisps[i] )      
                } else {
                    allps1=dnbinom( 0:ks, mu= norm.mA1,  
                        size=sA1/sumDisps[i])*dnbinom(ks:0, 
                        mu=norm.mB1, size=sB1/sumDisps[i] ) 
                    pobs1=dnbinom(kA1, mu=norm.mA1, 
                        size=sA1/sumDisps[i])*dnbinom(kB1, 
                        mu=norm.mB1, size=sB1/sumDisps[i] )      
                }
                sumsel=sum(allps1[allps1<=pobs1], na.rm=T)
                sumall=sum(allps1, na.rm=T) 
                log.p1=log(min(1, sumsel/sumall))
                if(is.finite(log.p1)) {
                    log.p=log.p + log.p1
                } 
            }
            if(useFisher){
                chi1=-2 * log.p
                pval[i]=1-pchisq(chi1, length(cl.nm)*2)        
            } else {
                pval[i]=exp(log.p)
            }
        } 
        pval[is.na(pval)]=1 
        return(list(pval=pval, dispMM=disp, dispSH=sumDisps, mu=norm.m))        
    }  
    
    if(is.null(contrast)){    
        for(i in 1:ng){
            if(progress.id<length(progress) & 
                i==progress[progress.id] 
            ){
                progress.id=progress.id + 1
                print(paste(perc[progress.id]*10,"% processed.", sep=''))
            }        
            sumsel=sumall=0  
            sign.diff=0 
            for(L in 1:length(cl.nm)){
                ks=(ll[[L]]$kAs[i] + ll[[L]]$kBs[i])
                nks=length(0:ks) 
                s=ll[[L]]$s 
                rA=ll[[L]]$rA 
                rB=ll[[L]]$rB     
                if(ll[[L]]$nA>1 & ll[[L]]$nB>1){    
                    norm.mA1=s_L.m*norm.m[i]*ll[[L]]$sA
                    norm.mB1=s_L.m*norm.m[i]*ll[[L]]$sB
                    sA1=s_L.m*ll[[L]]$sA*2
                    sB1=s_L.m*ll[[L]]$sB*2
                } else {
                    norm.mA1=s_L.m*norm.m[i]*ll[[L]]$s*ll[[L]]$rA
                    norm.mB1=s_L.m*norm.m[i]*ll[[L]]$s*ll[[L]]$rB              
                    sA1=s_L.m*ll[[L]]$s.tilde
                    sB1=s_L.m*ll[[L]]$s.tilde
                }
                allps1=dnbinom(0:ks, mu=norm.mA1, size=sA1/sumDisps[i])*
                    dnbinom(ks:0, mu=norm.mB1, size=sB1/sumDisps[i] )  
                pobs1=dnbinom(ll[[L]]$kAs[i], mu=norm.mA1, 
                    size=sA1/sumDisps[i])*dnbinom(ll[[L]]$kBs[i], 
                    mu=norm.mB1, size=sB1/sumDisps[i]) 
                sumsel=sumsel+ sum(allps1[allps1<=pobs1], na.rm=T)
                sumall=sumall+ sum(allps1, na.rm=T)
                sign.diff=sign.diff + 
                    (mean(ll[[L]]$normcA[i])<mean(ll[[L]]$normcB[i]))
            }
            if(keepLevelsConsistant & sign.diff>0&
                sign.diff<length(cl.nm)
            ){
                pval[i]=1
            } else {
                pval[i]=min(1, sumsel/sumall)        
            }
        }
        pval[is.na(pval)]=1
        if(adj.getT){
            rk=rank(pval, ties.method="min")
            wt=rk/max(rk)
            pval=pval*wt
        }       
    } else { 
        for(i in 1:ng){
            if (progress.id<length(progress)&i==progress[progress.id]){
                progress.id=progress.id + 1
                print(paste(perc[progress.id]*10,"% processed.", sep=''))
            }
            norm.mA1 =norm.mB1=kA1=kB1=sA1.sz=sB1.sz=
                sA1=sB1=sumDisps1= 0
            for(L in 1:length(cl.nm)){
                kA1=ll[[L]]$kAs[i] + kA1
                kB1=ll[[L]]$kBs[i] + kB1
                norm.mA1=ll[[L]]$mus[i]*ll[[L]]$sA+norm.mA1 
                norm.mB1=ll[[L]]$mus[i]*ll[[L]]$sB+norm.mB1 
                sA1=s_L.m*ll[[L]]$s.tilde + sA1
                sB1=s_L.m*ll[[L]]$s.tilde + sB1
            }
            sumDisps1=sumDisps[i] 
            ks= kA1 + kB1
            nks=length(0:ks)         
            allps1=dnbinom(0:ks, mu=norm.mA1, size=sA1/sumDisps1)*
                dnbinom(ks:0, mu=norm.mB1, size=sB1/sumDisps1)   
            pobs1=dnbinom(kA1, mu=norm.mA1, size=sA1/sumDisps1)*
                dnbinom(kB1, mu=norm.mB1, size=sB1/sumDisps1)  
            sumsel=sum(allps1[allps1<=pobs1], na.rm=T)
            sumall=sum(allps1, na.rm=T)
            pval[i]=min(1, sumsel/sumall)
        }
    }
    pval[is.na(pval)]=1    
    return(list(pval=pval, dispMM=disp, dispSH=sumDisps, mu=norm.m))
}
