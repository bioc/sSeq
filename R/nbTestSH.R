nbTestSH <- function( 
    countsTable, 
    conds,
    condA="A", 
    condB="B",
    numPart=1, 
    SHonly=FALSE, 
    propForSigma=c(0, 1), 
    shrinkTarget=NULL, 
    shrinkQuantile=NULL, 
    plotASD=FALSE, 
    coLevels=NULL, 
    contrast=NULL, 
    keepLevelsConsistant=FALSE, 
    useMMdisp=FALSE, 
    addRawData=FALSE, 
    shrinkVariance=FALSE,
    pairedDesign=FALSE,
    pairedDesign.dispMethod="per-pair", 
    useFisher=FALSE,
    Dispersions=NULL,
    eSlope=0.05,
    lwd_ASD=4.5, 
    cex_ASD=1.2
){
    if(is.null(coLevels)){
        coLevels=data.frame(exp=rep(1,ncol(countsTable)))
    }
    for(j in 1:ncol(coLevels)){
        cL=apply(coLevels, 1, paste, collapse="_")
    }
    cL.tb=table(cL)       
    if( sum(cL.tb<2)>0 ){
       stop("Errors in 'coLevels'. All the levels must have paired comparisons.")
    }
    if(!is.null(contrast)){
        if(length(contrast)!=ncol(countsTable)){
            stop(paste("Error: the length of contrast vector must equal", 
                "to the number of columns in countsTable."))
        }
        if(length(unique(conds[contrast>0]))!=1){
            stop(paste("Error: this package is currently only available", 
                "for the contrast between conditions, not the contrast", 
                "within conditions. Please revise the contrast vector,", 
                "such as c(1,1,-1,-1) for cond=c('A','A', 'B, 'B'),", 
                "instead of c(1,-1,1,-1)."))
        }
        countsTable=countsTable[, contrast!=0]
        conds=conds[contrast!=0]
        cL=cL[contrast!=0]  
        contrast= contrast[contrast!=0]
        contrast.mat= t(matrix(contrast, nrow=length(contrast), 
            ncol=nrow(countsTable)))
        if(sum(abs(contrast)!=1)>0){
            countsTable=round(countsTable * abs(contrast.mat))
        }
    }
    sf=getNormFactor(countsTable)

    colA=conds == condA
    colB=conds == condB

    if(sum(colB[1:sum(colA)])>0){
        stop(paste("re-order the columns of the countsTable so that all", 
            "samples in the same condition are in the adjacent columns.", 
            "For example, 'A A B B' is good, but 'A B A B' is bad."))
    }
    
    counts=as.matrix(countsTable)
    ng=nrow(counts)
    cntA=matrix(counts[, colA], ncol=sum(colA)) 
    cntB=matrix(counts[, colB], ncol=sum(colB))
        
    if(SHonly){
        disp=nbinomTestForMatricesSH(countsA=cntA, countsB=cntB, 
            sizeFactorsA=sf[colA], sizeFactorsB=sf[colB],
            numPart=numPart, SHonly=TRUE, propForSigma= propForSigma, 
            shrinkTarget=shrinkTarget, shrinkQuantile=shrinkQuantile,
            cLA=cL[colA], cLB=cL[colB], keepLevelsConsistant, useMMdisp,
            shrinkVariance=shrinkVariance, eSlope=eSlope, plotASD=plotASD, 
            lwd_ASD=lwd_ASD, cex_ASD=cex_ASD)
        return(disp)
    } else {
        t1=Sys.time()
        pval0=nbinomTestForMatricesSH(countsA=cntA, countsB=cntB, 
            sizeFactorsA=sf[colA], sizeFactorsB=sf[colB],
            numPart=numPart, SHonly=FALSE, propForSigma= propForSigma, 
            shrinkTarget=shrinkTarget, shrinkQuantile=shrinkQuantile, 
            cLA=cL[colA], cLB=cL[colB], contrast=contrast,  
            keepLevelsConsistant=keepLevelsConsistant, useMMdisp,
            shrinkVariance=shrinkVariance, pairedDesign=pairedDesign, 
            pairedDesign.dispMethod=pairedDesign.dispMethod, 
            useFisher=useFisher, Dispersions=Dispersions, 
            eSlope=eSlope, plotASD=plotASD,  lwd_ASD=lwd_ASD, 
            cex_ASD=cex_ASD)
        print(Sys.time()-t1)
        dispMM=pval0$dispMM
        dispSH=pval0$dispSH
        Mean=pval0$mu
        pval=pval0$pval
    }
    cl.nm=sort(unique(c(cL[colA], cL[colB])))
    rM.A=rowMeans(as.matrix(counts[,colA])) 
    rM.B=rowMeans(as.matrix(counts[,colB]))
    l2f=log2(rM.A/rM.B)
    
    rs=data.frame(  
        Mean=Mean,
        rawMeanA=rM.A,
        rawMeanB=rM.B, 
        rawLog2FoldChange=l2f, 
        dispMM=dispMM,
        dispSH=dispSH,
        pval=pval, 
        stringsAsFactors=FALSE)
    rownames(rs)=rownames(countsTable)
    if(addRawData){
        rs=data.frame(rs, countsTable)
    }
    return(rs)
}
