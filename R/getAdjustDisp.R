getAdjustDisp <- function( 
    obs, 
    propForSigma=c(0.5, 1), 
    shrinkTarget=NULL, 
    shrinkQuantile=NULL, 
    verbose=TRUE)
{
    obs[is.na(obs)]=0
    if(is.null(shrinkTarget)){
        upBound=quantile(obs, prob= shrinkQuantile, na.rm=T)
        if(verbose){
            print(paste("shrink toward ", shrinkTarget, " (", shrinkQuantile,
                "th quantile).", sep=''))
        } 
    } else {
        upBound=shrinkTarget
        if(verbose){
            print(paste("shrink toward ", shrinkTarget,".", sep=''))
        }
    }

    if(is.null(propForSigma)){
        subobs=obs[obs>=upBound & obs<=quantile(obs,prob=0.999)]
        S.mat=var(subobs, na.rm=T)
    } else if (length(propForSigma)==2){
        subobs=obs
        rg=quantile(subobs[is.finite(subobs)], na.rm=T, prob= propForSigma)
        subobs=subobs[subobs>=rg[1] & subobs<=rg[2]]
        S.mat=var(subobs[is.finite(subobs)], na.rm=T)
    } else if (length(propForSigma)==1 & is.numeric(propForSigma)){
        S.mat=propForSigma
    } else if ( is.na(propForSigma) ){
        subobs=obs[is.finite(obs)]
        S.mat=var(subobs[is.finite(subobs)], na.rm=T)
    } else { 
        stop(paste("if don't know the empirical value on the variance of", 
            "dispersion, please set it as NULL."))
    }

    cmp=data.frame(mean=mean(obs, na.rm=T), sigmasq.part=S.mat)
    mean.mat=rep(upBound, length(obs))
    dif.mat=obs - mean.mat
    dif2.mat=sum(dif.mat^2)
    deta=1- ((length(obs)-2)*S.mat/(dif2.mat))   
    jsDiff=pmax(0, deta) * dif.mat 
    jsest=jsDiff + mean.mat
    return(list( adj=jsest, cmp=cmp))
}
