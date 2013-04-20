\name{getAdjustDisp}
\alias{getAdjustDisp}
\title{
Calculate Shrinkage (SH) Estimates for Dispersion
}
\description{
In this shrinkage approach, the per-gene dispersion is considered as a 
variable in large dimensions. For example, if sequences of 30,000 genes are 
read in a RNA-seq experiment, then the dispersion variable is distributed 
in 30,000 dimensions. Firstly method-of-moment (MM) estimates on dispersion 
are calculated under the Negative Binomial (NB) modeling respectively for 
each gene. Those initial estimates are independently obtained in each 
dimension. Since RNA-seq experiments typically includes small number of 
samples (such as 1,2,3,4), the per-gene MM estimates are not reliable due 
to the limitation of sample size. We believe that there is a common 
variation shared across genes. The shrinkage approach regularizes per-gene 
dispersion estimates toward the common variation and produces robust 
estimates. Therefore in the second step, the MM estimates are shrunk 
towards an estimated target that minimizes the average squared difference 
(ASD) between the initial estimates and the shrinkage estimates. 
}
\usage{
getAdjustDisp(obs, propForSigma=c(0.5, 1), shrinkTarget=NULL, 
              shrinkQuantile=NULL, verbose=TRUE)
}
\arguments{
  \item{obs}{
A vector of initial estimates that are used to obtain the shrinkage (SH) estimates. The length of this vector must equal to the number of rows in the counts table. For example, the method-of-moment (MM) estimates for dispersion based on the Negative Binomial (NB) distribution are the initial estimates.
}
  \item{propForSigma}{
A range of percentiles that is used to identify a subset of data. It helps users to make a flexible choice on the subset of data when calculating variance of initial estimates among per-gene dispersion. A default input propForSigma=c(0, 1) is recommended. It means that we want to use all the data to estimate the variance.
}
   \item{shrinkTarget}{
A value that represents the targeted point of stabilization for shrinkage estimates on dispersion. When shrinkTarget=NULL, the point of stabilization will be calculated according to the input of shrinkQuantile. If a numeric value is input for shrinkTarget, the shrinkQuantile argument will be ignored.
}
   \item{shrinkQuantile}{
A value between 0 and 1 that represents the target quantile point of stabilization for shrinkage estimates on dispersion. When a numeric value is not provided for shrinkTarget, the shrinkQuantile argument is used. The default value is NULL and means that the function will automatically estimate the point of stabilization based on the pattern of the average squared difference (ASD) between the initial method of moment (MM) estimates and the shrinkage (SH) estimates on dispersion.
}
  \item{verbose}{
A logic value. When verbose=TRUE, the detail information will be printed in the console.
}
}
\value{
 \item{adj }{The SH estimates that shrink the input vector of obs toward the common information.}
 \item{cpm }{A data.frame that includes several summary statistics, such as the average and the variance of values in obs based on the subset controlled by the propForSigma argument.}
}
\examples{
data("countsTable");

#calculate the row means;
rM <- rowMeans(countsTable);

#calculate the row variances;
rV <- rowVars(countsTable);

#obtain an initial estimates;
disp <- (rV - rM)/rM^2;

#calculate the shrinkage estimates that shrink the initial estimates toward the common information;
dispSH <- getAdjustDisp(disp);
head(dispSH);
}
\keyword{ getAdjustDisp }
