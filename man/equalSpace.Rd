\name{equalSpace}
\alias{equalSpace}
\title{
Calculate Grouped Shrinkage Estimates
}
\description{
This is an internal function. When the local mean-dispersion dependence is 
present, data can be separated into groups based on the means. The windows 
used to partition groups have equal width upon each other. The shrinkage 
(SH) estimates for dispersion will be calculated within each group. For 
example, when range of the per-gene mean is 1 and 3000, if data will be 
separated into 3 groups, then group 1 includes the genes having mean values 
between 1 and 1000, group 2 includes the genes having mean values between 
1001 and 2000, and group 3 includes the genes having mean values between 
2001 and 3000. The SH estimates will be calculated within each of the 3 
groups, respectively.
}
\usage{
equalSpace(y, x, numcls=1, propForSigma=c(0, 1), shrinkTarget=NULL, 
          shrinkQuantile=0.975, vb=TRUE)
}
\arguments{
  \item{y}{
A vector including the initial values that will be regularized. For example, it can be the per-gene method of moment (MM) estimates for dispersion based on the Negative Binomial distribution for the counts table.
}
  \item{x}{
A vector that will be used to separate data into groups. For example, it can be the per-gene averages for the counts table. 
}
  \item{numcls}{
An integer that indicates the number of groups to be considered. The default value is 1. 
}
  \item{propForSigma}{
A range vector between 0 and 1 that is used to select a subset of data. It helps users to make a flexible choice on the subset of data when they believe only part of data should be used to estimate the variation among per-gene dispersion. A default input propForSigma=c(0, 1) is recommended. It means that we want to use all the data to estimate the variation.
}
   \item{shrinkTarget}{
A value that represents the targeted point of stabilization for shrinkage estimates on dispersion. When shrinkTarget=NULL, the point of stabilization will be calculated according to the input of shrinkQuantile. If a numeric value is input for shrinkTarget, the shrinkQuantile argument will be ignored.
}
   \item{shrinkQuantile}{
A value between 0 and 1 that represents the target quantile point of stabilization for shrinkage estimates on dispersion. When a numeric value is not provided for shrinkTarget, the shrinkQuantile argument is used. The default value is NULL and means that the function will automatically estimate the point of stabilization based on the pattern of the average squared difference (ASD) between the initial method of moment (MM) estimates and the shrinkage (SH) estimates on dispersion.
}
  \item{vb}{
A logic value. When verbose=TRUE, the detail information will be printed in the console.
}
}
\value{
This function returns a vector of shrinkage estimate on the basis of y.
}
\author{
Danni Yu
}
\examples{
data("countsTable");

#calculate the row means;
rM <- rowMeans(countsTable);

#calculate the row variances;
rV <- rowVars(countsTable);

#calculate the method-of-moment estimation on dispersions;
disp <- (rV - rM)/rM^2;

#calculate SH estimates in 3 groups;
disp3 <- equalSpace(disp, rM, 3);
head(disp3);
}
\keyword{ kw1 }
