\name{getT}
\alias{getT}
\title{
Estimate the shrinkage target based on the initial estimates
}
\description{
This function is recommended to estimate the shrinkage target.
}
\usage{
getT(countsTable, sizeFactors = NULL, q.vec = NULL, plotASD = FALSE, 
    numPart = 1, propForSigma = c(0, 1), verbose = TRUE, 
    shrinkTarget = NULL, shrinkQuantile = NULL, shrinkVar = FALSE, 
    eSlope = 0.05, disp = NULL, dispXX = NULL, normalize = FALSE, 
    lwd1 = 4.5, cexlab1 = 1.2)
}
\arguments{
  \item{countsTable}{
A data.frame or a matrix of counts in which a row represents for a gene and 
a column represents for a sample. There must be at least two columns in 
countsTable.
}
  \item{sizeFactors}{
A vector of values around 1 which are used to normalize between samples or 
libraries. The length of this vector equals to the number of columns in 
countsTable.
}
 \item{q.vec}{
A vector of sequence defines the quantiles. When q.vec=NULL, this function 
will generate a sequence for q.vec using seq(0.05, 0.995, 0.005).
}
  \item{plotASD}{
A logic value. If plotASD=TRUE, then the plot of ASD versus target points 
will be drawn. The SH estimates are obtained by shrinking the MM estimates 
toward a target point. Different SH estimates are generated using different 
target points. The target point that helps produce a small and stable 
averaged squared difference (ASD) between the MM estimates and the SH 
estimates is the point that approximates the common information across per-
gene dispersion. This target point is termed as the point of stabilization.
}
  \item{numPart}{
An integer indicates the number of groups for dispersion estimation. 
`numPart=1' is the default value. It assumes that most of the genes share 
one point of stabilization (POS), and calculates the SH estimates without 
separating data into groups. When we assumes that genes can share different 
points of stabilization, the grouped SH estimates on dispersion can be be 
utilized. In this situation, users need to provide a number indicating the 
number of POS.
}
  \item{propForSigma}{
A range vector between 0 and 1 that is used to select a subset of data. It 
helps users to make a flexible choice on the subset of data when they 
believe only part of data should be used to estimate the variation among 
per-gene dispersion. A default input propForSigma=c(0, 1) is recommended. 
It means that we want to use all the data to estimate the variation.
}
  \item{verbose}{
A logic value. When verbose=TRUE, the detail information will be printed in 
the console.
}
  \item{shrinkTarget}{
A value for the shrinkage target of dispersion estimates. If 
``shrinkTarget=NULL" and ``shrinkQuantile" is a value instead of NULL, then 
the quantile value for ``shrinkQuantile" is converted into the scale of 
dispersion estimates and used as the target. If both of them are NULL, then 
a value that is small and minimizes the average squared difference is 
automatically used as the target value. If both of them are not NULL, then 
the value of ``shrinkTarget" is used as the target.
}
  \item{shrinkQuantile}{
A quantile value for the shrinkage target of dispersion estimates. If 
``shrinkTarget=NULL" and ``shrinkQuantile" is a value instead of NULL, then 
the quantile value for ``shrinkQuantile" is converted into the scale of 
dispersion estimates and used as the target. If both of them are NULL, then 
a value that is small and minimizes the average squared difference is 
automatically used as the target value. If both of them are not NULL, then 
the value of ``shrinkTarget" is used as the target.
}
  \item{shrinkVar}{
A logic value. When ``shrinkVariance=TRUE", the testing is based on the 
shrinkage estimates for variance instead of dispersion.
}
  \item{eSlope}{
A positive value near to zero. When selecting the shrinkage target that is 
small and minimizing the average squared difference (ASD), the value of 
``elope" is a threshold to stop the selection steps if the absolute value 
of a local slope for the ASD is less than the threshold. The default value 
is 0.05.
}
  \item{disp}{
A vector of initial estimates of dispersions. The length of this vector 
equals to the number of rows in countsTable.
}
  \item{dispXX}{
A vector of normalized mean expression. The length of this vector equals to 
the number of rows in countsTable. 
}
  \item{normalize}{
A logic value. When estimating the shrinkage target based on the average 
squared difference (ASD) between the shrinkage estimates and the initial 
estimates, the initial estimates and ASD are normalized when 
``normalize=TRUE".
}
  \item{lwd1}{
A value specifying the width of the curve shown in the plot for the average 
squared difference when ``plotASD=TRUE". The default value is 4.5.
}
  \item{cexlab1}{
A value specifying the size of label text shown in the plot for the average 
squared difference when ``plotASD=TRUE". The default value is 1.2.
}
}
\value{
  \item{target }{The estimated point for stabilization that represents the 
common in formation across per-gene dispersion. }
  \item{q }{A value that shows the quantile of the target value across per-
gene dispersion. }
}
%\author{
%Danni Yu
%}
\examples{
#load a simulated data that includes a count table
data("countsTable")
conds <- c("A",  "B")
getT(countsTable, plotASD=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
