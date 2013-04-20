\name{getQ}
\alias{getQ}
\title{
Estimate the shrinkage target based on the quantiles of initial targets
}
\description{
The shrinkage target is estimated.
}
\usage{
getQ(countsTable, sizeFactors=NULL, q.vec=NULL, plotASD=FALSE, 
    numPart=1, propForSigma=c(0, 1), verbose=TRUE, shrinkTarget=NULL, 
    shrinkQuantile=NULL)
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
gene dispersion. 
}
  \item{numPart}{
An integer indicates the number of groups for dispersion estimation. 
`numPart=1' is the default value. It assumes that most of the genes share 
one point of stabilization (POS), and calculates the SH estimates without 
separating data into groups. When we assumes that genes can share different 
targets, the grouped SH estimates on dispersion can be be utilized. In this 
situation, users need to provide a number indicating the number of POS.
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
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{target }{The estimated point for stabilization that represents the 
common in formation across per-gene dispersion. }
  \item{q }{A value that shows the quantile of the target value across 
per-gene dispersion. }
}
\examples{
#load a simulated data that includes a count table
data("countsTable")
conds <- c("A",  "B")
getQ(countsTable, plotASD=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
