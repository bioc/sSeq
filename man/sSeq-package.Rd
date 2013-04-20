\name{sSeq-package}
\alias{sSeq-package}
\alias{sSeq}
\docType{package}
\title{
Shrinkage estimation of dispersion in Negative Binomial models for RNA-seq 
experiments with small sample size
}
\description{
This package is to discover the genes that differentially expressed between 
two conditions based on RNA-seq experiments. Gene expression is measured in 
counts of transcripts and modeled with the Negative Binomial (NB) 
distribution using a shrinkage approach for dispersion estimation. The 
method of moment (MM) estimates for dispersion are simply shrunk toward a 
target, which minimizes the average squared difference between the 
shrinkage estimates and the initial estimates. The exact per-gene 
probability under the NB model is calculated, and used to test the 
hypothesis that the expected expression of a gene in two conditions are not 
different.
}
\details{
\tabular{ll}{
Package: \tab sSeq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-02-25\cr
}
}
\author{
Danni Yu <dyu@purdue.edu>, Wolfgang Huber <whuber@embl.de> and Olga Vitek <ovitek@purdue.edu> 
}
\references{
Shrinkage estimation of dispersion in Negative Binomial models for RNA-seq 
experiments with small sample size
}
\keyword{ package }
\examples{
#load a simulated data that includes a count table
data("countsTable")

#calculate the p-values using the shrinkage approach.
conds <- c("A", "B")
resJS <- nbTestSH( countsTable, conds, "A", "B")
}
