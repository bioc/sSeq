\name{nbTestSH}
\alias{nbTestSH}
\title{
Differential Analysis based on RNA-seq experiments using Negative Binomial
 (NB) Model with Shrinkage Approach of Dispersion Estimation.
}
\description{
This is the main function calculating the exact per-gene probabilities for 
p-values. It tests the null hypothesis that the expected expression of a 
gene under two conditions are the same.
}
\usage{
nbTestSH(countsTable, conds, condA = "A", condB = "B", 
    numPart = 1, SHonly = FALSE, propForSigma = c(0, 1), 
    shrinkTarget = NULL, shrinkQuantile = NULL, plotASD = FALSE, 
    coLevels = NULL, contrast = NULL, keepLevelsConsistant = FALSE, 
    useMMdisp = FALSE, addRawData = FALSE, shrinkVariance = FALSE, 
    pairedDesign = FALSE, pairedDesign.dispMethod = "per-pair", 
    useFisher = FALSE, Dispersions = NULL, eSlope = 0.05, lwd_ASD = 4.5, 
    cex_ASD = 1.2)
}
\arguments{
  \item{countsTable}{
A data.frame or a matrix of counts in which a row represents for a gene and 
a column represents for a sample. There must be at least two columns in 
countsTable.
}
  \item{conds}{
A vector of characters representing the two conditions (or two groups). It 
must be matchable to the columns in countsTable. For example, c("A", "A", 
"B", "B") matches to a countsTable that has four columns (or samples) in 
which the first two columns are samples under condition A and the last two 
columns are samples under condition B.
}
  \item{condA}{
A character specifying the first condition in countsTable, e.g. condA="A".
}
  \item{condB}{
A character specifying the second condition in countsTable, e.g. condB="B".
}
  \item{numPart}{
An integer indicating the number of targets for the shrinkage dispersion 
estimates. ``numPart=1" is the default value. It assumes that all the genes 
share one common target, and then the method of moment estimates are shrunk 
toward one single target. When it is assumed that the genes share multiple 
targets, the value for ``numPart" is the number of targets and the grouped 
shrinkage estimates for dispersion are calculated. 
}
  \item{SHonly}{
If `SHonly' is TRUE, then the function outputs the shrinkage estimates for 
dispersion without testing the differentiation between conditions. If 
FALSE, then the function outputs a data frame including the per-gene p-
values of tests.
}
  \item{propForSigma}{
A range vector between 0 and 1 that is used to select a subset of data. It 
helps users to make a flexible choice on the subset of data when they 
believe only part of data should be used to estimate the variation among 
per-gene dispersion. An input ``propForSigma=c(0.1, 0.9)" means that the 
genes having method of moment estimates for dispersion greater than the 
10th quantile and less than the 90th quantile are used to estimate the 
dispersion variation. The default input ``propForSigma=c(0, 1)" is 
recommended. It means that we want to use all the data to estimate the 
dispersion variation.
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
  \item{plotASD}{
A logic value. If plotASD=TRUE, then the plot of average squared difference 
(ASD) versus target points is produced. The shrinkage (SH) estimates are 
obtained by shrinking the method of moment (MM) estimates toward a point 
target. In the figure, the vertical axis are ASD values when the shrinkage 
target (represented by the horizontal axis) varies within the range of 
dispersion estimates. The selected target is a small value minimizing ASD. 
}
  \item{coLevels}{
A data.frame specifying the additional factors for testing in complex 
experiments. The number of row in ``coLevels" matches the number of columns 
in countsTable. It describes the extra features or factors other than the 
two basic conditions. For example, ``conds=c("A","A","B","B")" and 
``coLevels=data.frame(sample=c(1,2,1,2))" indicate a paired design 
experiment. Column 1 and 3 in countsTable are a paired observations for 
sample 1 in two different conditions. 
}
  \item{contrast}{
A contrast vector for testing in complex experiments. The length of this 
vector equals to the number of columns in countsTable. 
}
  \item{keepLevelsConsistant}{
A logic TRUE/FALSE value. When ``coLevels" is used to indicate a paired 
design experiment, ``keepLevelsConsistant=TRUE" silences the genes that 
have different changing directions (i.e. positive and negative test 
statistics) among individual samples by setting their p-values as 1. 
}
  \item{useMMdisp}{
A logic value. When ``useMMdisp=TRUE" the method of moment (MM) estimates 
for dispersion without any shrinkage approach are used for testing the 
differentiation of genes between two conditions. 
}
  \item{addRawData}{
A logic value. When ``addRawData=TRUE", this function also outputs the 
original values of countsTable.  
}
  \item{shrinkVariance}{
A logic value. When ``shrinkVariance=TRUE", the testing is based on the 
shrinkage estimates for variance instead of dispersion.
}
  \item{pairedDesign}{
A logic value. When pairedDesign=TRUE is specified, the tests are performed 
specifically for the paired design experiment. The Null hypotheses 
\eqn{\sum_l(\mu_{gA,l}-\mu_{gB,l})=0} will be tested.  
}
  \item{pairedDesign.dispMethod}{
A character specifying the method of selecting data used for the paired 
design experiment. When the input is ``per-pair" (the default input), the 
dispersion estimates are shrunk within each pair of samples. The shrinkage 
target is different in different pair of samples. When the input is 
"pooled", firstly method of moment estimates for dispersion are obtained 
within each pair of samples, and then the average estimates across all 
pairs of samples are shrunk toward a common targets among genes. 
}
  \item{useFisher}{
A logic value specifying whether Fisher's method of combining multiple p-
values for a gene is used in the paired design experiment. In detail the 
formula of calculating the Fisher's combined p-values is 
pval\eqn{_g=\chi^2_{df=2k}(X>x)} where \eqn{k} is the number of pairs and 
\eqn{x=-2*\sum^k_{l=1}log_{e}(p_{l})}. The default input is FALSE and the 
formulae pval\eqn{_{g}=exp(\sum^k_{l=1}log_{e}(p_{l}))} is used.
}
  \item{Dispersions}{
If it is not null, then the input is a vector of known dispersion values. 
The length of the vector equals to the number of genes in the counts table. 
The default value is ``NULL".
}
  \item{eSlope}{
A positive value near to zero. When selecting the shrinkage target that is 
small and minimizing the average squared difference (ASD), the value of 
``elope" is a threshold to stop the selection steps if the absolute value 
of a local slope for the ASD is less than the threshold. The default value 
is 0.05.
}
  \item{lwd_ASD}{
A value specifying the width of the curve shown in the plot for the average 
squared difference when ``plotASD=TRUE". The default value is 4.5.
}
  \item{cex_ASD}{
A value specifying the size of label text shown in the plot for the average 
squared difference when ``plotASD=TRUE". The default value is 1.2.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{Mean }{The row per-gene averages over the values in countsTable.}
  \item{log2FoldChange }{The per-gene fold Changes between condition A and 
        B in the log2 scale.}
  \item{dispMM }{The per-gene method of moment (MM) estimates on 
        dispersion.}
  \item{dispSH }{The per-gene shrinkage (SH) estimates on dispersion.}
  \item{pval}{The per-gene p-values based on the exact tests. Smaller 
p-value indicates a higher chance of rejecting the null hypothesis that the 
expected gene expression distributes identically between the two 
conditions.}
}
\references{
Yu, D., Huber, W. and Vitek O. (2013). Shrinkage estimation of dispersion in Negative Binomial models for RNA-seq experiments with small sample size. Bioinformatics.
}
\examples{
#load a simulated data that includes a count table
data("countsTable")

#Differential analysis in sSeq.
conds <- c("A",  "B")
resSH <- nbTestSH( countsTable, conds, "A", "B")

#If users only want to calculate the SH dispersion estimates and 
#draw a mean-dispersion plot, the following scripts can be used.
library('RColorBrewer')
dispSH <- nbTestSH( countsTable, conds, "A", "B", SHonly=TRUE)
plotDispersion(dispSH)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
