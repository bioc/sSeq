\name{sim}
\alias{sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Simulated Data
}
\description{
This function is used to approximate the real experiment and to generate 
simulated counts table based on Negative Binomial distribution.
}
\usage{
sim(ngenes, true_mean1, conds, 
  alpha = function(m) {rep(0.1, length(m))}, 
  mean_DE = 0, sd_DE = 2, s0 = NULL, s0_mean = 2, s0_sd = 3, 
  true_isDE_proportion = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ngenes}{
The total number of genes or rows in the simulated counts table.
}
  \item{true_mean1}{
The expected gene expression (\eqn{\mu_{g}}) in a library or a sample. The 
length of this vector equals `ngenes'. It can generated from either random 
distributions or averages of counts table from a real experiment.
}
  \item{conds}{
A vector of characters representing the two conditions (or two groups). It 
must be matchable to the columns in countsTable, e.g., c("A", "A", "B", 
"B") matches to a countsTable that has four columns (or samples) in which 
the first two columns are samples under condition A and the last two 
columns are samples under condition B. 
}
  \item{alpha}{
A function used to generate the true dispersion values. The default 
function generates a constant 0.1 for all the genes. It can also be a 
function specifying the dependence between dispersion and mean.
}
  \item{mean_DE}{
A true mean value of \eqn{\epsilon} in \eqn{\mu_{gB}=\mu_{g}/exp(\epsilon)} 
where \eqn{\epsilon} follows a Normal distribution.
}
  \item{sd_DE}{
A true standard deviation of \eqn{\epsilon} in 
\eqn{\mu_{gB}=\epsilon\mu_{g}} where \eqn{\epsilon} follows a Normal 
distribution.
}
  \item{s0}{
The true size factors for samples. The length of this vector equals to the 
length of the vector `conds'. 
}
  \item{s0_mean}{
If the true size factors for samples are not defined for `s0', then the 
true size factors are assumed to follow a Normal distribution with mean as 
the value for `s0_mean'.
}
  \item{s0_sd}{
If the true size factors for samples are not defined for `s0', then the 
true size factors are assumed to follow a Normal distribution with standard 
deviation as the value for `s0_sd'.
}
  \item{true_isDE_proportion}{
The proportion of genes that are truly different. The default value is 0.3.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The function outputs a list including the simulated counts table, a vector 
with TRUE of FALSE values indicating the truly differentiating genes, the 
true mean values, the true variance values, and the true dispersion values. 
}
\references{
Yu, D., Huber, W. and Vitek O. (2013). Shrinkage estimation of dispersion in Negative Binomial models for RNA-seq experiments with small sample size. Bioinformatics.
}
\note{
We acknowledge Dr. Simon Anders since he provided the details for 
simulation in the manual of DESeq package.
}
\examples{
ng = 10000;
sim1 <- sim(ngenes=ng, conds=c("A","A","B","B"),
   true_mean1=round(rexp(ng, rate=1/200)), alpha=function(m){1/(m+100)},
   mean_DE=2, sd_DE=1, s0=runif(4, 0, 2) ); 
true_isDE <- sim1$true_isDE;    
countsTable <- sim1$countsTable;
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
