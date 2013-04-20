\name{getNormFactor}
\alias{getNormFactor}
\title{
Estimate size factors
}
\description{
Calculate the size factor.
}
\usage{
getNormFactor(countsTable1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countsTable1}{
A data.frame or a matrix of counts in which a row represents for a gene and 
a column represents for a sample. There must be at least two columns in 
countsTable.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Anders, S. and Huber, W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11, R106.
}
\examples{
#load a simulated data that includes a count table
data("countsTable");
getNormFactor(countsTable);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
