\name{countsTable}
\alias{countsTable}
\docType{data}
\title{
An Example Simulation Data
}
\description{
A subset of simulated data. It is used as an example for running some 
functions in this package.
}
\usage{data(countsTable)}
\format{
  The format is:
 num [1:10000, 1:2] 90 155 13347 254 228 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:10000] "1_FALSE" "2_FALSE" "3_TRUE" "4_FALSE" ...
  ..$ : chr [1:2] "A1" "B1"
}
\details{
A simulation counts table.
}
%\source{
%%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(countsTable)
head(countsTable)
}
\keyword{countsTable}
