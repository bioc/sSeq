\name{Sultan}
\alias{Sultan}
\docType{data}
\title{
An example of real experiment. 
}
\description{
A subset of the real experiment Sultan et al. It is used as an example for 
running some functions in this package.
}
\usage{data(Sultan)}
\format{A data.frame containing 4 columns and 52580 rows.}
\details{
It compares two biological replicates of human cell lines Ramos B and 
HEK293T with the Illumina platform.
}
\source{
\url{http://bowtie-bio.sourceforge.net/recount/}
}
\references{
Sultan, M. et al. (2008). A global view of gene activity and alternative 
splicing by deep sequencing of the human transcriptome. Science, 321, 956.

Frazee, A. et al. (2011). ReCount: A multi-experiment resource of 
analysis-ready RNA-seq gene count datasets. BMC Bioinformatics, 12, 449.
}
\examples{
data(Sultan);
head(countsTable);
}
\keyword{Sultan}
