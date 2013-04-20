\name{Hammer2months}
\alias{Hammer2months}
\docType{data}
\title{
An example of real experiment. 
}
\description{
A subset of the real experiment Hammer et al. It is used as an example for 
running some functions in this package.
}
\usage{data(Hammer2months)}
\format{A data.frame containing 4 columns and 29516 rows.}
\details{
It compares gene expression in rat strains Sprague Dawley and L5 SNL 
Sprague Dawley 2 at the end of two months in a factorial design. Two 
distinct biological libraries per condition were quantified using the 
Illumina platform.
}
\source{
\url{http://bowtie-bio.sourceforge.net/recount/}
}
\references{
Hammer, P. et al. (2010). mRNA-seq with agnostic splice site discovery for 
nervous system transcriptomics tested in chronic pain. Genome Res., 20, 
847-860.

Frazee, A. et al. (2011). ReCount: A multi-experiment resource of 
analysis-ready RNA-seq gene count datasets. BMC Bioinformatics, 12, 449.
}
\examples{
data(Hammer2months);
head(countsTable);
}
\keyword{Hammer2months}
