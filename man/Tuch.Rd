\name{Tuch}
\alias{Tuch}
\docType{data}
\title{
An example of real experiment. 
}
\description{
A subset of the real experiment Tuch et al. It is used as an example for 
running some functions in this package.
}
\usage{data(Tuch)}
\format{A data.frame containing 6 columns and 10453 rows.}
\details{
It compares the expression of genes in normal human tissues and in tissues 
with oral squamous cell carcinoma. The experiment had a paired design in 
that pairs of normal and tumor samples were obtained from three patient. 
The six libraries were sequenced using the SOLiD platform.
}
\source{
The table of read counts was downloaded from GEO (accession GSE20116).
}
\references{
Tuch, B. et al. (2010).	Tumor transcriptome sequencing reveals allelic 
expression imbalances associated with copy number alterations.	PloS One, 
5, e9317.
}
\examples{
data(Tuch);
head(countsTable);
}
\keyword{Tuch}
