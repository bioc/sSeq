\name{exactNBtest1}
\alias{exactNBtest1}
\title{
Perform only one exact test under the Negative Binomial modeling.
}
\description{One exact test for only one gene.}
\usage{
exactNBtest1(kA, kB, mu, disp, sA=1, sB=1, rA=0.5, rB=0.5)
}
\arguments{
  \item{kA}{
An integer matrix under condition A.
}
  \item{kB}{
An integer under matrix condition B.
}
  \item{mu}{
The expectation.
}
  \item{disp}{
The dispersion.
}
  \item{sA}{
The size factors under condition A.
}
  \item{sB}{
The size factors under condition B.
}
  \item{rA}{
Proportion of samples that are under condition A.
}
  \item{rB}{
Proportion of samples that are under condition B.
}
}
\value{
 \item{pval }{P-value.}
}
\examples{
exactNBtest1(100, 150, 125, 1.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
