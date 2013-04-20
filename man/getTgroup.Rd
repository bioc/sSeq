\name{getTgroup}
\alias{getTgroup}
\title{
This is an internal function used to calculate the shrinkage estimation when multiple shrinkage targets are considered.
}
\description{
Internal function where there are multiple shrinkage targets.
}
\usage{
getTgroup
}
\seealso{
  \code{\link{nbTestSH}}.
}
\examples{
data("countsTable")
conds <- c("A",  "B")
resSH <- nbTestSH( countsTable, conds, "A", "B", numPart=10)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
