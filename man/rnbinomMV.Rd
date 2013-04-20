\name{rnbinomMV}
\alias{rnbinomMV}
\title{
Randomly Generate Negative Binomial Variable with parameters mean and 
variance.
}
\description{
This function is based on the re-parameterized Negative Binomial 
distribution to generate random observations.
}
\usage{
rnbinomMV(n, mu, v)
}
\arguments{
  \item{n}{
The number of values that will be randomly generated.
}
  \item{mu}{
The expectation of the Negative Binomial distribution.
}
  \item{v}{
The variance of the Negative Binomial distribution.
}
}
\examples{
x <- rnbinomMV(50, 10, 15)
hist(x)
}
\keyword{ rnbinomMV }