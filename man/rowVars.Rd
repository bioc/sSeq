\name{rowVars}
\alias{rowVars}
\title{
Calculating the sample variance within each row of A matrix
}
\description{
This function helps to obtain row-wise estimation across columns.
}
\usage{
rowVars(x)
}
\arguments{
  \item{x}{
A matrix or data.frame that includes multiple columns.
}
}
\value{
A vector showing the per-row variance estimates for the matrix or 
data.frame.
}
\examples{
x <- matrix(rnorm(10), 5)
rowVars(x)
}
\keyword{ rowVars }
