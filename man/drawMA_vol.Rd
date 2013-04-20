\name{drawMA_vol}
\alias{drawMA_vol}
\title{
Draw MA Plot and Volcano Plot
}
\description{
Based on the count table and the p-values, this function can be used to 
draw a MA plot of the log2 ratios versus the log2 averages upon means of 
gene expression in condition A and B, and a volcano plot of negative log2 
p-values versus the log2 ratios. 
}
\usage{
drawMA_vol(y, groups2, pv, cutoff=NULL, xlab1="(log2(A)+log2(B))/2",
           ylab1="log2(A)-log2(B)", tt1="MA plot", tt2="volcano plot", 
           log2FoldChange =NULL, col1=c("black","red"))
}
\arguments{
  \item{y}{
A count table in which row represents genes and column represents samples.
}
  \item{groups2}{
A vector indicates the two groups information of samples. It must match to 
the column in the count table, which is the input for y. For example, 
groups2=c("A","A","B","B") when the first two columns in the count table 
are the two samples from condition A, and the second two columns in the 
count table are the two samples from condition B.
}
  \item{pv}{
A vector of per-gene p-values based on the count table. The order of genes 
in pv does matter. It must be the same as the order of genes in the count 
table.
}
  \item{cutoff}{
A value used as a threshold for per-gene p-values to decide the genes that 
are differentially expressed between two conditions. If NULL, the cutoff 
value is calculated so that the red dots in the MA plot and volcano plot 
represent the first 5% of genes ranked by p-values.
}
  \item{xlab1}{
A character indicating the label of x axis in MA plot.
}
  \item{ylab1}{
A character indicating the label of y axis in MA plot.
}
  \item{tt1}{
A character indicating the title of the MA plot.
}
  \item{tt2}{
A character indicating the title of the volcano plot.
}
  \item{log2FoldChange}{
A vector of fold changes in log2 scale. It will be calculated automatically 
when ``log2FoldChange=NULL". 
}
  \item{col1}{
A vector with two values including the colors of points. The first color in 
``col1" is the color for the points that are non-differentially changed. 
The second value in ``col1" is the color for the points that are 
differentially changed. The default is c(``black", ``red").
}
}
\examples{
x <- matrix(rnorm(4000, 10), ncol=4)
px <- apply(x, 1, function(y){t.test(y[1:2], y[3:4])$p.value})
drawMA_vol(x, c("A","A","B","B"), px, cutoff=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
