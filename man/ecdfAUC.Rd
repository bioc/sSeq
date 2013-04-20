\name{ecdfAUC}
\alias{ecdfAUC}
\title{
Draw Empirical Cumulative Density Function (ECDF) plot
}
\description{
This function is used to draw Empirical CDF plot. It relies on the trapz 
function in the caTools package. A user needs to install the caTool library first. 
}
\usage{
ecdfAUC(dd, col.line=NULL, main="ECDF", cex.leg=1, drawRef=FALSE, 
       rm1=FALSE, lineType=NULL, addLeg=TRUE, xlab="p-value", 
       ylab="ECDF", cex.axis=1.5, cex.main=1.8, cex.lab=1.2, 
      axis.padj=c(-1, 1), lab.padj=c(-1.5, 1), lwd=1, box.lwd=1.2)
}
\arguments{
  \item{dd}{
A data frame of p-values in which a column represents the p-values or 
posterior probabilities resulted by a method.
}
   \item{col.line}{
A vector of color characters. The default is NULL and this function 
automatically assigns the color for each cover shown in the ECDF plot.
}
  \item{main}{
The title of the plot.
}
  \item{cex.leg}{
An integer specifying size of the legend in the the plot.
}
  \item{drawRef}{
If TRUE, then a gray 45 degree line will be added in the plot.
}
  \item{rm1}{
If users believe that the p-values equal to 1 belong to the different group 
of the others, and want to exclude them from the calculation of empirical 
CDF, then use rm1=TRUE. 
}
  \item{lineType}{
A vector of integers indicating the type of lines used for the methods.
}
  \item{addLeg}{
If ``TRUE" then a legend box with legend is added to the figure.
}
  \item{xlab}{
Label of x axis.
}
  \item{ylab}{
Label of y axis.
}
  \item{cex.axis}{
The size of labels on the axes.
}
  \item{cex.main}{
A characteristic string indicating the size of the main.
}
  \item{cex.lab}{
The size for the labels on x and y.
}
  \item{axis.padj}{
The perpendicular adjustment of ticks.
}
  \item{lab.padj}{
The perpendicular adjustment of labels for an axix.
}
  \item{lwd}{
The width of the line shown in a figure.
}
  \item{box.lwd}{
The width of the box line in a figure.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%}
%\references{
%}
%\author{
%  ~~who you are~~
%}
%\note{
%  ~~further notes~~
\examples{
x<-data.frame(A=runif(100), B=rbeta(100, 0.5, 1.2))
ecdfAUC(x);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
