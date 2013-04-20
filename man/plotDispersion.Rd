\name{plotDispersion}
\alias{plotDispersion}
\title{
Drawing Dispersion-Mean plot.
}
\description{
This function is used to draw a scatter plot of dispersion versus mean of 
count table. It helps to visually inspect the dependence between the 
dispersion estimates and the mean estimates.
}
\usage{
plotDispersion(DispSH, extraOutput=NULL, plotMethod="logDisp",
    ylim1=NULL, legPos="topleft", myCol=brewer.pal(9, "Set1"), 
    tt=NULL)
}
\arguments{
  \item{DispSH}{
A data frame includes `SH', `raw', and `mus'. They are the shrinkage 
estimates of dispersion, the method of moment estimates of dispersion, and 
the estimates of mean.  This data frame is obtained using the function 
`nbTestSH' and specifying `SHonly=TRUE'. 
}
  \item{extraOutput}{
A data.frame including dispersion estimates and expectation estimates using 
another method. When users want to compare the dispersion estimates using 
two different method, this argument can be used to include the result from 
the second method. The default value is NULL. This means that no extra 
method is compared.
}
  \item{plotMethod}{
If plotMethod=``logDisp" which is the default, then both dispersion and 
mean estimates are shown in the log scale. If plotMethod=``Disp", then only 
mean estimates are shown in the log scale. 
}
  \item{ylim1}{
A vector of two values that specifies the minimum and maximum values of the 
vertical y axis in the plot. It is used to limit the presenting range of y 
axis in the plot. If ylime1=NULL then the range of the shrinkage estimates 
of dispersion is used. 
}
  \item{legPos}{
A character indicating the position of legend in the plot. The value of 
this argument can be "topleft", "topright", "bottomleft" and "bottomright".
}
  \item{myCol}{
A vector of colors corresponding to the dispersions estimated using 
different methods.
}
  \item{tt}{
A character representing the title of the plot, which is shown on the top 
in the plot.
}
}
\examples{
data("countsTable")
conds <- c("A", "B")
dispSH <- nbTestSH( countsTable, conds, "A", "B", SHonly=TRUE)

library('RColorBrewer')
plotDispersion(dispSH, legPos="topleft")
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
