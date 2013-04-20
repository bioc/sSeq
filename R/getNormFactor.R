getNormFactor <- function(
    countsTable1
){
    countsTable1.log=log(countsTable1) 
    row.mean1.log=rowMeans(countsTable1.log)
    geo.dev1.log=countsTable1.log-row.mean1.log
    apply(geo.dev1.log, 2, function(x){ exp( median(x[is.finite(x)]) ) })        
}
