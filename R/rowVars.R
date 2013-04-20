rowVars <- function(
    x
){
    apply(x, 1, var, na.rm=T)
}
