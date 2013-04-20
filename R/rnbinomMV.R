rnbinomMV <- function( 
    n, 
    mu, 
    v 
){
    rnbinom( n, prob=mu/v, size = mu^2/(v-mu))
}
