#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

 log_gampdf <-function( x,param_a,param_b ) {

if(param_a>0 && param_b>0 && x>0) {
    ln_f <- -param_a*log(param_b) - lgamma(param_a) + (param_a-1)*log(x) -x/param_b
} else {
    ln_f <- -inf
   } #
return(ln_f)
 }
