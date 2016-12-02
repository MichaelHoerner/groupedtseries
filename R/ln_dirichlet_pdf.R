#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


 ln_dirichlet_pdf <-function( X,alpha ) {

log_dens <- 0
dimension <- dim( as.data.frame(alpha))
dimension <- max(dimension)

for (  i in 1 : dimension ) {
    log_dens <- log_dens + (alpha[i]-1)*log(X[i]) - lgamma( alpha[i])
   } #
log_dens <- log_dens + lgamma( sum(alpha))

return(log_dens)
}
