#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

DP_gamma <-function( etha,m_bar,regime ) {

vect_a <- etha/regime*matrix(1,regime,1)

for (  i in 1 : regime ) {
    vect_a[i] <- vect_a[i] + sum(m_bar[,i])
   } #
gamma <- dirichlet_sample(vect_a)

return(gamma)

}
