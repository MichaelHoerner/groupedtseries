#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


considered_table <-function( m,r,regime ) {

m_bar <- m
if (regime>0){
for (  i in 1 : regime ) {
    m_bar[i,i] <- m_bar[i,i]-r[i]
   } #
}
return(m_bar)
}
