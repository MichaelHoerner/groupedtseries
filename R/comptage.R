#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

comptage <-function( taille,sn,regime ) {

d_t <- matrix(1,taille,1)
n_ii <- matrix(0,regime,regime)
for (  i in 1 : taille ) {
    if(i==1) {
        etat_prev <- sn[i]
       } #

    etat <- sn[i]
    n_ii[etat_prev,etat] <- n_ii[etat_prev,etat] +1

    if(etat_prev!=etat) {
        d_t[i] <- 1
    } else if(i>1) {
        d_t[i] <- d_t[i-1]+1
       } #
    etat_prev <- etat

   } #

returnlist <- list(n_ii, d_t)
return(returnlist)
}
