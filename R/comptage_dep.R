#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

 comptage_dep <-function( taille,sn,d_t,regime,dist ) {

n_ii <- matrix(0,regime,regime)
n_ii_lambda <- matrix(0,regime,2)
for (  i in 1 : taille ) {
    if(i==1) {
        etat_prev <- sn[i]
       } #

    etat <- sn[i]
    if(d_t[i]>dist) {
        n_ii[etat_prev,etat] <- n_ii[etat_prev,etat] +1
    } else if(etat==etat_prev) {
        n_ii_lambda[etat_prev,1] <- n_ii_lambda[etat_prev,1]+1
    } else {
        n_ii_lambda[etat_prev,2] <- n_ii_lambda[etat_prev,2]+1
       } #
    etat_prev <- etat

   } #

returnlist <- list(n_ii, n_ii_lambda)
}
