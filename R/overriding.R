#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

overriding <-function( param,regime,etat_init ) {
###########################
##### Calcule le nombre de repas impos? par kappa
##### Dep   } # de m, alpha, kappa
###########################

r <- matrix(0,regime,1)
rho <- param$kappa/(param$alpha+param$kappa)
if (regime>0) {
for (  i in 1 : regime ) {
    r[i] <- rbinom(1, param$m[i,i], rho/(rho+(1-rho)*param$gamma[i]))
   } #
if(r[etat_init]>0) {
    r[etat_init] <- r[etat_init]-1
   } #
}

return(r)
}
