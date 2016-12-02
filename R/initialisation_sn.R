#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


 initialisation_sn <-function( nb_init,nb_obs,regime,CP ) {
# On tire al?atoirement un certain nombre de vecteur sn



# initialisation des variables
init_sn <- matrix(1,nb_obs,nb_init)
init_tau <- matrix(0,regime,nb_init)
n_ii <- array(0,dim = c(regime,regime,nb_init))
interval <- nb_obs/(regime)

if(regime>1) {

    if(CP==1) {
    # Tirage al?atoire de cassure pour un CP
        tau <- matrix(0,regime,1)
        for ( i in 1 : nb_init ) {
            for ( q in 1 : regime ) {
                if(q==1) {
                    tau[q] <- round(interval) + round((runif(1, 0, 1)-0.5)*interval)
                } else {
                    tau[q] <- tau[q-1] + round(interval) + round((runif(1, 0, 1) -0.5)*interval)
                    if(tau[q]>=nb_obs && q<regime) {
                        tau[q] <- nb_obs-q
                    } else if(q==regime) {
                        tau[q] <- nb_obs
                    } #
                } #
             } #
            etat <- 1
            for (  q in 1 : nb_obs ) {
                if(q>tau[etat]) {
                    etat <- etat +1
                   } #
                init_sn[q,i] <- etat
               } #
            init_tau[,i] <- tau
           } #
    } else {
        interval <- interval/4
        for (  i in 1 : nb_init ) {
            etat <- 1

            unif <- ceiling((runif(1, 0, 1) +0.3)*interval)
            duree <- 0
            for (  q in 1 : nb_obs ) {
                if(duree>unif) {
                    duree <- 0
                    unif <- ceiling((runif(1, 0, 1) +0.3)*interval)
                    etat <- ceiling(runif(1, 0, 1) *regime)
                   } #
                init_sn[q,i] <- etat
                if(q>1) {
                    n_ii[init_sn[q-1,i],etat,i] <- n_ii[init_sn[q-1,i],etat,i] +1
                } else {
                    n_ii[1,etat,i] <- 1
                   } #
                duree <- duree +1
               } #
           } #
       } #
} else {
    for (  i in 1 : nb_init ) {
        init_tau[1,i] <- nb_obs
        n_ii[1,i] <- nb_obs
       } #
   } #

returnvalues <- list(init_sn, init_tau, n_ii)

return(returnvalues)

}
