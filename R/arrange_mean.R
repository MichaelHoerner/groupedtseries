#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

arrange_mean <-function( regime_max,n_ii,taille,sn,theta,taille_ARMA ) {

exist <- matrix(0,regime_max,1)

if((sum(n_ii[1,]) + sum(n_ii[,1]))<=1) {
    exist[1] <- 0
} else {
    exist[1] <- 1
   } #

for (  i in 2 : regime_max ) {
    if(sum(n_ii[i,])+sum(n_ii[,i])==0) {
        exist[i] <- 0
    } else {
        exist[i] <- 1
       } #
   } #

it <- 0

for (  i in 1 : regime_max ) {
    if(exist[i]==1) {
        it <- it+1
        ## Changer sn
         if(i!=it) {
            for (  j in 1 : taille ) {
                if(sn[j]==i) {
                    sn[j]<-it
                   } #
               } #
            } #


       } #
   } #

regime <- sum(exist)
beta_new  <- matrix(0,regime_max*taille_ARMA,1)
ind <- 1
ind_off <- 1
for (  i in 1 : regime_max ) {
    if(exist[i]==1) {
        beta_new[((ind-1)*taille_ARMA+1):(ind*taille_ARMA)] <- theta[((i-1)*taille_ARMA+1):(i*taille_ARMA)]
        ind <- ind+1
    } else {
        beta_new[((regime+ind_off-1)*taille_ARMA+1):((regime+ind_off)*taille_ARMA)] <- theta[((i-1)*taille_ARMA+1):(i*taille_ARMA)]
        ind_off <- ind_off+1
       } #

   } #

returnlist <- list(beta_new, regime, sn)
return(returnlist)
}


