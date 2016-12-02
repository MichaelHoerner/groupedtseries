#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

 arrange_sig <-function( regime_max,n_ii,taille,sn,sigma ) {

exist <- matrix(0,regime_max,1)

if(sum(n_ii[1,])+ sum(n_ii[,1])<=1) {
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
sigma_new <- matrix(0,regime_max,1)
regime <- sum(exist)

ind <- 1
ind_off <- 1

for (  i in 1 : regime_max ) {
    if(exist[i]==1) {
        sigma_new[ind] <- sigma[i]
        ind <- ind+1
    } else {
        sigma_new[ind_off+regime] <- sigma[i]
        ind_off <- ind_off+1
       } #

   } #

returnlist <- list(sigma_new, regime, sn)
}


