#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


table_IHMM <-function( param,n_ii,regime ) {
######################
### Calcule le nombre de table
### Dep   } # de pi, alpha, kappa, n_ii
######################


m <- matrix(0,regime,regime)
if (regime>0) {
  for (  i in 1 : regime ) {
      for (  j in 1 : regime ) {
        if(n_ii[i,j] > 0){
          for (  k in 1 : n_ii[i,j] ) {
              if(i==j) {
                  prob <- (param$alpha*param$gamma[j]+param$kappa)/(k-1+param$alpha*param$gamma[j]+param$kappa)

              } else {
                  prob <- (param$alpha*param$gamma[j])/(k-1+param$alpha*param$gamma[j])
                 } #
              if(runif(1, 0, 1) < prob) {
                  m[i,j]<-m[i,j]+1
                 } #

             } #
          }
         } #
     } #
}

return(m)
}

