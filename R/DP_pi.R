#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

 DP_pi <-function( param,regime,n_ii ) {

P <- matrix(0,regime,regime)
if(regime>0){
for (  i in 1 : regime ) {
    vect_a <- matrix(0,regime,1)
    for (  j in 1 : regime ) {
        if(i==j) {
            vect_a[j] <- param$alpha*param$gamma[j] + n_ii[i,j] + param$kappa
        } else {
            vect_a[j] <- param$alpha*param$gamma[j] + n_ii[i,j]
           } #
       } #

    P[i,] <- dirichlet_sample(t(vect_a))

#    param$pi[i,] <- dirichlet_sample[vect_a']
   } #
}

return(P)
}

# P <- matrix(0,regime,regime)
# for (  i in 1 : regime ) {
#     vect_a <- 0.5*matrix(1,regime,1)
#     #vect_a[i] <- 1000
#     for (  j in 1 : regime ) {
#         if(i==j) {
#             vect_a[j] <- vect_a[j]+ n_ii[i,j]
#         } else {
#             vect_a[j] <- vect_a[j] + n_ii[i,j]
#            } #
#        } #
#
#     P[i,] <- dirichlet_sample[vect_a']
#
# #    param.pi[i,] <- dirichlet_sample[vect_a']
#    } #




