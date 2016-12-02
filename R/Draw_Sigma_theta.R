#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

Draw_Sigma_theta <-function( theta,mu_theta,V_prior_inv,dv_prior,regime,taille_ARMA ) {

dv_post <- dv_prior + regime
V_post <- V_prior_inv
dimension <- dim(V_post)


for (  q in 1 : regime ) {
    V_post <- V_post + matrix((theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)]-mu_theta)%*%t(theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)]-mu_theta),dimension[1],dimension[2]) #traverse was defined the other way around
   } #

V_post <- solve(V_post)

Sigma_theta_inv <- rWishart(1, dv_post, V_post)[,,1]

return(Sigma_theta_inv)
}
