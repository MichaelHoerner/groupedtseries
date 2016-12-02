#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

Draw_mu_theta <-function( theta,Sigma_theta_inv,taille_ARMA,regime,mean_mu_prior_theta,Sigma_mu_prior_theta_inv ) {

Sig_mu <- solve(regime*Sigma_theta_inv + Sigma_mu_prior_theta_inv)

mean_mu <- Sigma_mu_prior_theta_inv%*%mean_mu_prior_theta
if (regime > 0) {
  for (  q in 1 : regime ) {
      mean_mu <- mean_mu + Sigma_theta_inv%*%(theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)])
     } #
}
mean_mu <- Sig_mu%*%mean_mu

mu_theta <- as.matrix(mvrnorm(1, mean_mu,Sig_mu))

return(mu_theta)

}
