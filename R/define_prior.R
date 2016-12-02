#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


 define_prior <-function( taille_ARMA,CP_or_MS ) {

  if(nargs()<2) {
    CP_or_MS <- "MS"
   } #

Struct_Prior <- list()

######
### Prior de sigma : sigma^[-2] ~ G[e_prior,f_prior] et donc E[.] <-
### e_prior*f_prior et V[.] <- e_prior*f_prior^2
######
Struct_Prior$e_prior <- 1/4
Struct_Prior$f_prior <- 4
Struct_Prior$c_0 <- 10#24
Struct_Prior$d_0 <- 1/5#1/3
Struct_Prior$MH_e <- 1
Struct_Prior$rho_e <- 2#5


######
### Prior de beta [coefficients de la moyenne] : beta ~
### N[mu_theta,Sigma_theta] et mu_theta ~
### N[mean_mu_prior_theta,Sigma_mu_prior_theta] et Sigma_theta^[-1] ~
### W[V_prior,dv_prior]
######

dv_prior <- 5#taille_ARMA+2# +20
V_prior_inv <- diag(1,taille_ARMA)#*dv_prior
V_prior_inv[1,1]  <- 5#2#5
for (  q in 2 : taille_ARMA ) {
    V_prior_inv[q,q]  <- 5#20#0.5
   } #
Struct_Prior$V_prior_inv <- V_prior_inv*dv_prior
Struct_Prior$dv_prior <- dv_prior
Struct_Prior$mean_mu_prior_theta <- matrix(0,taille_ARMA,1)
Struct_Prior$Sigma_mu_prior_theta_inv <- diag(1,taille_ARMA)*10#*100
Struct_Prior$Sigma_theta_inv <- diag(1,taille_ARMA)


Struct_Prior$prior_eta_a <- 1#1#10#5
Struct_Prior$prior_eta_b <- 10#10#0.5#1
switch(CP_or_MS,
       MS={
        Struct_Prior$prior_rho_a <- 10#10#10#10#15#10#10000#10#15 #10 Jochmann
       },
       random={
        Struct_Prior$prior_rho_a <- 100
       },
       {
        Struct_Prior$prior_rho_a <- 1000#10#10#10#15#10#10000#10#15 #10 Jochmann
       })

Struct_Prior$prior_rho_b <- 1#1#6#1
Struct_Prior$prior_alpha_kappa_a <- 1#1#1000#125
Struct_Prior$prior_alpha_kappa_b <- 10#10#5#1/5#1#5#5 #1/5

Struct_Prior$hier_rho_a <- 1
Struct_Prior$hier_rho_b <- 100

return(Struct_Prior)
}
