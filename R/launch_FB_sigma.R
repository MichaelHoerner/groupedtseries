#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

launch_FB_sigma <-function( y,X,taille,regime_sig,theta,sigma,sn,sn_sig,p_sig,AR_lags,MA_lags,temper,dependance,d_t,dist_in,lambda ) {

eps_t <- double(length=taille)
dens <- 0
eps_t <- .C("dens_CP_ARMA_c",
              as.vector(y, mode = "double"),
              as.vector(X),
              as.integer(AR_lags),
              as.integer(MA_lags),
              as.vector(theta, mode = "double"),
              as.vector(sigma, mode = "double"),
              as.vector(sn, mode = "double"),
              as.vector(sn_sig, mode = "double"),
              as.integer(1),
              as.integer(taille),
              dens = as.double(dens),
              eps_out = as.vector(eps_t, mode="double"),
              PACKAGE="groupedtseries")$eps_out

p_lambda <- matrix(0,regime_sig, regime_sig)

if(dependance==1) {
    for (  j in 1 : regime_sig ) {
        prob <- (1-lambda[j])*matrix(1,1,regime_sig)/(regime_sig-1)
        prob[j] <- lambda[j]
        p_lambda[j,] <- prob
       } #
   } #

forward_out <- matrix(0, nrow=taille, ncol=regime_sig)
sn_move <- double(length=taille)
log_like <- 0

templist <- .C("FB_sigma",
               as.vector(eps_t, mode = "double"),
               as.integer(regime_sig),
               as.vector(sigma, mode="double"),
               as.vector(p_sig, mode = "double"),
               as.vector(runif(taille,0,1), mode="double"),
               as.double(temper),
               as.integer(dependance),
               as.vector(d_t, mode="double"),
               as.double(dist_in),
               as.double(p_lambda, mode="double"),
               as.integer(taille),
               sn_move = as.vector(sn_move, mode="double"),
               log_like_out = as.double(log_like),
               forward_out = as.vector(forward_out, mode="double"),
               PACKAGE="groupedtseries")

sn_sig <- templist$sn_move
log_like <- templist$log_like_out
forward <- templist$forward_out
forward <- matrix(forward, nrow=taille, ncol=regime_sig)

sn_sig <- round(sn_sig)+1

returnlist <- list(sn_sig, log_like, forward)
}
