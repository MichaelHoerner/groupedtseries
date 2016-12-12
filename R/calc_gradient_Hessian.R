#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


calc_gradient_Hessian <-function( y,X,theta,sigma,sn,sn_sig,regime,AR_lags,MA_lags,inv_Sig_prior,temper, ts_segments) {

  number_ts_segments <- nrow(ts_segments)

##if(MA_lags>1 || AR_lags==0) {
##    warning("sMALA does not work for a number of error lags greater than 1")
##} else {
  taille <- max( dim( as.data.frame(y)))
    dens <- 0
    eps_t <- double(length=taille)
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
                  as.integer(number_ts_segments),
                  as.vector(ts_segments$index_AR_start, mode="integer"),
                  as.vector(ts_segments$index_AR_end, mode="integer"),
                  as.vector(ts_segments$length_AR, mode="integer"),
                  dens = as.double(dens),
                  eps_out = as.vector(eps_t),
                  PACKAGE="groupedtseries")$eps_out

    ##############################################
    #### Derivative of epsilon_t ####
    ##############################################

    taille_ARMA <- 1+AR_lags+MA_lags
    derivative <- matrix(0,taille,regime*taille_ARMA)
    Jacobien <- matrix(0,taille,regime*taille_ARMA)
    aide_Hessian <- matrix(0,taille,regime*taille_ARMA)
    etat <- sn[1]
    derivative[1,etat] <- -1
    Jacobien[1,etat] <- -eps_t[1]*derivative[1,1]/sigma[sn_sig[1]]
    std_sigma <- sqrt(sigma)
    for (  t in 2 : taille ) {
        etat <- sn[t]
        etat_sig <- sn_sig[t]
        for (  q in 1 : regime ) {
            ind <- (q-1)*taille_ARMA
            ind2 <- (etat-1)*taille_ARMA
            derivative[t,ind+1] <- -1*(etat==q) - theta[ind2+3]*derivative[(t-1),(ind+1)]
            Jacobien[t,ind+1] <-  -temper*eps_t[t]*derivative[t,ind+1]/sigma[etat_sig]
            aide_Hessian[t,ind+1] <- -sqrt(temper)*derivative[t,ind+1]/std_sigma[etat_sig]
            for (  z in 1 : AR_lags ) {
                if(t-z>0) {
                    derivative[t,ind+1+z] <- -y[t-z]*(etat==q) - theta[ind2+1+AR_lags+1]*derivative[t-1,ind+1+z]
                    Jacobien[t,ind+1+z] <-  -temper*eps_t[t]*derivative[t,ind+1+z]/sigma[etat_sig]
                    aide_Hessian[t,ind+1+z] <- -sqrt(temper)*derivative[t,ind+1+z]/std_sigma[etat_sig]
                   } #
               } #
            derivative[t,ind+1+AR_lags+1] <- -eps_t[t-1]*(etat==q) - theta[ind2+1+AR_lags+1]*derivative[t-1,ind+1+AR_lags+1]
            Jacobien[t,ind+1+AR_lags+1] <-  -temper*eps_t[t]*derivative[t,ind+1+AR_lags+1]/sigma[etat_sig]
            aide_Hessian[t,ind+1+AR_lags+1] <- -sqrt(temper)*derivative[t,ind+1+AR_lags+1]/std_sigma[etat_sig]


           } #
       } #
    mat_prior <- kronecker(diag(1,regime),inv_Sig_prior)
    gradient <- as.matrix(colSums(Jacobien)) + mat_prior%*%theta[1:(taille_ARMA*regime)]
    Hessian_approx <- t(aide_Hessian)%*%aide_Hessian + mat_prior
##   } #

   returnlist <- list(gradient, Hessian_approx, Jacobien)
   return(returnlist)
 }


