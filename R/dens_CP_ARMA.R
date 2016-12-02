 dens_CP_ARMA <-function( y,X,AR_lags,MA_lags,theta,sigma,tau,taille ) {
taille_ARMA <- 1+AR_lags+MA_lags
ind  <-1
log_dens <- 0
for (  i in 1 : (taille-1) ) {
    if(i>tau[ind]) {
        ind <- ind +1
       } #
    log_dens <- log_dens-0.5*log(2*pi*sigma[ind]) - 0.5*(y[i,1]-X[i,]%*%theta[((ind-1)*taille_ARMA+1):((ind)*taille_ARMA)])^2/sigma[ind]
 } #


return(log_dens)
}
