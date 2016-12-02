#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

create_X_AR <- function(y,taille,AR_lags,MA_lags) {


X_AR <- matrix(0,taille-AR_lags,1+AR_lags+MA_lags)
X_AR[,1] <- matrix(1,taille-AR_lags,1)

y_AR <- vector(length = taille-1, mode = "double")

y_AR[1:length(y_AR)] <- y[(AR_lags+1):max(dim(as.data.frame(y))),1]

for (  i in 1 : AR_lags ) {
    X_AR[,(AR_lags+2-i)] <- y[i:(taille-AR_lags+i-1),1]
   } #

returnvalues <- list(y_AR, X_AR)

return(returnvalues)

}

