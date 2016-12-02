#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

 prior <-function( param,regime,n_ii ) {

## sampling alpha + kap
v <- matrix(0,regime,1)
w <- matrix(0,regime,1)
for (  i in 1 : regime ) {
    n <- sum(n_ii[i,])
    if(runif(1, 0, 1) < n/(n+param$alpha+param$kappa)) {
        v[i] <- 1
       } #
    w[i] <- rbeta(1, param$alpha+param$kappa+1,n)
   } #
m <- sum(sum(param$m))
AlphaKapa <- rgamma(1, param$prior_alpha_kappa_a + m - sum(v),1/(1/param$prior_alpha_kappa_b-sum(log(w))) )

## Sampling rho

rho <- rbeta(1, param$prior_rho_a + sum(param$r),param$prior_rho_b+m-sum(param$r))

## calculate alpha and kapa
alpha <- (1-rho)*AlphaKapa
kappa  <- rho*AlphaKapa

## Sample etha
marg <- matrix(0,regime,1)
indice <- matrix(0,regime,1)
for (  i in 1 : regime ) {
    marg[i] <- sum(param$m_bar[,i])
    if(marg[i]>0) {
        indice[i] <- 1
    } else {
        indice[i] <- 0
       } #
   } #
m_bar <- sum(marg)
if(runif(1, 0, 1) <(m_bar/(m_bar+param$etha))) {
    t <- 1
} else {
    t <- 0
   } #
lambda <- rbeta(1, param$etha+1,m_bar)


K_bar <- sum(indice)
etha <- rgamma(1, param$prior_eta_a + K_bar - t,1/(1/param$prior_eta_b-log(lambda)) )

returnlist <- list(alpha, kappa, etha)

return(returnlist)
}

