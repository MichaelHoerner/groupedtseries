#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


 multinomialrnd <-function( vect ) {

unif <- runif(1, 0, 1)
ind <- 0
Q <- 0
while(Q<unif) {
    ind <- ind+1
    Q <- Q + vect[ind]
   } #

return(ind)
}
