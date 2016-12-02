#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


calc_mode_sn <- function( post_sn ) {

dimension <- dim(post_sn)
regime <- max( max( post_sn))
sn_mode <- matrix(0,dimension[1],1)
for (  t in 1 : dimension[1] ) {
    count <- matrix(0,regime,1)
   for (  j in 1 : dimension[2] ) {
       etat <- post_sn[t,j]
       count[etat] <- count[etat] +1
      } #
   sn_mode[t] = which.max(count)

   } #

return(sn_mode)
}
