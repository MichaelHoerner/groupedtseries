#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


dirichlet_sample <-function( alpha ) {
# /*
# Fonction qui tire al?atoirement un ?l?ment dans une distribution dirichlet
#     Input :
#           alpha : repr?sente le vecteur de param?tre de la dirichlet D[alpha]
#           taille : dimension de la distribution dirichlet
#           al?atoire : vecteur de taille "size' qui contient des ?l?ments al?atoires tir?s d" une uniforme [0,1[
#     Output :
#            Vecteur al?atoire d'une distribution dirichlet
# */
dimension <- dim(as.data.frame(alpha))
dimension <- max(dimension)

sample <- matrix(0,dimension,1)
Q <- 0
for (  i in 1 : dimension ) {
    if(alpha[i]!=0) {
        sample[i] <- rgamma(1, alpha[i],1)
    } else {
        sample[i] <- 0
       } #
    Q <- Q + sample[i]
   } #
sample <- sample/Q

return(sample)

}


