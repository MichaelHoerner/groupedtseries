#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

Haas_FB_ARMA <-function( y,X,taille,regime,theta,sigma,sn,sn_sig,p,AR_lags,MA_lags,deb,fin,temper ) {

val_sigma <- sigma[sn_sig]
taille_ARMA <- 1+AR_lags+MA_lags
forward <- matrix(0,taille,regime)
eps_prev <- matrix(0,regime,1)
X_aide <- X
eps <- matrix(0,regime,1)
p_aide <- p[1:regime,1:regime]

err_out <- double(length=taille)
X_out <- double(length=taille_ARMA*taille)

eps_Haas <- matrix(0,taille,regime)

for (  q in 1 : regime ) {

    eps_reg <- .C("create_X_MA_term_c",
                  as.vector(y, mode = "double"),
                  as.vector(X, mode = "double"),
                  as.integer(AR_lags),
                  as.integer(MA_lags),
                  as.vector(theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)], mode = "double"),
                  as.vector(matrix(1,taille,1), mode = "double"),
                  as.integer(taille),
                  err_out = as.vector(err_out, mode="double"),
                  X_out = as.vector(X_out, mode="double"),
                  PACKAGE="groupedtseries")$err_out

    eps_Haas[2:nrow(eps_Haas),q] <- eps_reg[1:(length(eps_reg)-1)]
   } #


for (  i in deb : fin ) {
    constante <- 0


    #p_aide <- compute_p_aide_GRAY[p,u_t[i],slice,k,k_1]

   if(i==deb) {
       if(i==1) {
           a <- matrix(0,regime,1)
           a[sn[i]] <- 1
           transition <- t(a)%*%p_aide
           for (  q in 1 : regime ) {
               eps[q] <- y[i]- t(as.matrix(X_aide[i,]))%*%as.matrix(theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)])

               if(transition[q]!=0) {
                   forward[i,q] <- temper*(-0.5*log(2*pi*val_sigma[i]) - 0.5*(eps[q])^2/val_sigma[i]) + log(transition[q])# + log[p[q,q]]
                   if(constante==0) {
                       constante <- forward[i,q]
                   } else if(constante<forward[i,q]) {
                       constante <- forward[i,q]
                      } #
                  } #


              } #
       } else {
           a <- matrix(0,regime,1)
           a[sn[i-1]] <- 1
           transition <- t(a)%*%p_aide
           for (  q in 1 : regime ) {
               eps[q] <- y[i]-t(as.matrix(X_aide[i,]))%*%as.matrix(theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)])
               eps_prev[q] <- X_aide[i,1+AR_lags+1]
               if(transition[q]!=0) {
                   forward[i,q] <- temper*(-0.5*log(2*pi*val_sigma[i]) - 0.5*(eps[q])^2/val_sigma[i]) + log(transition[q])
                   if(constante==0) {
                       constante <- forward[i,q]
                   } else if(constante<forward[i,q]) {
                       constante <- forward[i,q]
                      } #
                  } #

              } #
          } #

   } else {
       transition <- t(as.matrix(forward[i-1,]))%*%p_aide

       for (  q in 1 : regime ) {
           if(transition[q]!=0) {
               X_aide[i,1+AR_lags+1] <- eps_Haas[i,q]


               if(MA_lags>1) {
                   if(i-1<MA_lags-1) {
                       iter <- i-1
                   } else {
                       iter <- MA_lags-1
                      } #
                   for (  z in 1 : iter ) {
                       X_aide[i,1+AR_lags+1+z] <- eps_Haas[i-z,q]
                      } #
                  } #


               eps[q] <- y[i]-t(as.matrix(X_aide[i,]))%*%as.matrix(theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)])

               forward[i,q] <- (temper*(-0.5*log(2*pi*val_sigma[i]) - 0.5*(eps[q])^2/val_sigma[i]) + log(transition[q]))
               if(constante==0) {
                   constante <- forward[i,q]
               } else if(constante<forward[i,q]) {
                   constante <- forward[i,q]
                  } #
              } #


          } #
      } #



   for (  q in 1 : regime ) {
       if(forward[i,q]!=0) {
         forward[i,q] <- exp(forward[i,q]-constante)
          } #
      } #
   forward[i,] <- forward[i,]/sum(forward[i,])
   } #
# if(regime>1) {
#      plot[forward]
#      theta'
#      sigma
#      pause
#    } #

sn_new <- sn
prior_stay <- 0
prior_move <- 0
q_move <- 0
q_stay <- 0
back <- matrix(0,regime,1)

if(fin==taille) {
    sn_new[taille] <- multinomialrnd(forward[taille,])
    q_move <- log(forward[taille,sn_new[taille]])
    q_stay <- log(forward[taille,sn[taille]])
    incr <- fin-1
} else {
    incr <- fin-1
#     for (  q in 1 : regime ) {
#         for (  z in 1 : regime ) {
#             #p_aide[q,z] <- [p[q,z]>u_t[i+1]]*[[u_t[i+1]<slice*p[q,z]]*k+[u_t[i+1]>=slice*p[q,z]]*k_1]
#             if(p[q,z]>u_t[i+1]) {
#                 if(u_t[fin+1]<slice*p[q,z]) {
#                     p_aide[q,z] <- k
#                 } else {
#                     p_aide[q,z] <- k_1
#                    } #
#                 #                 p_aide[q,z] <- 1#/(1-p[q,z])
#             } else {
#                 p_aide[q,z] <- 0
#                } #
#            } #
#        } #
#     etat <- sn_new[fin+1]
#    for (  q in 1 : regime ) {
#       back[q] <- forward[i,q]*p_aide[q,etat]
#       } #
#    back <- back/sum(back]
#    etat_current <-  multinomialrnd[back]
#
#    sn_new[fin] <- etat_current
#    if(MA_lags>0) {
#        q_move <- log(back[etat_current]]#+ log[p[sn_new[fin],sn_new[fin]]]
#        prior_move <- prior_move + log(p_aide[etat_current,etat]]
#
#        for (  q in 1 : regime ) {
#            back[q] <- forward[fin,q]*p_aide[q,sn[fin+1]]
#           } #
#        back <- back/sum(back]
#
#        q_stay <-  log(back[sn[fin]]]# + log[p[sn[fin],sn[fin+1]]]
#        prior_stay <- prior_stay + log(p_aide[sn[fin],sn[fin+1]]]
#       } #

   } #


 for (  i in seq(from=incr, to=deb, by=-1 ) ) {
    #p_aide <- compute_p_aide_GRAY[p,u_t[i+1],slice,k,k_1]
#     p_aide <- [p>u_t[i+1]].*[[u_t[i+1]<slice*p]*k+[u_t[i+1]>=slice*p]*k_1]
#    for (  q in 1 : regime ) {
#         for (  z in 1 : regime ) {
#             p_aide[q,z] <- [p[q,z]>u_t[i+1]]*[[u_t[i+1]<slice*p[q,z]]*k+[u_t[i+1]>=slice*p[q,z]]*k_1]
# #             if(p[q,z]>u_t[i+1]) {
# #                 if(u_t[i+1]<slice*p[q,z]) {
# #                     p_aide[q,z] <- k
# #                 } else {
# #                     p_aide[q,z] <- k_1
# #                    } #
# # #                 p_aide[q,z] <- 1#/(1-p[q,z])
# #             } else {
# #                 p_aide[q,z] <- 0
# #                } #
#            } #
#        } #

   etat <- sn_new[i+1]

   back <- as.matrix(forward[i,])*as.matrix(p_aide[,etat])

#    for (  q in 1 : regime ) {
#       back[q] <- forward[i,q]*p_aide[q,etat]
#
#       } #
#    back
#    pause
   back <- back/sum(back)
   etat_current <-  multinomialrnd(back)

   sn_new[i] <- etat_current


       q_move <- q_move  + log(back[etat_current])
       prior_move <- prior_move + log(p_aide[etat_current,etat])
#       prior_move <- prior_move + log(p[etat_current,etat]]

      back <- as.matrix(forward[i,])*as.matrix(p_aide[,sn[i+1]])
#
#        for (  q in 1 : regime ) {
#            back[q] <- forward[i,q]*p_aide[q,sn[i+1]]
#           } #
       back <- back/sum(back)

       q_stay <- q_stay  + log(back[sn[i]])
       prior_stay <- prior_stay + log(p_aide[sn[i],sn[i+1]])
#       prior_stay <- prior_stay + log(p[sn[i],sn[i+1]]]


   } #

if(deb!=1) {
    prior_stay <- prior_stay + log(p_aide[sn[deb-1],sn[deb]])
    prior_move <- prior_move + log(p_aide[sn_new[deb-1],sn_new[deb]])
   } #


if(MA_lags>0) {
  dens <- 0
  eps_out <- double(length=taille)

    dens_move <- .C("dens_CP_ARMA_c",
                    as.vector(y, mode = "double"),
                    as.vector(X),
                    as.integer(AR_lags),
                    as.integer(MA_lags),
                    as.vector(theta, mode = "double"),
                    as.vector(sigma, mode = "double"),
                    as.vector(sn_new, mode = "double"),
                    as.vector(sn_sig, mode = "double"),
                    as.integer(deb),
                    as.integer(taille),
                    dens = as.double(dens),
                    eps_out = as.vector(eps_t, mode="double"),
                    PACKAGE="groupedtseries")$dens

    dens_stay <- .C("dens_CP_ARMA_c",
                    as.vector(y, mode = "double"),
                    as.vector(X),
                    as.integer(AR_lags),
                    as.integer(MA_lags),
                    as.vector(theta, mode = "double"),
                    as.vector(sigma, mode = "double"),
                    as.vector(sn, mode = "double"),
                    as.vector(sn_sig, mode = "double"),
                    as.integer(deb),
                    as.integer(taille),
                    dens = as.double(dens),
                    eps_out = as.vector(eps_t, mode="double"),
                    PACKAGE="groupedtseries")$dens


    prob <- exp(dens_move+prior_move+q_stay-dens_stay-prior_stay-q_move)
    if(runif(1, 0, 1) <prob) {
        sn <- sn_new
        accept <- 1
    } else {
        accept <- 0
       } #
} else {
    sn <- sn_new
    accept <- 1
   } #

 returnlist <- list(sn, accept, forward)
 return(returnlist)
 }



