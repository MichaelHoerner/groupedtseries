#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export


PSO_ARMA_sn <-function( y,X,regime_max,regime_sig_max,AR_lags,MA_lags,sn,sn_sig,nb_iter,nb_part, ts_segments) {


dimension <- dim(as.data.frame(y))
taille <- dimension[1]
number_ts_segments <- nrow(ts_segments)

sn <- as.matrix(sn)
sn_sig <- as.matrix(sn_sig)

dim_sn <- dim(sn)
# if(taille!=dim_sn[1]) {
#     sn <- t(sn)
#     sn_sig <- t(sn_sig)
#     dim_sn <- dim(as.data.frame(sn))
#    } #
if(dim_sn[2]!=1) {
    multi_sn <- 1
} else {
    multi_sn <- 0
} #

taille_ARMA <- AR_lags + MA_lags +1
# for (  q in 1 : regime ) {
#     var[ q].y <- matrix(0,sum(n_ii[,q] ] ,1)
#     var[ q].X <- matrix(0,sum(n_ii[,q] ] ,AR_lags+1)
#    } #
# count <- matrix(0,regime,1)
# for (  t in 1 : taille ) {
#     etat <- mode_sn[t]
#     count[etat] <- count[etat] +1
#     var[ etat].y[count[etat]] <- y[t ]
#     var[ etat].X[count[etat],] <- X[t,1:AR_lags+1 ]
#    } #
# for (  q in 1 : regime ) {
#     var[ q].beta <- solve(var[q].X"*var[q].X]*var[q].X" *var[q )  .y
#     var[ q].sigma <- mean( [var[q].y-var[q].X*var[q].beta].*[var[q].y-var[q].X*var[q].beta] ]
#    } #

theta <- matrix(0,taille_ARMA*regime_max,nb_part)
sigma <- matrix(0,regime_sig_max,nb_part)
best_chain <- matrix(0,nb_part,1)
best_theta_chain<- matrix(0,taille_ARMA*regime_max,nb_part)
best_sigma_chain<- matrix(0,regime_sig_max,nb_part)
best_theta <- matrix(0,taille_ARMA*regime_max,1)
best_sigma <- matrix(0,regime_sig_max,1)
vel <- matrix(0,taille_ARMA*regime_max,nb_part)
vel_sigma <- matrix(0,regime_sig_max,nb_part)
mu_theta <- matrix(0,taille_ARMA,1)
Sigma_theta <- diag(1,taille_ARMA)*0.25
best_dens <- 0
a <- 1.5
b <- 1.5


for (  i in 1 : nb_iter ) {
    w <- 0.9-0.7*(i/nb_iter)
    for (  j in 1 : nb_part ) {
        if(i==1) {
#             if(multi_sn==1) {
#                 it <- ceiling[dim_sn[1,2]*runif(1, 0, 1) ]
#                 mode_sn <- sn[,it]
#                 n_ii <- comtage(taille,mode_sn,regime)
#                 for (  q in 1 : regime ) {
#                     var[ q].y <- matrix(0,sum(n_ii[,q] ] ,1)
#                     var[ q].X <- matrix(0,sum(n_ii[,q] ] ,AR_lags+1)
#                    } #
#                 count <- matrix(0,regime,1)
#                 for (  t in 1 : taille ) {
#                     etat <- mode_sn[t]
#                     count[etat] <- count[etat] +1
#                     var[ etat].y[count[etat]] <- y[t ]
#                     var[ etat].X[count[etat],] <- X[t,1:AR_lags+1 ]
#                    } #
#                 for (  q in 1 : regime ) {
#                     var[ q].beta <- solve(var[q].X"*var[q].X]*var[q].X" *var[q )  .y
#                     var[ q].sigma <- mean( [var[q].y-var[q].X*var[q].beta].*[var[q].y-var[q].X*var[q].beta] ]
#                    } #
#                } #

            for (  q in 1 : regime_max ) {
                test <- 0
                while(test==0) {
                    theta_prop <- t(mvrnorm(1, mu_theta,Sigma_theta)) ## Normal multivari? RESTRICTION A METTRE ?
                    test <- 1
                    if(AR_lags>0) {
                        root_AR <- polyroot(rev(c(1, -t(theta_prop[2:(1+AR_lags)]))))
                        if(sum(abs( Re(root_AR))>=1)!=0 || any(abs(theta_prop[2:(1+AR_lags)])>=1)) {
                            test <- 0
                           } #
                       } #
                    if(MA_lags>0) {
                        root_MA <- polyroot(rev(c(1, -t(theta_prop[(2+AR_lags):taille_ARMA]))))
                        if(sum(abs( Re( root_MA))>=1)!=0 || any(abs(theta_prop[(2+AR_lags):taille_ARMA])>=1)) {
                            test <- 0
                           } #
                       } #
                    if(test == 1) {
                        theta[((q-1)*taille_ARMA+1):(q*taille_ARMA),j] <- theta_prop
                       } #
                   } #

                while(sigma[q,j]<=0) {
                    sigma[q,j] <- runif(1, 0, 1)*2#normrnd[var[ q].sigma,0.1 ]
                   } #
               } #

            for (  q in 1 : regime_sig_max ) {
                while(sigma[q,j]<=0) {
                    sigma[q,j] <- runif(1, 0, 1) *2#normrnd[var[ q].sigma,0.1 ]
                   } #
               } #

            if(multi_sn==0) {
             dens <- 0
             eps_t <- double(length=taille)
             dens <- .C("dens_CP_ARMA_c",
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
                 eps_out = as.vector(eps_t, mode="double"),
                 PACKAGE="groupedtseries")$dens

            } else {
                dens <- 0
                eps_t <- double(length=taille)
                for (  z in 1 : dim_sn[2] ) {
                    val_dens <- .C("dens_CP_ARMA_c",
                                 as.vector(y, mode = "double"),
                                 as.vector(X),
                                 as.integer(AR_lags),
                                 as.integer(MA_lags),
                                 as.vector(theta[,j], mode = "double"),
                                 as.vector(sigma[,j], mode = "double"),
                                 as.vector(sn[,z], mode = "double"),
                                 as.vector(sn_sig[,z], mode = "double"),
                                 as.integer(1),
                                 as.integer(taille),
                                 as.integer(number_ts_segments),
                                 as.vector(ts_segments$index_AR_start, mode="integer"),
                                 as.vector(ts_segments$index_AR_end, mode="integer"),
                                 as.vector(ts_segments$length_AR, mode="integer"),
                                 dens = as.double(dens),
                                 eps_out = as.vector(eps_t, mode="double"),
                                 PACKAGE="groupedtseries")$dens

                    dens <- dens + val_dens
                   } #
                dens <- dens/dim_sn[2]
               } #


            best_chain[j] <- dens
            best_theta_chain[,j] <- theta[,j]
            best_sigma_chain[,j] <- sigma[,j]
            if(j==1) {
                best_dens <- dens
                best_theta <- theta[,j]
                best_sigma <- sigma[,j]
            } else if(dens>best_dens) {
                best_dens <- dens
                best_theta <- theta[,j]
                best_sigma <- sigma[,j]
               } #
        } else {
            vel[,j] <- w*vel[,j] + a*runif(regime_max*taille_ARMA,0,1)*(best_theta-theta[,j]) + b*runif(regime_max*taille_ARMA,0,1)*(best_theta_chain[,j]-theta[,j])
            vel_sigma[,j] <- w*vel_sigma[,j] + a*runif(regime_sig_max,0,1)*(best_sigma-sigma[,j]) + b*runif(regime_sig_max,0,1)*(best_sigma_chain[,j]-sigma[,j])

            theta_prop_all <- theta[,j] + vel[,j]
            for (  q in 1 : regime_max ) {
                theta_prop <- theta_prop_all[((q-1)*taille_ARMA+1):(q*taille_ARMA)]
                test <- 1
                if(AR_lags>0) {
                    root_AR <- polyroot(rev(c(1, -t(theta_prop[2:(1+AR_lags)]))))
                    if(sum(abs(Re(root_AR))>=1)!=0 || any(abs(theta_prop[2:(1+AR_lags)])>=1)) {
                        test <- 0
                       } #
                   } #
                if(MA_lags>0) {
                    root_MA <- polyroot(rev(c(1, -t(theta_prop[(2+AR_lags):taille_ARMA]))))
                    if(sum(abs( Re( root_MA))>=1)!=0 || any(abs(theta_prop[(2+AR_lags):taille_ARMA])>=1)) {
                        test <- 0
                       } #
                   } #
                if(test == 0) {
                    theta_prop_all[((q-1)*taille_ARMA+1):(q*taille_ARMA)] <- theta[((q-1)*taille_ARMA+1):(q*taille_ARMA),j]
                   } #
               } #

            sigma_prop <- sigma[,j] + vel_sigma[,j]
            test <- 1
            for ( q in 1 : regime_sig_max ) {
                if(sigma_prop[q]<0) {
                    test <- 0
                   } #
               } #
            if(test==1) {
                theta[,j] <- theta_prop_all
                sigma[,j] <- sigma_prop
                dens <- 0
                eps_t <- double(length=taille)
                if(multi_sn==0) {
                    dens <- .C("dens_CP_ARMA_c",
                             as.vector(y, mode = "double"),
                             as.vector(X),
                             as.integer(AR_lags),
                             as.integer(MA_lags),
                             as.vector(theta[,j], mode = "double"),
                             as.vector(sigma[,j], mode = "double"),
                             as.vector(sn, mode = "double"),
                             as.vector(sn_sig, mode = "double"),
                             as.integer(1),
                             as.integer(taille),
                             as.integer(number_ts_segments),
                             as.vector(ts_segments$index_AR_start, mode="integer"),
                             as.vector(ts_segments$index_AR_end, mode="integer"),
                             as.vector(ts_segments$length_AR, mode="integer"),
                             dens = as.double(dens),
                             eps_out = as.vector(eps_t, mode="double"),
                             PACKAGE="groupedtseries")$dens
                } else {
                    dens <- 0
                    eps_t <- double(length=taille)
                    for ( z in 1 : dim_sn[2] ) {
                        val_dens <- .C("dens_CP_ARMA_c",
                                       as.vector(y, mode = "double"),
                                       as.vector(X),
                                       as.integer(AR_lags),
                                       as.integer(MA_lags),
                                       as.vector(theta[,j], mode = "double"),
                                       as.vector(sigma[,j], mode = "double"),
                                       as.vector(sn[,z], mode = "double"),
                                       as.vector(sn_sig[,z], mode = "double"),
                                       as.integer(1),
                                       as.integer(taille),
                                       as.integer(number_ts_segments),
                                       as.vector(ts_segments$index_AR_start, mode="integer"),
                                       as.vector(ts_segments$index_AR_end, mode="integer"),
                                       as.vector(ts_segments$length_AR, mode="integer"),
                                       dens = as.double(dens),
                                       eps_out = as.vector(eps_t, mode="double"),
                                       PACKAGE="groupedtseries")$dens

                        dens <- dens + val_dens
                       } #
                    dens <- dens/dim_sn[2]
                   } #

                if(is.na(dens)==0 &&  is.infinite(dens)==0) {
                    if(dens>best_chain[j]) {
                        best_chain[j] <- dens
                        best_theta_chain[,j] <- theta[,j]
                        best_sigma_chain[,j] <- sigma[,j]
                        if(dens>best_dens) {
                            best_theta <- theta[,j]
                            best_sigma <- sigma[,j]
                            best_dens <- dens
                           } #
                       } #
                   } #
               } #
           } #
       } #
   } #


returnlist <- list(best_dens, best_theta, best_sigma)
return(returnlist)

}

