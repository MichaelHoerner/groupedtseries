#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

launch_FB_Klaassen <- function( y,X,regime,theta,sigma,sn,sn_sig,p_mat,AR_lags,MA_lags,deb,fin,temper) {



accept <- 0
lags_q <- 1
sigma_t <- sigma[sn_sig]

taille <- max(dim(as.data.frame(y)))
#print(taille)
sn_move <- matrix(0, nrow=taille, ncol=1)
log_like <- double(length=1)
log_q_prior <- double(length=4)
forward_out <- matrix(0, nrow=(fin-deb+1), ncol=regime)
X_MA <- matrix(0, nrow=(fin+1), ncol=regime)

if(lags_q == 1) {

    # ###################################DEBUG#################################
    # debug_dataframe <- data.frame(X_AR = X, Regime = regime)
    # debug_count <- 120
    # xtable_debug_dataframe <- xtable(debug_dataframe)
    # debug_filename <- paste0("debug/debug_data_", debug_count, ".csv")
    # write.csv2(xtable_debug_dataframe, file = debug_filename)
    # #########################################################################
    #
    # ###################################DEBUG#################################
    # debug_dataframe <- data.frame(theta = theta)
    # debug_count <- 121
    # xtable_debug_dataframe <- xtable(debug_dataframe)
    # debug_filename <- paste0("debug/debug_data_", debug_count, ".csv")
    # write.csv2(xtable_debug_dataframe, file = debug_filename)
    # #########################################################################
    #
    # ###################################DEBUG#################################
    # debug_dataframe <- data.frame(sn = sn)
    # debug_count <- 1201
    # xtable_debug_dataframe <- xtable(debug_dataframe)
    # debug_filename <- paste0("debug/debug_data_", debug_count, ".csv")
    # write.csv2(xtable_debug_dataframe, file = debug_filename)
    # #########################################################################



    templist <- .C("FB_Klaassen",
                   as.vector(y, mode = "double"),
                   as.integer(regime),
                   as.vector(sigma_t, mode="double"),
                   as.vector(p_mat, mode = "double"),
                   as.vector(X, mode = "double"),
                   as.vector(theta, mode = "double"),
                   as.vector(sn-1, mode = "double"),
                   as.integer(AR_lags),
                   as.integer(MA_lags),
                   as.integer(deb),
                   as.integer(fin),
                   as.integer(1),
                   as.vector(runif(fin-deb+1,0,1), mode="double"),
                   as.double(temper),
                   as.integer(taille),
                   sn_move = as.vector(sn_move, mode="double"),
                   log_like_out = as.double(log_like),
                   log_q_prior_out = as.vector(log_q_prior, mode="double"),
                   forward_out = as.vector(forward_out, mode="double"),
                   epsilon_Haas = as.vector(X_MA, mode="double"),
                   PACKAGE="groupedtseries")

    sn_move <- templist$sn_move
    log_q_prior <- templist$log_q_prior_out
    forward_out <- templist$forward_out
    forward_out <- matrix(forward_out, nrow=(fin-deb+1), ncol=regime)
    X_MA <- templist$epsilon_Haas
    X_MA <- matrix(X_MA, nrow = (fin+1), ncol = regime)

    # ###################################DEBUG#################################
    # debug_dataframe <- data.frame(log_q_prior = log_q_prior)
    # debug_count <- 122
    # xtable_debug_dataframe <- xtable(debug_dataframe)
    # debug_filename <- paste0("debug/debug_data_", debug_count, ".csv")
    # write.csv2(xtable_debug_dataframe, file = debug_filename)
    # #########################################################################
    #
    # ###################################DEBUG#################################
    # debug_dataframe <- data.frame(forward_out = forward_out)
    # debug_count <- 123
    # xtable_debug_dataframe <- xtable(debug_dataframe)
    # debug_filename <- paste0("debug/debug_data_", debug_count, ".csv")
    # write.csv2(xtable_debug_dataframe, file = debug_filename)
    # #########################################################################
    #
    # ###################################DEBUG#################################
    # debug_dataframe <- data.frame(X_MA = X_MA)
    # debug_count <- 124
    # xtable_debug_dataframe <- xtable(debug_dataframe)
    # debug_filename <- paste0("debug/debug_data_", debug_count, ".csv")
    # write.csv2(xtable_debug_dataframe, file = debug_filename)
    # #########################################################################

    } else {
    # templist <- FB_Maciej(y,regime,sigma_t,p_mat,X,theta,sn-1,AR_lags,MA_lags,deb,fin,1,rand[fin-deb+1,1],temper,lags_q)

    } #

sn_move <- round(sn_move)+1

if(MA_lags==0) {
    sn <- sn_move
    accept <- 1
} else {
    dens <- 0
    eps_out <- double(length=taille)
    dens_move <- .C("dens_CP_ARMA_c",
                    as.vector(y, mode = "double"),
                    as.vector(X),
                    as.integer(AR_lags),
                    as.integer(MA_lags),
                    as.vector(theta, mode = "double"),
                    as.vector(sigma, mode = "double"),
                    as.vector(sn_move, mode = "double"),
                    as.vector(sn_sig, mode = "double"),
                    as.integer(deb),
                    as.integer(taille),
                    as.integer(1),
                    as.integer(1),
                    as.integer(taille),
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
                    as.integer(1),
                    as.integer(1),
                    as.integer(taille),
                    as.integer(taille),
                    dens = as.double(dens),
                    eps_out = as.vector(eps_t, mode="double"),
                    PACKAGE="groupedtseries")$dens

  if (any(is.na(log_q_prior)) ){ #|| any(is.infinite(log_q_prior))
    ##print(log_q_prior)
    ##print(dens_stay)
    ##print(dens_move)
    warning("log_q_prior shows an na-value")
  } else {

    if(exp(temper*(dens_move-dens_stay)+log_q_prior[2]+log_q_prior[3]-log_q_prior[1]-log_q_prior[4])>runif(1, 0, 1) ) {
        sn <- sn_move
        accept <- 1
       }
    }#
   } #

returnlist <- list(sn, accept, forward_out, X_MA)

return(returnlist)
}
