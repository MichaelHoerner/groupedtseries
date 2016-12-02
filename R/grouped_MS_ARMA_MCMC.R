#' Fit grouped GARCH model to time series
#'
#' Fit a Generalized Autoregressive Conditional Heteroscedastic GARCH(p,q)
#' time series model to the data by computing the maximum-likelihood
#' estimates of the conditionally normal model
#' to a group of univariate time series by finding coefficients
#' that serve best the whole group.
#' @param x a two-dimensional numeric vector of a grouped time series with accompanying group index in second dimension.
#' One-dimensional if input shall be only a single time series.
#' @param order a vector c(p,q) that indicates the order p of ARCH and q of GARCH part
#' @param series name for the series. Defaults to \code{deparse(substitute(x))}
#' @param control a list of control parameters as set up by \code{garch.control}
#' @param maxiter gives the maximum number of log-likelihood function
#' evaluations \code{maxiter} and the maximum number of iterations
#' \code{2*maxiter} the optimizer is allowed to compute.
#' @param trace logical. Trace optimizer output?
#' @param start If given this numeric vector is used as the initial estimate
#' of the GARCH coefficients.  Default initialization is to set the
#' GARCH parameters to slightly positive values and to initialize the
#' intercept such that the unconditional variance of the initial GARCH
#' is equal to the variance of \code{x}
#' @param grad character indicating whether analytical gradients or a numerical
#' approximation is used for the optimization.
#' @param abstol absolute function convergence tolerance.
#' @param reltol relative function convergence tolerance.
#' @param xtol coefficient-convergence tolerance.
#' @param falsetol false convergence tolerance.
#' @param \dots additional arguments for \code{\link{qr}} when computing
#' the asymptotic covariance matrix.
#'
#' @details \code{garch} uses a Quasi-Newton optimizer to find the maximum
#' likelihood estimates of the conditionally normal model.  The first
#' max(p, q) values are assumed to be fixed.  The optimizer uses a hessian
#' approximation computed from the BFGS update.  Only a Cholesky factor
#' of the Hessian approximation is stored.  For more details see Dennis
#' et al. (1981), Dennis and Mei (1979), Dennis and More (1977), and
#' Goldfarb (1976).  The gradient is either computed analytically or
#' using a numerical approximation.
#' @return Object of class \code{GARCH}.
#' \item{order}{the order of the fitted model.}
#' \item{coef}{coef estimated GARCH coefficients for the fitted model across all grouped time series.}
#' \item{n.likeli}{the negative log-likelihood function evaluated at the
#' coefficient estimates (apart from some constant).}
#' \item{n.used}{the number of observations of \code{x}.}
#' \item{residuals}{the series of residuals.}
#' \item{fitted.values}{the bivariate series of conditional standard
#' deviation predictions for \code{x}.}
#' \item{series}{the name of the series \code{x}.}
#' \item{frequency}{the frequency of the series \code{x}.}
#' \item{call}{the call of the \code{garch} function.}
#' \item{vcov}{outer product of gradient estimate of the asymptotic-theory
#' covariance matrix for the coefficient estimates.}
#'
#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

MS_ARMA_MCMC <- function(y,Expl_X=NULL,nb_MCMC=10000,regime=1,AR_lags=1,MA_lags=1,CP_or_MS="MS",temper=1,force_reg_1=0,sn_dep=0,sn_sig_dep=0,theta_dep=0,sigma_dep=0) {
#################################################################################################################?
##### Function that draws realizations of the posterior distribution of an
##### IHMM-ARMAX model with different break structures for the mean
##### parameters and the variance.
#################################################################################################################?
#################################################################################################################?
##### Inputs :
##### y - Dep   } #ent variable (dimension : T*1)
##### nb_MCMC - Number of MCMC iterations (default : 10000)
##### regime - Initial number of regimes (default : 1)
##### AR_lags - Number of autoregressive components (default : 1)
##### MA_lags - Number of lagged error terms (default : 1) - Cannot be
##### higher than 1.
##### temper - Useful for the Marginal likelihood computation (default : 1)
##### force_reg_1 - if <- 1, only one regime exists for the mean parameters
##### and the variance
##### sn_dep and sn_sig_dep - Initial values of state vectors (default :
##### none)
#################################################################################################################?
#############
##### Output :
##### Result : Struture that contains many information on the posterior
##### distribution.
##### result.post.* are related to the realizations of the parameters after
##### burn-in
##### result.theta.t stands for a 3D matrix with the mean parameter over
##### time per MCMC iteration
##### result.sigma.t is the obtained variance over time for each draw of
##### the MCMC.
##### result.mean.theta denotes the mean of result.theta.t. It is therefore
##### the posterior expectations of the mean parameters over time.
##### result.std.theta is the standard deviation of result.theta.t
##### result.mean.sigma and result.std.sigma are similar to mean.theta and
##### std.theta but for the variance
##### result.prob.reg summarizes the frequencies of observing a specific
##### number of regimes in the mean and in the variance.
#################################################################################################################?
#################################################################################################################?
##### Examples to run the program
##### [result] <- MS_ARMA_MCMC(y,[],10000,3,2,1,'MS') : Estimation of an
##### IHMM-ARMA(2,1) with 10000 MCMC iterations and with state vectors exhibiting 3 regimes as starting values
##### for the MCMC and with Markov-switching-type of hyperparameters.
##### [result] <- MS_ARMA_MCMC(y,[],10000,1,1,1,'MS',1,1) : Estimation of an ARMA
##### without any structural breaks.
##### [result] <- MS_ARMA_MCMC(y,X,10000,1,2,1,'MS',1) : Estimation of an
##### IHMM-ARMAX(2,1)
#################################################################################################################?
#################################################################################################################?
sample_all <- 0 ## If <- 1, the state vector is sampled in one pass (not recommended when MA_lags <- 1)
if(MA_lags==0) {
    sample_all<-1
   } #
sampling_sig <- 1 ## If <- 1, Indep   } #ent structural break in the variance
choix_Haas <- 0 ## If <- 1, sampling of the mean state vector by the approximate model of HMP instead of Klaassen.
dependance <- 0 ### Si on veut de la d?pendance entre les vecteurs ?tats (si la moyenne change, la variance change p-e aussi)
dist_sig <- 0
if(dependance==1) {
    if(max(size(y))<1000) {
        if(max(size(y))<400) {
            dist_sig <- 4
            }                 ### Horizon de dependance
        else {
            dist_sig <- 12 ### Horizon de dependance
           } #
    }
    else {
        dist_sig <- 25 ### Horizon de dependance
       } #
   } #

y <- as.data.frame(y)

display_graph <- 0 ## If <- 1, at the end of a run, a Figure that summarizes the posterior realizations appears.
nb_forecast <- 16 ## Number of forecast horizon
seed_id <- 0 ## if <-1, it fixes the seed in order to reproduce the results.
if(seed_id==1) {
    warning('SEEDS IDENTIQUES') ##ok<WNTAG>
    rnorm(100) ##ok<RAND>
    runif(100,0,1) ##ok<RAND>
   } #

###### test issues
Expl_X<-NULL
nb_MCMC<-10
regime<-1
AR_lags<-1
MA_lags<-1
CP_or_MS<-"MS"
temper<-1
force_reg_1<-0
sn_dep<-NULL
sn_sig_dep<-NULL
theta_dep<-NULL
sigma_dep<-NULL
###################


dep_MCMC <- 1
Rolling_MCMC <- 1
if(is.null(theta_dep) || is.null(sigma_dep)) {
    Rolling_MCMC <- 0
    if(is.null(sn_dep) || is.null(sn_sig_dep)) {
        dep_MCMC <- 0

       } #
   } #



if(force_reg_1!=0 && regime>1) {
    regime <- 1
    ##disp('Initial number of regime set to one as the model does not have a structure to handle breaks')
   } #

dimension <- dim(y)
taille <- dimension[1]
if(taille<dimension[2]) {
    taille <- dimension[2]
    y <- t(y)
   } #
##### y_AR : variable d?p   } #ante
### x_AR : variables explicatives de y_AR. x_AR <- [cst lags.de.y lags.erreur]
AR_lags.init <- AR_lags
templist <- create_X_AR(y,taille,AR_lags,MA_lags) ## OwnFunction ## Cr?ation des donn?es en prenant en compte le nombre de lag en AR (on ?limine toutes les valeurs jusqu'? max(p,q))

y_AR <- templist[[1]]
X_AR <- templist[[2]]

dimension <- dim(as.data.frame(y_AR))
taille <- dimension[1]
if(taille<dimension[2]) {
    taille <- dimension[2]
    y_AR <- t(y_AR)
   } #
taille_Expl_X <- dim(as.data.frame(Expl_X))
if(min(taille_Expl_X)!=0){
    if(taille_Expl_X[1]<taille_Expl_X[2] || taille_Expl_X[1]!=taille) {
        if(taille_Expl_X[1]-AR_lags==taille) {
            x_AR <- c(x_AR,Expl_X[(AR_lags+1):length(Expl_X),])
            AR_lags <- AR_lags + min(taille_Expl_X)
        } else {
            error('Dimension issue related to the explanatory variables')
           } #
    } else {
        x_AR <- c(x_AR,Expl_X)
        AR_lags <- AR_lags + min(taille_Expl_X)
       } #
   } #
taille_ARMA <- 1 + AR_lags + MA_lags ## taille_ARMA est une constante p   } #ant tout le programme et repr?sente le nombre d'inconnue dans la moyenne (dans un r?gime)

eps_t <- double(length=taille)

# if(Rolling_MCMC==0)
#     ##disp('********************************************************************')
#     if(taille_Expl_X>0)
#         ##disp(['*****************  ARMAX(' num2str(AR_lags.init) ',' num2str(MA_lags) ')'])
#     else
#         ##disp(['*****************  ARMA(' num2str(AR_lags) ',' num2str(MA_lags) ')'])
#        } #
#     if(force_reg_1==0)
#         ##disp('*****************  Including the IHMM structure into the model ')
#     else
#         ##disp('*****************  No break in the model ')
#        } #
#     if(dependance==1)
#         ##disp('*****************  Dep   } #ence in the state vectors ')
#     else
#         ##disp('*****************  No dep   } #ence in the state vectors ')
#        } #
#     switch CP_or_MS
#         case 'CP'
#             ##disp('*****************  IHMM hyper-parameters to behave like a Change-point model')
#         otherwise
#             ##disp('*****************  IHMM hyper-parameters are uninformative (MS behaviour) ')
#        } #
#     ##disp(['*****************  Time Series of ' num2str(taille) ' observations'])
#     ##disp(['*****************  Number of MCMC iterations : ' num2str(nb_MCMC)])
#     ##disp('********************************************************************')
#     ##disp('Finding Initial values for the MCMC...')
#    } #

###############
### Initilization of the hyper-parameters
###############
Struct_Prior <- define_prior(taille_ARMA,CP_or_MS) ## OwnFunction
V_prior_inv <- Struct_Prior$V_prior_inv
dv_prior <- Struct_Prior$dv_prior
mean_mu_prior_theta <- Struct_Prior$mean_mu_prior_theta
Sigma_mu_prior_theta_inv <- Struct_Prior$Sigma_mu_prior_theta_inv
Sigma_theta_inv <- Struct_Prior$Sigma_theta_inv
e_prior <- Struct_Prior$e_prior
f_prior <- Struct_Prior$f_prior
c_0 <- Struct_Prior$c_0
d_0 <- Struct_Prior$d_0
MH_e <- Struct_Prior$MH_e
rho_e <- Struct_Prior$rho_e
accept_e <- 0


###############
### Finding Initial values for the MCMC
###############
dens_chain <- matrix(0,nb_MCMC,1)
nb_init <- 1
if(regime>1) {
  if(regime>4 && taille>1000) {
    nb_init <- 60
  } else if (regime>4 && taille<1000) {
    nb_init <- 50 }
  else {
    nb_init <- 10
  } #
} #
 regime_sn_max <- 20 ## Maximum number of regimes in the mean
 regime_sn_sig_max <- 20 ## Maximum number of regimes in the variance
 templist <- initialisation_sn(nb_init,taille,regime,0) ## OwnFunction ## Generation of starting values of the mean state vector
 init_sn <- templist[[1]]
 templist <- initialisation_sn(nb_init,taille,regime,0) ## OwnFunction ## Generation of starting values of the variance state vector
 init_sn_sig <- templist[[1]]
 init_reg <- matrix(0,nb_init,2)
 dens <- matrix(0,nb_init,1)
 theta_init <- matrix(0,taille_ARMA*regime_sn_max,nb_init)
 sigma_init <- matrix(0,regime_sn_sig_max,nb_init)
 if(regime==1) {
  nb_part <- 20
  nb_iter <- 20
 } else {
  nb_part <- 40
  nb_iter <- 30
 } #

 ###############
 ### Soit recherche d'un point de d?part par particle swarm : Optimisation
 ### avec structural break ou bien on donne directement les valeurs de
 ### d?part
 ###############
 if(Rolling_MCMC==1) {
  # dens_chain[1] <- .C("dens_CP_ARMA_c", as.vector(y_AR), X_AR,
  #               AR_lags, MA_lags, theta_dep, sigma_dep, sn_dep,
  #               sn_sig_dep, 1, PACKAGE = "tseries")
  # dens_chain[1] <- dens_CP_ARMA(y_AR,X_AR,AR_lags,MA_lags,theta_dep,sigma_dep,sn_dep,sn_sig_dep) #OwnFunction
  # theta <- theta_dep
  # sigma <- sigma_dep
  # sn <- double(sn_dep)
  # sn_sig <- double(sn_sig_dep)
  # regime <-  max(sn)
 } else {
  if(regime==1) {
    #disp('Initialisation by Particle Swarm Optimization')
    templist <- PSO_ARMA_sn(y_AR,X_AR,regime_sn_max,regime_sn_sig_max,AR_lags,MA_lags,init_sn[,1],init_sn[,1],nb_iter,nb_part) ## OwnFunction
    dens_chain[1] <- templist[[1]]
    theta <- templist[[2]]
    sigma <- templist[[3]]
    sn <- matrix(1,taille,1)
    sn_sig <- sn
  } else {
    if(dep_MCMC==0) {
      for (i in 1 : nb_init) { #here parfor
        if((i%%20) == 0) {
         #disp(['Iteration ' num2str(i) ' on ' num2str(nb_init)])
        } #self function call
        res_dep <- MS_ARMA_MCMC(y,Expl_X,30,regime,AR_lags_init,MA_lags,CP_or_MS,temper,0,init_sn[,i],init_sn_sig[,i]) ## OwnFunction
        maxi <- max(res_dep$post_dens)
        ind <- which.max(res_dep$post_dens)
        init_reg[i,] <- res_dep$post_regime[ind,]
        sigma_init[,i] <- res_dep$post_sigma[,ind]
        theta_init[,i] <- res_dep$post_theta[,ind]
        init_sn[,i] <- res_dep$post_sn[,ind]
        init_sn_sig[,i] <- res_dep$post_sn_sig[,ind]
        dens[i] <- maxi
      } #
      ind <- which.max(dens)
      regime <- init_reg[ind,1]
      theta <- theta_init[1:(taille_ARMA*regime_sn_max),ind]
      sigma <- sigma_init[1:regime_sn_sig_max,ind]
      sn <- init_sn[,ind]
      sn_sig <- init_sn_sig[,ind]
    } else {
    #disp('Initialisation by Particle Swarm Optimization')
    #[dens_chain[1],theta,sigma] <- PSO_ARMA_sn(y_AR,X_AR,regime_sn_max,regime_sn_sig_max,AR_lags,MA_lags,sn_dep,sn_sig_dep,nb_iter,nb_part) ## OwnFunction
    sn <- sn_dep
    sn_sig <- sn_sig_dep
    } #
  } #
  } #
 regime_sig <- max(sn_sig)


 #############################################################
 ### The starting values are stored in theta (mean parameters, dimension
 ### regime*taille_ARMA x 1), sigma (variances, dimension regime x 1), sn
 ### and sn_sig that are the corresponding state vectors.
 ### n_ii : a matrix (regime x regime) that counts the number of entries in
 ### the state i,j according to the state vector sn
 ### P : Transition Matrix (regime_sn_max x regime_sn_amx)
 ### n_ii_sig and P_sig have the same representation but for the variance.
 #############################################################


 ###############
 ### Initialisation of all the variables for the MCMC. The variables related
 ### to the Hierarchical Dirichlet processes are stored in IHMM and IHMM_sig
 ###############
 P <- matrix(0,regime_sn_max,regime_sn_max)
 alpha <- 0.5*matrix(1,regime_sn_max,regime_sn_max)
 n_ii <- comptage(taille,sn,regime_sn_max)[[1]] ## OwnFunction
 n_ii_sig <- comptage(taille,sn_sig,regime_sn_sig_max)[[1]] ## OwnFunction
 for (q in 1 : regime_sn_max) {
   P[q,] <- dirichlet_sample(alpha[q,] + n_ii[q,]) ## OwnFunction
 } #
 mu_theta <- Draw_mu_theta(theta,Sigma_theta_inv,taille_ARMA,regime,mean_mu_prior_theta,Sigma_mu_prior_theta_inv) ## OwnFunction
 Sigma_theta_inv <- Draw_Sigma_theta(theta,mu_theta,V_prior_inv,dv_prior,regime,taille_ARMA) ## OwnFunction
 Sigma_theta <- solve(Sigma_theta_inv) ## MatlabFunction

 theta_aide <- matrix(0,taille_ARMA*regime_sn_max,1)
 theta_aide[1:(regime_sn_max*taille_ARMA),1] <- theta
 theta <- theta_aide
 if(regime<regime_sn_max) {
   for (q in regime+1 : regime_sn_max) {
     test <- 0
     while(test==0) {
      theta_prop <- as.matrix(mvrnorm(1, mu_theta,Sigma_theta)) ## OwnFunction ## MatlabFunction ## Normal multivari? RESTRICTION A METTRE ?
      test <- 1
      if(AR_lags>0) {
        root_AR <- polyroot(rev(1 - t(theta_prop[2:(1+AR_lags)]))) ## OwnFunction ## MatlabFunction ## Check if complex numbers are allowed
        if(sum(abs(Re(root_AR))>=1)!=0) { ## MatlabFunction
          test <- 0
        } #
      } #
      if(MA_lags>0) {
        root_MA <- polyroot(rev(1 - t(theta_prop[(2+AR_lags):taille_ARMA])))
        if(sum(abs(Re(root_MA))>=1)!=0) {
          test <- 0
        } #
      } #
      if(test == 1) {
        theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)] <- theta_prop
      } #
    } #
    } #
  } #
  sigma_aide <- matrix(0,regime_sn_sig_max,1)
  sigma_aide[1:regime_sn_sig_max] <- sigma
  sigma <- sigma_aide
  nb_sig_tot <- 30
  if(regime_sig<regime_sn_sig_max) {
    for (q in (regime_sig+1) : regime_sn_sig_max ) {
      nb_it <- 0
      sigma[q] <- 1/rgamma(1, e_prior,1/f_prior) ## MatlabFunction #check if distribution is correct
      while(sigma[q]>10) {
        nb_it <- nb_it +1
        if(nb_it>nb_sig_tot) {
          sigma[q] <- runif(1, 0, 1)*10
        } #
      } #
    } #
  } #
  dens <- 0
  dens_chain[1] <- .C("dens_CP_ARMA_c",
                      as.vector(y_AR, mode = "double"),
                      as.vector(X_AR),
                      as.integer(AR_lags),
                      as.integer(MA_lags),
                      as.vector(theta, mode = "double"),
                      as.vector(sigma, mode = "double"),
                      as.vector(sn, mode = "double"),
                      as.vector(sn_sig, mode = "double"),
                      as.integer(1),
                      as.integer(taille),
                      dens = as.double(dens),
                      eps_out = as.vector(eps_t, mode="double"),
                      PACKAGE="groupedtseries")$dens

  burn <- round(0.25*nb_MCMC)

  if(Rolling_MCMC==0){
    #disp('********************************************************************')
    #disp('############ MCMC Starting point : ARMA parameters')
    #disp(['############ ' num2str(theta[1:regime)')])
    #disp(['############ Density of the initial values : ' num2str(dens_chain[1])])
    #disp('############ MCMC Starting point : variance parameters')
    #disp(['############ ' num2str(sigma[1:regime_sig)')])
    #disp(['############ Burn-in period : ' num2str(burn)])
    #disp('********************************************************************')
   } #



   post_theta <- matrix(0,regime_sn_max*taille_ARMA,(nb_MCMC-burn))
   post_diag_n_ii <- matrix(0,regime_sn_max,(nb_MCMC-burn))
   post_sn <- matrix(0,taille,(nb_MCMC-burn)) #these were 8-bit integers in Matlab
   post_sn_sig <- matrix(0,taille,(nb_MCMC-burn)) #these were 8-bit integers in Matlab
   post_sigma <- matrix(0,regime_sn_sig_max,(nb_MCMC-burn))
   post_dens <- matrix(0,(nb_MCMC-burn),1)
   post_mu_theta <- matrix(0,taille_ARMA,(nb_MCMC-burn))
   post_Sigma_theta <- matrix(0,taille_ARMA*taille_ARMA,(nb_MCMC-burn))
   post_Sigma_prior <- matrix(0,2,(nb_MCMC-burn))
   post_alp_kap_eth <- matrix(0,6,(nb_MCMC-burn))
   post_P <- matrix(0,regime_sn_max*regime_sn_max,(nb_MCMC-burn))
   post_P_sig <- matrix(0,regime_sn_sig_max*regime_sn_max,nb_MCMC-burn)
   post_lambda <- matrix(0,regime_sn_sig_max,nb_MCMC-burn)
   post_gamma <- matrix(0,regime_sn_max,(nb_MCMC-burn))
   post_regime <- matrix(0,nb_MCMC-burn,2)
   post_hier_rho <- matrix(0,nb_MCMC-burn,1)
   theta_t <- array(0,dim=c(taille_ARMA,taille,(nb_MCMC-burn)))
   sigma_t <- matrix(0,taille,(nb_MCMC-burn))
   y_for <- matrix(0,nb_forecast,(nb_MCMC-burn-1))
   ################################################
   #### Initialisation IHMM
   ################################################
   IHMM <- list()
   IHMM$pi <- matrix(0,regime+1,1)
   IHMM$alpha <- 1.2
   IHMM$kappa <- 1.8
   IHMM$etha <- 3

   IHMM$prior_eta_a <- Struct_Prior$prior_eta_a
   IHMM$prior_eta_b <- Struct_Prior$prior_eta_b
   IHMM$prior_rho_a <- Struct_Prior$prior_rho_a
   IHMM$prior_rho_b <- Struct_Prior$prior_rho_b
   IHMM$prior_alpha_kappa_a <- Struct_Prior$prior_alpha_kappa_a
   IHMM$prior_alpha_kappa_b <- Struct_Prior$prior_alpha_kappa_b
   IHMM$hier_rho_a <- Struct_Prior$hier_rho_a
   IHMM$hier_rho_b <- Struct_Prior$hier_rho_b

   vect_a <- (IHMM$etha/regime_sn_max)*matrix(1,regime_sn_max,1) #unclear whether matrix multiplication

   for (i in 1 : regime) {
    vect_a[i] <- sum(n_ii[,i])
    } #
   IHMM$gamma <- dirichlet_sample(vect_a)  ## OwnFunction
   IHMM$m <- table_IHMM(IHMM,n_ii,regime_sn_max) ## OwnFunction
   IHMM$r <- overriding(IHMM,regime_sn_max,sn[1]) ## OwnFunction
   IHMM$m_bar <- considered_table(IHMM$m,IHMM$r,regime_sn_max) ## OwnFunction
   IHMM$gamma <- DP_gamma(IHMM$etha,IHMM$m_bar,regime_sn_max) ## OwnFunction
   P <- DP_pi(IHMM,regime_sn_max,n_ii) ## OwnFunction

   templist <- prior(IHMM,regime_sn_max,n_ii) ## OwnFunction
   IHMM$alpha <- templist[[1]]
   IHMM$kappa <- templist[[2]]
   IHMM$etha <- templist[[3]]

   IHMM_sig <- IHMM
   IHMM_sig$gamma <- dirichlet_sample(vect_a) ## OwnFunction
   IHMM_sig$m<- table_IHMM(IHMM_sig,n_ii_sig,regime_sn_sig_max) ## OwnFunction
   IHMM_sig$r <- overriding(IHMM_sig,regime_sn_sig_max,sn_sig[1]) ## OwnFunction
   IHMM_sig$m_bar <- considered_table(IHMM_sig$m,IHMM_sig$r,regime_sn_sig_max) ## OwnFunction
   IHMM_sig$gamma <- DP_gamma(IHMM_sig$etha,IHMM_sig$m_bar,regime_sn_sig_max) ## OwnFunction
   P_sig <- DP_pi(IHMM_sig,regime_sn_sig_max,n_ii_sig) ## OwnFunction
   p_lambda <- matrix(0,regime_sn_sig_max)
   lambda <- matrix(0,regime_sn_sig_max,1)
   if(dependance==1) {
    for (j in 1 : regime_sn_sig_max) {
      lambda[j] <- rbeta(1, IHMM_sig$alpha*IHMM_sig$gamma[j] + IHMM_sig$kappa, IHMM_sig$alpha*(1-IHMM_sig$gamma[j]))
      #not clear whether result shall be a matrix or scalar (in matlab function can give out matrix)
    } #
   } #

   templist <- prior(IHMM_sig,regime_sn_sig_max,n_ii_sig) ## OwnFunction
   IHMM_sig$alpha <- templist[[1]]
   IHMM_sig$kappa <- templist[[2]]
   IHMM_sig$etha <- templist[[3]]
   accept_rho <- 0
   adapt_rho <- 5

 accept_MA <- 0
 count_MA <- 0
 accept_MA_test <- matrix(1,taille_ARMA*regime_sn_max,1)
 count_MA_test <- matrix(1,taille_ARMA*regime_sn_max,1)
 rate_ARMA <- matrix(1,taille_ARMA*regime_sn_max,1)
 adapt_rate_ARMA <- 0.1*matrix(1,taille_ARMA*regime_sn_max,1)
 count_GRAY <- 0
 GRAY_MH <- 0

 #############################
 ### DEBUT DU MCMC
 #############################
 i <- 1
 while(i<=nb_MCMC) {
  #progressbar(i/nb_MCMC)
  if(i%%500==0) {
    #disp('********************************************************************')
    #disp('********************************************************************')
    #disp(['************* TEMPERING : ' num2str(temper) ' ******'])
    #disp(['************* MCMC iteration ' num2str(i) ' on ' num2str(nb_MCMC)])
    #disp('************* Acceptance rate for the ARMA parameters ')
    #disp(rate_ARMA[1:taille_ARMA*regime)') ## OwnFunction
    #disp('************* Adapt rate for theta ')
    #disp(adapt_rate_ARMA[1:taille_ARMA*regime)') ## OwnFunction

    if(CP_or_MS == 'random')  {
    #disp(['************* Acceptance rate for the hierarchical IHMM parameters  : ' num2str(accept_rho/(i-1))])
    #disp(['************* Adapt rate for the hierarchical IHMM parameters  : ' num2str(adapt_rho)])
    } #

    #disp(['************* Number of active regimes for the mean : ' num2str(regime) ' and for the variance : ' num2str(regime_sig)])
    if(MA_lags>0) {
      Taux <- accept_MA/count_MA
      #disp(['************* Average acceptance rate for the ARMA parameters :   ' num2str(Taux)])
      if(sample_all==1) {
      #disp(['************* Acceptance rate for the mean state vector  : ' num2str(GRAY_MH/i)] )
      } else {
      #disp(['************* Acceptance rate for the mean state vector  : ' num2str(GRAY_MH/count_GRAY)])
      } #
    } #
    #disp('************* Number of observations in each regime : ')
    #disp(n_ii[1:regime,1:regime))
    #disp('************* Number of observation in each regime of the variance : ')
    #disp(n_ii_sig[1:regime_sig,1:regime_sig))
    #disp('********************************************************************')
    #disp('********************************************************************')
  } #





  ##############################
  #### MCMC sub-kernel : ARMA parameters
  #### If MA_lags>0, sampling by Metropolis-Hastings using the sMALA
  #### Otherwise, sampling by Gibbs (conjugate distribution)
  #############################
  if(MA_lags!=0) {
    dens <- 0
    dens_stay <- .C("dens_CP_ARMA_c",
                    as.vector(y_AR, mode = "double"),
                    as.vector(X_AR),
                    as.integer(AR_lags),
                    as.integer(MA_lags),
                    as.vector(theta, mode = "double"),
                    as.vector(sigma, mode = "double"),
                    as.vector(sn, mode = "double"),
                    as.vector(sn_sig, mode = "double"),
                    as.integer(1),
                    as.integer(taille),
                    dens = as.double(dens),
                    eps_out = as.vector(eps_t, mode="double"),
                    PACKAGE="groupedtseries")$dens

    inv_Sig_prior <- solve(Sigma_theta)
    det_Sig_prior <- det(Sigma_theta)

    theta_prop <- theta
    saut_MH <- ceiling(runif(1, 0, 1)*regime*taille_ARMA)
    saut_MH_dep <- saut_MH
    dim_ARMA <- regime*taille_ARMA
    for (z in seq(from=1, to=dim_ARMA, by=saut_MH)) {
      if(z+saut_MH>dim_ARMA) {
        fin <- dim_ARMA
        saut_MH <- fin-z+1
      } else {
        fin <- z+saut_MH-1
      } #
      indice <- saut_MH
      templist <- calc_gradient_Hessian(y_AR,X_AR,theta,sigma,sn,sn_sig,regime,AR_lags,MA_lags,inv_Sig_prior,temper) ## OwnFunction ##  Computation of the gradient
      gradient_move <- templist[[1]]
      Hessian_approx_move <- templist[[2]]
      C <- chol(Hessian_approx_move[z:fin,z:fin]) ## MatlabFunction #not solved
      inv_C <- solve(C)
      inv_G_move <- inv_C%*%t(inv_C) #unsure if matrix calculation needed
      ###### New candidate for the ARMA parameters
      mu_prop <- theta[z:fin] + adapt_rate_ARMA[indice]^2/2*inv_G_move%*%gradient_move[z:fin] ## OwnFunction
      theta_prop[z:fin] <- mu_prop + adapt_rate_ARMA[indice]*inv_C*rnorm(saut_MH) ## randn(saut_MH, 1) ok<*MINV> ## OwnFunction

      test <- 1
      deb_etat <- ceiling(z/taille_ARMA)
      fin_etat <- ceiling(fin/taille_ARMA)
      for (q in deb_etat : fin_etat) {
        ind <- (q-1)*taille_ARMA
        if(AR_lags>0) {
          root_AR <- polyroot(rev(1 - t(theta_prop[(ind+2):(ind+1+AR_lags)])))
          if(sum(abs(Re(root_AR))>=1)!=0) {
            test <- 0
          } #
        } #
        if(MA_lags>0) {
          root_MA <- polyroot(rev(1 -t(theta_prop[(ind+2+AR_lags):(ind+taille_ARMA)])))
          if(sum(abs(Re(root_MA))>=1)!=0){
            test <- 0
          } #
        } #
      } #
      count_MA <- count_MA +1
      count_MA_test[indice] <- count_MA_test[indice] +1

      if(test==1) {
        ### Density of the new parameters
        dens <- 0
        dens_move <- .C("dens_CP_ARMA_c",
                        as.vector(y_AR, mode = "double"),
                        as.vector(X_AR),
                        as.integer(AR_lags),
                        as.integer(MA_lags),
                        as.vector(theta_prop, mode = "double"),
                        as.vector(sigma, mode = "double"),
                        as.vector(sn, mode = "double"),
                        as.vector(sn_sig, mode = "double"),
                        as.integer(1),
                        as.integer(taille),
                        dens = as.double(dens),
                        eps_out = as.vector(eps_t, mode="double"),
                        PACKAGE="groupedtseries")$dens

        prior_stay <- 0
        prior_move <- 0
        #### Prior evaluations for the candidate and the current
        #### parameters.
        for (q in deb_etat : fin_etat) {
          prior_stay <- prior_stay -taille_ARMA*0.5*log(2*pi) -0.5*log(det_Sig_prior) -0.5*t(theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)]-mu_theta)%*%inv_Sig_prior%*%(theta[(q-1)*taille_ARMA+1:q*taille_ARMA]-mu_theta)
          prior_move <- prior_move -taille_ARMA*0.5*log(2*pi) -0.5*log(det_Sig_prior) -0.5*t(theta_prop[((q-1)*taille_ARMA+1):(q*taille_ARMA)]-mu_theta)%*%inv_Sig_prior%*%(theta_prop[(q-1)*taille_ARMA+1:q*taille_ARMA]-mu_theta)
        } #

       templist <- calc_gradient_Hessian(y_AR,X_AR,theta_prop,sigma,sn,sn_sig,regime,AR_lags,MA_lags,inv_Sig_prior,temper) ## OwnFunction
        gradient_stay <- templist[[1]]
        Hessian_approx_stay <- templist[[2]]

       C <- chol(Hessian_approx_stay[z:fin,z:fin])
       inv_C <- solve(C)
       inv_G_stay <- inv_C*t(inv_C)
       mu_stay <- as.matrix(theta_prop[z:fin]) + adapt_rate_ARMA[indice]^2/2*inv_G_stay%*%as.matrix(gradient_stay[z:fin])
       inv_sig_stay <- Hessian_approx_stay[z:fin,z:fin]/adapt_rate_ARMA[indice]^2
       det_stay <- det(inv_G_stay*adapt_rate_ARMA[indice]^2)
       prop_stay <- -saut_MH*0.5*log(2*pi) - 0.5*log(det_stay) -0.5*t(as.matrix((theta[z:fin]-mu_stay)))%*%inv_sig_stay%*%as.matrix(theta[z:fin]-mu_stay)


       inv_Sig_move <- Hessian_approx_move[z:fin,z:fin]/adapt_rate_ARMA[indice]^2
       det_move <- det(inv_G_move*adapt_rate_ARMA[indice]^2)
       prop_move <- -saut_MH*0.5*log(2*pi) -0.5*log(det_move) -0.5*t(as.matrix(theta_prop[z:fin]-mu_prop))%*%inv_Sig_move%*%as.matrix(theta_prop[z:fin]-mu_prop)

       if(runif(1, 0, 1)<exp(temper*(dens_move-dens_stay)+prior_move+prop_stay-prior_stay-prop_move)) { ## Metropolis-Hastings ratio
         theta[z:fin] <- theta_prop[z:fin]
       dens_stay <- dens_move
       accept_MA <- accept_MA +1
       accept_MA_test[indice] <- accept_MA_test[indice] +1
       } else {
         theta_prop[z:fin] <- theta[z:fin]
       } #
      } else {
        theta_prop[z:fin] <- theta[z:fin]
      } #
    } #

    ##### Adaptation of the scale (see Atchade and Rosenthal (2001))
    rate_ARMA[saut_MH] <- accept_MA_test[saut_MH]/(count_MA_test[saut_MH])
    adapt_rate_ARMA[saut_MH] <- max(adapt_rate_ARMA[saut_MH] + (rate_ARMA[saut_MH]-0.4)/(i^(0.6)),1e-8)
    if(saut_MH!=saut_MH_dep) {
      rate_ARMA[saut_MH_dep] <- accept_MA_test[saut_MH_dep]/(count_MA_test[saut_MH_dep])
      adapt_rate_ARMA[saut_MH_dep] <- max(adapt_rate_ARMA[saut_MH_dep] + (rate_ARMA[saut_MH_dep]-0.4)/(i^(0.6)),1e-8)
    }

  } else {
    sigma_over_time <- sigma[sn_sig]
    for (q in 1 : regime) {
      ind_etat <- as.logical(sn==q)
      y_current <- y_AR[ind_etat]
      sigma_current <- sigma_over_time[ind_etat]
      X_current <- X_AR[ind_etat,]
      X_current_for_var <- X_current
      X_current_for_mu <- X_current
      for (z in 1 : taille_ARMA) {
        X_current_for_var[,z] = X_current[,z]/sqrt(sigma_current)
        X_current_for_mu[,z] = X_current[,z]/sigma_current
      } #
      Sigma <- solve(temper*(t(X_current_for_var)%*%X_current_for_var) + Sigma_theta_inv)
      mu <- t(t(temper*t(X_current_for_mu)%*%y_current + Sigma_theta_inv%*%mu_theta)%*%Sigma)
      test <- 0
      nb_max_count <- 100
      count <- 1

      while(test==0 && count<nb_max_count) {
        theta_prop <- as.matrix(mvrnorm(1,mu,Sigma)) ## Tirage dans la distribution conditionelle ## MatlabFunction
        root_AR <- polyroot(rev(1 -t(theta_prop[2:taille_ARMA])))
        if(sum(abs(Re(root_AR))>=1)==0) {
          test <- 1
        } #
      } #
      if(test==1) {
        theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)] <- theta_prop
      } #
    } #
  } #

  for (q in (regime+1) : regime_sn_max) {
    theta_prop <- as.matrix(mvrnorm(1,mu_theta,Sigma_theta))
    test <- 1
    if(AR_lags>0) {
      root_AR <- polyroot(rev(1 -t(theta_prop[2:(1+AR_lags)])))
      if(sum(abs(Re(root_AR))>=1)!=0){
        test <- 0
      } #
    } #
    if(MA_lags>0) {
      root_MA <- polyroot(rev(1 - t(theta_prop[(2+AR_lags):taille_ARMA])))
      if(sum(abs(Re(root_MA))>=1)!=0){
        test <- 0
      } #
    } #
    if(test == 1) {
      theta[((q-1)*taille_ARMA+1):(q*taille_ARMA)] <- theta_prop
    } #
  } #


   ###############
   #### Sampling of the hierarchical parameters of the ARMA
   ###############
   mu_theta <- Draw_mu_theta(theta,Sigma_theta_inv,taille_ARMA,regime_sn_max,mean_mu_prior_theta,Sigma_mu_prior_theta_inv) ## OwnFunction
   Sigma_theta_inv <- Draw_Sigma_theta(theta,mu_theta,V_prior_inv,dv_prior,regime_sn_max,taille_ARMA) ## OwnFunction
   Sigma_theta <- solve(Sigma_theta_inv)


   ###############
   #### Sampling of the variance by Gibbs as the conditional distribution is conjugate when the prior is  1/sigma^2 ~
   #### G(e_prior,f_prior)
   ###############
   dens <- 0
   eps_for_sigma <- .C("dens_CP_ARMA_c",
                       as.vector(y_AR, mode = "double"),
                       as.vector(X_AR),
                       as.integer(AR_lags),
                       as.integer(MA_lags),
                       as.vector(theta, mode = "double"),
                       as.vector(sigma, mode = "double"),
                       as.vector(sn, mode = "double"),
                       as.vector(sn_sig, mode = "double"),
                       as.integer(1),
                       as.integer(taille),
                       dens = as.double(dens),
                       eps_out = as.vector(eps_t, mode="double"),
                       PACKAGE="groupedtseries")$eps_out

   e_post <- matrix(0,regime_sn_sig_max,1)
   f_post <- matrix(0,regime_sn_sig_max,1)
   for (q in 1 : regime_sig ) {
     indice <- as.logical(sn_sig==q)
     eps_sig <- eps_for_sigma[indice]
     f_post[q] <- t(eps_sig)%*%eps_sig
   } #
   for (q in 1 : regime_sn_sig_max) {
     count <- sum(n_ii_sig[,q])
     e_post[q] <- 0.5*temper*count + e_prior
     f_post[q] <- 1/(0.5*temper*f_post[q] + 1/f_prior)


     val_sig <- 1/rgamma(1, e_post[q],f_post[q])
     if(val_sig<50) {
      sigma[q] <- val_sig
     } #
  } #

 ###############
 #### Sampling of the hierarchical parameters of the variance
 ###############
 c_post <- regime_sn_sig_max*e_prior + c_0
 d_post <- 1/d_0
 for (q in 1 : regime_sn_sig_max) {
    d_post <- d_post + 1/sigma[q]
 } #
 d_post <- 1/d_post
 f_prop <- rgamma(1, c_post,d_post) ## MatlabFunction
 if (f_prop>0.1) {
  f_prior <- f_prop
 } #

 e_prop <- e_prior + sqrt(MH_e)*rnorm(1) ## MatlabFunction
 if(e_prop>0) {
    prior_move <- -rho_e*e_prop
    prior_stay <- -rho_e*e_prior
    dens_stay <- e_prior*regime_sn_sig_max*log(f_prior) - regime_sn_sig_max*lgamma(e_prior) -(e_prior-1)*sum(log(sigma))
    dens_move <- e_prop*regime_sn_sig_max*log(f_prior) - regime_sn_sig_max*lgamma(e_prop) -(e_prop-1)*sum(log(sigma))

    if(runif(1, 0, 1)<exp(dens_move+prior_move-dens_stay-prior_stay)) {
      e_prior <- e_prop
      accept_e <- accept_e+1
    } #
  } #



  if(force_reg_1==0) {
    ###############
    #### Sampling of the mean state vector using a metropolis-Hastings
    #### algorithm with an approximate model as proposal distribution
    ###############
    if(regime_sn_max>1 || regime_sig>1) {
      if(MA_lags>0) {
        err_out <- double(length=taille_ARMA*taille)
        X_out <- double(length=taille_ARMA*taille)
        #X_AR
        X_AR_out  <- .C("create_X_MA_term_c",
                   as.vector(y_AR, mode = "double"),
                   as.vector(X_AR),
                   as.integer(AR_lags),
                   as.integer(MA_lags),
                   as.vector(theta, mode = "double"),
                   as.vector(sn, mode = "double"),
                   as.integer(taille),
                   err_out = as.vector(err_out, mode="double"),
                   X_out = as.vector(X_out, mode="double"),
                   PACKAGE="groupedtseries")$X_out
       X_AR <- matrix(X_AR_out, nrow=taille, ncol=taille_ARMA)

      } #
      if(MA_lags==0) {
        templist <- launch_FB_Klaassen(y_AR,X_AR,regime_sn_max,theta,sigma,sn,sn_sig,P,AR_lags,MA_lags,1,taille,temper) ## OwnFunction
        sn <- templist[[1]]
      } else {
        if(taille<250) {
          iter_GRAY <- round(5+runif(1, 0, 1)*taille)
        } else {
          iter_GRAY <- round(40 + runif(1, 0, 1)*(150))
        } #

        for (q in seq(from=1, to=taille, by=iter_GRAY)) {

          if(q+iter_GRAY>taille) {
            fin <- taille
          } else {
            fin <- q+iter_GRAY
          } #

          if(regime_sn_max>1) {
            if(choix_Haas==0) {
              templist <- launch_FB_Klaassen(y_AR,X_AR,regime_sn_max,theta,sigma,sn,sn_sig,P,AR_lags,MA_lags,q,fin,temper) ## OwnFunction
              sn <- templist[[1]]
              accept <- templist[[2]]
            } else {
              templist <- Haas_FB_ARMA(y_AR,X_AR,taille,regime_sn_max,theta,sigma,sn,sn_sig,P,AR_lags,MA_lags,q,fin,temper) ## OwnFunction
              sn <- templist[[1]]
              accept <- templist[[2]]
            } #
          } else {
            accept <- 1
          } #

          if((accept>0) && (MA_lags>0)) {
            GRAY_MH <- GRAY_MH+1
          } #
          count_GRAY <- count_GRAY+1
        } #

      } #

    } else {
      count_GRAY <- count_GRAY+1
      GRAY_MH <- GRAY_MH+1
    } #
    templist <- comptage(taille,sn,regime_sn_max) ## OwnFunction
    n_ii <- templist[[1]]
    d_t <- templist[[2]]
    if(sampling_sig==1) {

      sn_sig <- launch_FB_sigma(y_AR,X_AR,taille,regime_sn_sig_max,theta,sigma,sn,sn_sig,P_sig,AR_lags,MA_lags,temper,dependance,d_t,dist_sig,lambda)[[1]] ## OwnFunction
    } else if (sampling_sig==0) {
      sn_sig <- matrix(1,taille,1)
    } #

    templist <- comptage(taille,sn_sig,regime_sn_sig_max) ## OwnFunction
    n_ii_sig <- templist[[1]]
    ###############
    #### Sorting of the variables depending on the number of active
    #### regimes.
    ###############
    templist <- arrange_mean(regime_sn_max,n_ii,taille,sn,theta,taille_ARMA) ## OwnFunction
    theta <- templist[[1]]
    regime <- templist[[2]]
    sn <- templist[[3]]

    templist <- arrange_sig(regime_sn_sig_max,n_ii_sig,taille,sn_sig,sigma) ## OwnFunction
    sigma <- templist[[1]]
    regime_sig <- templist[[2]]
    sn_sig <- templist[[3]]

  } #

  templist <- comptage(taille,sn,regime_sn_max) ## OwnFunction
  n_ii <- templist[[1]]
  d_t <- templist[[2]]
  templist <- comptage(taille,sn_sig,regime_sn_sig_max) ## OwnFunction
  n_ii_sig <- templist[[1]]

  if(dependance==1) {
    templist <- comptage_dep(taille,sn_sig,d_t,regime_sn_sig_max,dist_sig) ## OwnFunction
    n_ii_sig_P <- templist[[1]]
    n_ii_sig_lambda <- templist[[2]]
  } #
  dens <- 0
  dens_chain[i] <- .C("dens_CP_ARMA_c",
                      as.vector(y_AR, mode = "double"),
                      as.vector(X_AR),
                      as.integer(AR_lags),
                      as.integer(MA_lags),
                      as.vector(theta, mode = "double"),
                      as.vector(sigma, mode = "double"),
                      as.vector(sn, mode = "double"),
                      as.vector(sn_sig, mode = "double"),
                      as.integer(1),
                      as.integer(taille),
                      dens = as.double(dens),
                      eps_out = as.vector(eps_t, mode="double"),
                      PACKAGE="groupedtseries")$dens

  ###############
  #### Sampling the variables of the IHMM structure
  ###############
  if(force_reg_1==0) {
    IHMM$m <- table_IHMM(IHMM,n_ii,regime_sn_max) ## OwnFunction
    IHMM$r <- overriding(IHMM,regime_sn_max,sn[1]) ## OwnFunction
    IHMM$m_bar <- considered_table(IHMM$m,IHMM$r,regime_sn_max) ## OwnFunction

    IHMM$gamma <- DP_gamma(IHMM$etha,IHMM$m_bar,regime_sn_max) ## OwnFunction
    P <- DP_pi(IHMM,regime_sn_max,n_ii) ## OwnFunction
    templist <- prior(IHMM,regime_sn_max,n_ii) ## OwnFunction
    IHMM$alpha <- templist[[1]]
    IHMM$kappa <- templist[[2]]
    IHMM$etha <- templist[[3]]

    IHMM_sig$m <- table_IHMM(IHMM_sig,n_ii_sig,regime_sn_sig_max) ## OwnFunction
    IHMM_sig$r <- overriding(IHMM_sig,regime_sn_sig_max,sn_sig[1]) ## OwnFunction
    IHMM_sig$m_bar <- considered_table(IHMM_sig$m,IHMM_sig$r,regime_sn_sig_max) ## OwnFunction

    IHMM_sig$gamma <- DP_gamma(IHMM_sig$etha,IHMM_sig$m_bar,regime_sn_sig_max) ## OwnFunction
    if(dependance==1) {
      P_sig <- DP_pi_sig(IHMM_sig,regime_sn_sig_max,n_ii_sig_P) ## OwnFunction
      for (j in 1 : regime_sn_sig_max) {
        lambda[j] <- rbeta(1, IHMM_sig$alpha*IHMM_sig$gamma[j] + n_ii_sig_lambda[j,1] + IHMM_sig$kappa,IHMM_sig$alpha*(1-IHMM_sig$gamma[j]) + n_ii_sig_lambda[j,2]) ## OwnFunction
      } #
  } else {
    P_sig <- DP_pi_sig(IHMM_sig,regime_sn_sig_max,n_ii_sig) ## OwnFunction
  } #

  #[IHMM_sig$alpha,IHMM_sig$kappa,IHMM_sig$etha] <- prior(IHMM_sig,regime_sn_sig_max,n_ii_sig) ## OwnFunction

  if('random' == CP_or_MS) {
    prop_rho_a <- IHMM$prior_rho_a + rnorm(1)*adapt_rho
    if(prop_rho_a>0) {
      prior_move <- log_gampdf(prop_rho_a,IHMM$hier_rho_a,IHMM$hier_rho_b) ## OwnFunction
      prior_stay <- log_gampdf(IHMM$prior_rho_a,IHMM$hier_rho_a,IHMM$hier_rho_b) ## OwnFunction
      alpha_move <- cbind(prop_rho_a, IHMM$prior_rho_b)
      alpha_stay <- cbind(IHMM$prior_rho_a, IHMM$prior_rho_b)
      rho_sig <- IHMM_sig$kappa/(IHMM_sig$alpha+IHMM_sig$kappa)
      rho <- IHMM$kappa/(IHMM$alpha+IHMM$kappa)
      log_dens_move <- ln_dirichlet_pdf(cbind(rho,1-rho),alpha_move) + ln_dirichlet_pdf(cbind(rho_sig,1-rho_sig),alpha_move) ## OwnFunction
      log_dens_stay <- ln_dirichlet_pdf(cbind(rho,1-rho),alpha_stay) + ln_dirichlet_pdf(cbind(rho_sig,1-rho_sig),alpha_stay) ## OwnFunction

      if(exp(log_dens_move+prior_move-log_dens_stay-prior_stay)>runif(1, 0, 1)) {
        accept_rho <- accept_rho +1
        IHMM$prior_rho_a <- prop_rho_a
        IHMM_sig$prior_rho_a <- prop_rho_a
      } #
    } #
    adapt_rho <- max(adapt_rho + (accept_rho/i-0.4)/(i^(0.6)),1e-8)
  } #
  } #




  ###############
  #### After the burn-in period, the realizations are stored.
  ###############
  if(i>burn) {

    post_diag_n_ii[1:regime_sn_max,i-burn] <- t(diag(n_ii))
    post_regime[i-burn,1] <- regime
    post_dens[i-burn] <- dens_chain[i]
    post_theta[1:(taille_ARMA*regime_sn_max),i-burn] <- theta
    post_sigma[1:regime_sn_sig_max,i-burn] <- sigma
    post_sn[,i-burn] <- sn
    post_mu_theta[,i-burn] <- mu_theta
    post_Sigma_theta[,i-burn] <- as.vector(Sigma_theta) ## MatlabFunction
    post_Sigma_prior[2,i-burn] <- f_prior
    post_Sigma_prior[1,i-burn] <- e_prior
    post_alp_kap_eth[1,i-burn] <- IHMM$alpha
    post_alp_kap_eth[2,i-burn] <- IHMM$kappa
    post_alp_kap_eth[3,i-burn] <- IHMM$etha

    post_alp_kap_eth[4,i-burn] <- IHMM_sig$alpha
    post_alp_kap_eth[5,i-burn] <- IHMM_sig$kappa
    post_alp_kap_eth[6,i-burn] <- IHMM_sig$etha
    post_hier_rho[i-burn] <- IHMM$prior_rho_a

    post_gamma[1:regime_sn_max,i-burn] <- IHMM$gamma
    post_P[1:(regime_sn_max*regime_sn_max),i-burn] <- as.vector(P[1:regime_sn_max,1:regime_sn_max])
    post_P_sig[1:(regime_sn_sig_max*regime_sn_sig_max),i-burn] <- as.vector(P_sig[1:regime_sn_sig_max,1:regime_sn_sig_max])
    post_lambda[,i-burn] <- lambda



    for (t in 1 : taille) {
      etat <- sn[t]

      theta_t[1,t,i-burn] <- theta[(etat-1)*taille_ARMA+1]/(1-sum(theta[((etat-1)*taille_ARMA+2):((etat-1)*taille_ARMA+1+AR_lags)]))
      theta_t[-1,t,i-burn] <- theta[((etat-1)*taille_ARMA+2):(etat*taille_ARMA)]
    } #
    sigma_t[,i-burn] <- sigma[sn_sig]

    post_regime[i-burn,2] <- regime_sig
    post_sn_sig[,i-burn] <- sn_sig
    ###########################
    #### PREDICTIVE LIKELIHOOD : Forecasting of the time series
    ###########################
    if(nb_forecast>0) {
      etat <- sn[length(sn)]
      etat_sig <- sn_sig[length(sn_sig)]
      d_t_for <- d_t[length(d_t)]
      sigma_new <- sigma
      beta_new <- theta
      P_new <- P


      if(dependance==1) {
        for (j in 1 : regime_sn_sig_max) {
          prob <- (1-lambda[j])*matrix(1,1,regime_sn_sig_max)/(regime_sn_sig_max-1)
          prob[j] <- lambda[j]
          p_lambda[j,] <- prob
        } #
      } #


     if(MA_lags>0) {
       dens <- 0
      eps_for_MH <- .C("dens_CP_ARMA_c",
                       as.vector(y_AR, mode = "double"),
                       as.vector(X_AR),
                       as.integer(AR_lags),
                       as.integer(MA_lags),
                       as.vector(theta, mode = "double"),
                       as.vector(sigma, mode = "double"),
                       as.vector(sn, mode = "double"),
                       as.vector(sn_sig, mode = "double"),
                       as.integer(1),
                       as.integer(taille),
                       dens = as.double(dens),
                       eps_out = as.vector(eps_t, mode="double"),
                       PACKAGE="groupedtseries")$eps_out

      X_AR[-1,AR_lags+2] <- eps_for_MH[1:(length(eps_for_MH)-1)]
      eps <- eps_for_MH[length(eps_for_MH)]
     } #



    X_aide <- X_AR[length(X_AR[,1]),]
    if(AR_lags>0) {
      if(AR_lags>1) {
        for (r in seq(from=AR_lags, to=2, by=-1)) {
          X_aide[1+r] <- X_aide[r]
        } #
      } #
      X_aide[2] <- y_AR[length(y_AR)]
    } #

    if(MA_lags>0) {
      if(MA_lags>1) {
        for (r in seq(from=MA_lags, to=2, by=-1)) {
          X_aide[1+AR_lags+r] <- X_aide[AR_lags+r]
        } #
      } #
      X_aide[1+AR_lags+1] <- eps
    } #

    for (q in 1 : nb_forecast) {
      if(force_reg_1==0) {
        etat_prev<- etat
        etat <- multinomialrnd(P_new[etat,]) #OwnFunction
        if(sampling_sig==1) {
          d_t_for <- (d_t_for+1)*(etat_prev==etat) + 1*(etat_prev!=etat)
          if(dependance==1 && d_t_for<dist_sig) {
            etat_sig <- multinomialrnd(p_lambda[etat_sig,]) ## OwnFunction
          } else {
            etat_sig <- multinomialrnd(P_sig[etat_sig,]) ## OwnFunction
          } #
        } #

      } else {
        etat <- 1
        etat_sig <- etat
      } #

      if(q==1) {
        eps <- rnorm(1)*sqrt(sigma_new[etat_sig])
        y_for[q,i-burn] <- t(as.matrix(X_aide))%*%as.matrix(beta_new[((etat-1)*taille_ARMA+1):(etat*taille_ARMA)]) + eps
      } else {
        if(AR_lags>0) {
          if(AR_lags>1) {
            for (r in seq(from=AR_lags, to=2, by=-1)) {
              X_aide[1+r] <- X_aide[r]
            } #
          } #
          X_aide[2] <- y_for[q-1,i-burn]
        } #
        if(MA_lags>0) {
          if(MA_lags>1) {
            for (r in seq(from=MA_lags, to=2, by=-1)) {
              X_aide[1+AR_lags+r] <- X_aide[AR_lags+r]
            } #
          } #
          X_aide[1+AR_lags+1] <- eps
        } #
        eps <- rnorm(1)*sqrt(sigma_new[etat_sig])
        y_for[q,i-burn] <- t(as.matrix(X_aide))%*%as.matrix(beta_new[((etat-1)*taille_ARMA+1):(etat*taille_ARMA)]) + eps
      } #
    } #

  } #

  } #

 i <- i+1
} #

 ###############
 #### End of the MCMC. Everything is gathered in the result structure.
 ###############
 result <- list()
 result$post_sn <- post_sn
 result$post_sn_sig <- post_sn_sig
 result$post_theta <- post_theta
 result$post_sigma <- post_sigma
 result$post_dens <- post_dens
 result$dens_chain <- dens_chain
 result$post_mu_theta <- post_mu_theta
 result$post_Sigma_theta <- post_Sigma_theta
 result$post_sigma_prior <- post_Sigma_prior
 result$post_alp_kap_eth <- post_alp_kap_eth
 result$post_gamma <- post_gamma
 result$post_regime <- post_regime
 result$post_P <- post_P
 result$post_P_sig <- post_P_sig
 result$post_lambda <- post_lambda
 result$post_hier_rho <- post_hier_rho
 result$theta_t <- theta_t
 result$sigma_t <- sigma_t
 result$mean_theta <- apply(theta_t, 1:2, mean)
 result$y_for <- y_for
 result$post_diag_n_ii <- post_diag_n_ii
 dimension <- dim(theta_t)
 result$std_theta <- matrix(0,taille_ARMA,taille)

 for (j in 1 : dimension[2]) {
   for (i in 1 : dimension[3]) {
     result$std_theta[,j] <- result$std_theta[,j] + sqrt((theta_t[,j,i]-result$mean_theta[,j])*(theta_t[,j,i]-result$mean_theta[,j]))
   } #
   result$std_theta[,j] <- result$std_theta[,j]/dimension[3]
 } #
 result$mean_sigma <- mean(sigma_t)
 result$conf70_sigma <- matrix(0,3,taille)
 result$conf70_sigma[1,] <- quantile(t(result$sigma_t),0.15)
 result$conf70_sigma[2,] <- quantile(t(result$sigma_t),0.5)
 result$conf70_sigma[3,] <- quantile(t(result$sigma_t),0.85)

  result$conf70_theta <- array(0,dim=c(3,taille,taille_ARMA))
  a <- aperm(result$theta_t,c(3, 2, 1))
  for (q in 1 : taille_ARMA) {
    result$conf70_theta[1,,q] <- quantile(a[,,q], 0.15)
    result$conf70_theta[2,,q] <- quantile(a[,,q], 0.5)
    result$conf70_theta[3,,q] <- quantile(a[,,q], 0.85)
  } #
  result$std_sigma <- t(apply(sigma_t,2,sd))

  result$AR <- AR_lags
  result$MA <- MA_lags
  result$regime <- regime
  result$regime_sig <- regime_sig
  sn_mode <- calc_mode_sn(post_sn) #OwnFunction
  result$sn_mode <- sn_mode

  k <- max(max(post_regime))
  max_reg <- max(max(post_regime))
  dimension <- dim(t(result$post_sn))
  N <- dimension[1]
  forward <- matrix(0,taille,k)
  test_CP <- matrix(0,taille,1)
  prob_cass<- matrix(0,taille,1)
  prob_reg <- matrix(0,max_reg,2)
  for (z in 1 : N) {
    regime_current <- post_regime[z,1]
    regime_sig_current <- post_regime[z,2]
    struct_break <- matrix(0,regime_current,3)
    prob_reg[regime_current,1] <- prob_reg[regime_current,1]+1
    prob_reg[regime_sig_current,2] <- prob_reg[regime_sig_current,2]+1
    for (q in 1 : taille) {
      etat <- result$post_sn[q,z]
      forward[q,etat] <- forward[q,etat] + 1

      if(struct_break[etat,1]==0) {
        struct_break[etat,1] <- 1
        struct_break[etat,2] <- q
      } else if (struct_break[etat,2]==q-1) {
        struct_break[etat,1] <- struct_break[etat,1] +1
        struct_break[etat,2] <- q
      } #
      struct_break[etat,3] <- struct_break[etat,3] +1
      if (q>1) {
        if(etat!=etat_prev) {
          prob_cass[q] <- prob_cass[q] +1
        } #
      } #
      etat_prev <- etat

    } #
    isCP <- matrix(0,regime_current,1)
    for (q in 1 : regime_current) {
      if(struct_break[q,3]==struct_break[q,1]) {
        isCP[q] <- 1
      } else {
        isCP[q] <- 0
      } #
    } #
    for (q in 1 : taille) {
      if(isCP[result$post_sn[q,z]]==1) {
        test_CP[q] <- test_CP[q] +1
      } #
    } #
  } #
  forward <- forward/(N)
  test_CP <- test_CP/(N)
  prob_reg <- prob_reg/(N)
  result$prob_sn <- forward
  result$isCP <- test_CP
  result$prob_cass <- prob_cass/(N)
  result$prob_reg <- prob_reg
  if (MA_lags>0) {
    result$accept_MA <- accept_MA/count_MA
    result$accept_sn <- GRAY_MH/count_GRAY
  } #


  if(display_graph==1) {
    #subplot(3,3,1),plot(y_AR),title('Time series')
    #subplot(3,3,2),plot(result$prob_cass),set(gca,'YLim',[0 1]),title('Prob. break')
    #subplot(3,3,3),plot(result$isCP),set(gca,'YLim',[0 1]),title('Prob. Change-Point')

    # str <- 'Intercept '
    # for (  i in 1 : AR_lags ) {
    # str <- [str 'AR lags  ' num2str(i)] ##ok<AGROW>
    # } #
    # for (  i in 1 : MA_lags ) {
    # str <- [str 'MA lags  ' num2str(i)] ##ok<AGROW>
    # } #

    # subplot(3,3,4),plot(result$mean_theta(1,]),set(gca,'YLim',[-3 5]),hold on,plot(result$mean_theta(1,]+1.96*result$std_theta(1,],'--r'),plot(result$mean_theta(1,]-1.96*result$std_theta(1,],'--r'),legend(str(1,]),title(' Intercept ')
    # if(AR_lags>0) {
    # subplot(3,3,5),plot(t(result$mean_theta(2:AR_lags+1,])),set(gca,'YLim',[-1 1]),hold on,plot(t(result$mean_theta(2:AR_lags+1,])+1.96*t(result$std_theta(2:AR_lags+1,]),'--r'),plot(t(result$mean_theta(2:AR_lags+1,])-1.96*t(result$std_theta(2:AR_lags+1,]),'--r'),legend(str(2:AR_lags+1,]),title(' AR lags ')
    # } #
    # if(MA_lags>0) {
    # subplot(3,3,6),plot(t(result$mean_theta(AR_lags+2:taille_ARMA,])),set(gca,'YLim',[-1 1]),hold on,plot(t(result$mean_theta(AR_lags+2:taille_ARMA,])+1.96*t(result$std_theta(AR_lags+2:taille_ARMA,]),'--r'),plot(t(result$mean_theta(AR_lags+2:taille_ARMA,])-1.96*t(result$std_theta(AR_lags+2:taille_ARMA,]),'--r'),legend(str(AR_lags+2:taille_ARMA,]),title(' MA lags ')
    # } #


    # subplot(3,3,7),bar(result$prob_reg),title('Prob. Nb regime(s)')
    # subplot(3,3,8),plot(result$prob_sn),set(gca,'YLim',[0 1]),title('Prob. state')
    # subplot(3,3,9),plot(result$dens_chain),title('Density')
  } #
  if(Rolling_MCMC==0) {
    nb_reg <- apply(result$post_regime, 1, function(x) as.numeric(names(which.max(table(x)))))
    N <- max(dim(post_regime))
    sn_mode <- matrix(0,taille,2)
    max_dens <- 0
    for (i in 1 : N) {
      if(post_regime[i,1]==nb_reg[1] && post_regime[i,2]==nb_reg[2]) {
        if(max_dens==0) {
          sn_mode[,1] <- post_sn[,i]
          sn_mode[,2] <- post_sn_sig[,i]
          max_dens <- post_dens[i]
          theta_at_mode <- result$post_theta[1:(nb_reg[1]*taille_ARMA),i]
          sigma_at_mode <- result$post_sigma[1:nb_reg[2],i]
        } else {
          if(post_dens[i]>max_dens) {
            sn_mode[,1] <- post_sn[,i]
            sn_mode[,2] <- post_sn_sig[,i]
            max_dens <- post_dens[i]
            theta_at_mode <- result$post_theta[1:(nb_reg[1]*taille_ARMA),i]
            sigma_at_mode <- result$post_sigma[1:nb_reg[2],i]
          } #

        } #
      } #
    } #
    result$sn_at_mode <- sn_mode
    result$dens_at_mode <- max_dens
    result$mode_regime <- nb_reg
    result$theta_at_mode <- theta_at_mode
    result$sigma_at_mode <- sigma_at_mode

    err_out <- double(length=taille_ARMA*taille)
    X_out <- double(length=taille_ARMA*taille)

    eps_for_MH  <- .C("create_X_MA_term_c",
                    as.vector(y_AR, mode = "double"),
                    as.vector(X_AR),
                    as.integer(AR_lags),
                    as.integer(MA_lags),
                    as.vector(theta_at_mode, mode = "double"),
                    as.vector(sn_mode[,1], mode = "double"),
                    as.integer(taille),
                    err_out = as.vector(err_out, mode="double"),
                    X_out = as.vector(X_out, mode="double"),
                    PACKAGE="groupedtseries")$err_out

    result$eps_at_mode <- eps_for_MH
  } #

  result$choix_Haas <- choix_Haas
  result$dependance <- dependance
  result$dist_sig <- dist_sig


  if(nb_MCMC>50) {
    #result <- label_bivar_kmean(result,display_graph) ## OwnFunction
  } #

}
                     #################################
                     ######### Analyse au mode - A voir ensemble si on fait cette analyse.
                     #################################
                     # [~, nb_reg] <- max(result$prob_reg)
                     # indice <- 0
                     # for (  z in 1 : N ) {
                     #    if(nb_reg==result$post_regime[z))
                       #        if(indice==0)
                       #            max_dens <- result$post_dens[z)
                       #            indice <- z
                       #        else if(max_dens<result$post_dens[z))
                       #            max_dens <- result$post_dens[z)
                       #            indice <- z
                       #           } #
                       #       } #
                       #    } #
                       #
                       # mode_theta <- result$post_theta[1:nb_reg*taille_ARMA,indice]
                       # mode_sigma <- result$post_sigma[1:nb_reg,indice]
                       # mode_P <- reshape(result$post_P[1:nb_reg^2,indice],nb_reg,nb_reg)
                       # for (  q in 1 : nb_reg ) {
                       #     mode_P(q,] <- mode_P(q,]/sum(mode_P(q,])
                       #    } #
                       #
                       # result$mode_theta <- mode_theta
                       # result$mode_sigma <- mode_sigma
                       # result$mode_P <- mode_P
                       # result$mode_sn <- result$post_sn[,indice]
                       # sn <- result$mode_sn
                       # iter <- 500
                       # nb_save <- floor(iter/5)
                       # mode_sn <- int8(matrix(0,taille,nb_save))
                       # iter_burn <- 0
                       # ['Debut Backward ']
                       # GRAY_MH <- 0
                       # count_GRAY <- 0
                       # tic
                       # for (  i in 1 : iter ) {
                       #     if(taille<250)
                       #         iter_GRAY <- taille
                       #     else
                       #         iter_GRAY <- round(50 + runif(1, 0, 1)*(taille-50))
                       #        } #
                       #
                       #     for (  q in 1:iter_GRAY : taille ) {
                       #         if(q+iter_GRAY>taille)
                       #             fin <- taille
                       #         else
                       #             fin <- q+iter_GRAY
                       #            } #
                       #         if(MA_lags>0)
                       #             [X_AR[,AR_lags+2:taille_ARMA)] <- create_X_MA_term_sn(y_AR,X_AR,taille,AR_lags,MA_lags,mode_theta,sn)
                       #            } #
                       #         [sn accept] <- GRAY_ARMA_part_all_ss_beam(y_AR,X_AR,taille,nb_reg,mode_theta,mode_sigma,sn,mode_P,AR_lags,MA_lags,q,fin)
                       #         if(accept>0 && MA_lags>0)
                       #             GRAY_MH <- GRAY_MH+1
                       #            } #
                       #         count_GRAY <- count_GRAY+1
                       #        } #
                       #     if(mod(i,5)==0)
                       #         iter_burn <- iter_burn +1
                       #         mode_sn[,iter_burn) <- sn
                       #        } #
                       #    } #
                       # result$post_back_sn <- mode_sn
                       # result$MH_mode <- GRAY_MH/count_GRAY
                       # prob_cass_mode <- matrix(0,taille,1)
                       # prob_regime_mode <- matrix(0,taille,nb_reg)
                       #
                       # for (  i in 1 : iter_burn ) {
                       #     etat <- mode_sn(1,i)
                       #     prob_regime_mode(1,etat) <- prob_regime_mode(1,etat)+1
                       #     for (  t in 2 : taille ) {
                       #         etat <- mode_sn(t,i)
                       #         prob_regime_mode(t,etat) <- prob_regime_mode(t,etat)+1
                       #         if(mode_sn(t-1,i)!=etat)
                       #             prob_cass_mode(t) <- prob_cass_mode(t) +1
                       #            } #
                       #        } #
                       #    } #
                       # toc
                       # prob_cass_mode <- prob_cass_mode/iter_burn
                       # result$prob_cass_mode <- prob_cass_mode
                       # prob_regime_mode <- prob_regime_mode/iter_burn
                       # result$prob_regime_mode <- prob_regime_mode
