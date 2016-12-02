/*********************************************************************
 * Demo.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab.  This example shows how to use 2D matricies.  This may
 *
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 * For more information, see my site: www.shawnlankton.com
 * by: Shawn Lankton
 *
 ********************************************************************/
#include <math.h>
#include <float.h>


/*
  La fonction calcule le Forward-Backward d'un MSM model avec 3 ?tats differents dans la distribution M.
 *Pour l'appeler :
 * [sn_move,log_like,log_q_prior,forward_out,X_MA] = FB_Klaassen(y,regime,sigma_sorted,p_mat,X,theta_sorted,sn_current,(*AR_lags),MA_lags,deb,fin,path_dep,rand(fin-deb+1,1));
 *
 */
void FB_Klaassen(double *y, int *regime_sig, double *sigma, double *p_sig, double *X, double *theta,
                 double *sn_current, int *AR_lags, int *MA_lags, int *deb, int *fin, int *path_dep,
                 double *Unif_rand, double *temper, int *taille, double *sn_move, double *log_like_out,
                 double *log_q_prior_move, double *forward, double *epsilon_Haas)
{

//declare variables
    double y_sq,maxi,Q,sig_current,log_like,transition,unif;
    double error_tilde,eps_reg,q_move,q_stay,prior_move,prior_stay,val_prev;
    int taille_current;
    int t,taille_eps,incr,i,z,q,j,taille_ARMA,index;

    double PI = 3.141592;

    (*deb) = (*deb)-1;
    (*fin) = (*fin)-1;
    val_prev = 0;

    /*
     Creation des variables de sortie. On n'alloue la place sans en sp?cifier la valeur pour gagner du temps
    */


//     mexPrintf("Creation sn");
    Q = 0;
    for(t=0;t<(*taille);t++)
    {
        sn_move[t] = sn_current[t];
//         mexPrintf("   sn_move[%d] = %f   sn_stay = %f   \n",t,sn_move[t],sn_current[t]);
    }
//     mexPrintf("Creation prob init");

    double prob_init[(*regime_sig)];


    for(i=0;i<(*regime_sig);i++)
    {
        prob_init[i] = 0;
    }
    if((*deb)<=0)
    {
        prob_init[0] = 1;
    }
    else
    {
        prob_init[(int) sn_current[(*deb)-1]] = 1;
    }
    taille_current = (*fin)-(*deb)+1;

    taille_eps = (*fin)+1;

//     mexPrintf("Debut Haas");
    taille_ARMA = 1+(*AR_lags)+(*MA_lags);

    for(t=0;t<(*deb);t++)
    {
        //mexPrintf("Iteration %d",t);
        i = (int) sn_current[t];
        eps_reg = 0;
        index = i*taille_ARMA;
        for(z=0;z<(*AR_lags)+1;z++) {
            eps_reg = eps_reg + X[z*(*taille) + t]*theta[index+z];
        }
        if((*MA_lags)>0) {
            eps_reg = eps_reg + val_prev*theta[index+(*AR_lags)+1+z];
            for(z=(*MA_lags)-1;z>0;z--) {
                if(t-z-1>=0) {
                    i = (int) sn_current[t-z-1];
                    eps_reg = eps_reg + epsilon_Haas[t-z-1 + (taille_eps)*i]*theta[index+(*AR_lags)+1+z];
                }
            }
        }
        val_prev = eps_reg;
        epsilon_Haas[t + taille_eps*i] = y[t]-eps_reg;

    }


//     mexPrintf("Debut Forward");
    log_like = 0;
    transition = 0;
    for(t=0;t<taille_current;t++)
    {
       //mexPrintf("Iteration %d",t);
       maxi = -1e8;
       sig_current = sigma[t+(*deb)];
       for(i=0;i<(*regime_sig);i++)
       {


           if(t>0)
           {
               transition = 0;
               for(z=0;z<(*regime_sig);z++)
               {
                    transition = transition + p_sig[z + i*(*regime_sig)]*forward[(t-1) + z*taille_current];
               }
               if((*MA_lags)>0)
               {
                   eps_reg = 0;
                   if(transition!=0)
                   {
                       for(z=0;z<(*regime_sig);z++) {
                           eps_reg = eps_reg + epsilon_Haas[((*deb)+t-1) + taille_eps*z]*p_sig[z + i*(*regime_sig)]*forward[(t-1) + z*taille_current]/transition;
                       }
                   }
                   error_tilde = eps_reg;
               }
           }
           else
           {
               transition = p_sig[0 + i*(*regime_sig)];
               if((*deb)>0)
               {
                   z = (int) sn_current[(*deb)-1];
                   error_tilde = epsilon_Haas[((*deb)+t-1) + taille_eps*z];
               }
               else
               {
                   error_tilde = 0;
               }
           }


           eps_reg = 0;
           index = i*taille_ARMA;
           for(z=0;z<(*AR_lags)+1;z++)
           {
               eps_reg = eps_reg + X[z*(*taille) + ((*deb)+t)]*theta[index+z];
           }
           if((*MA_lags)>0)
           {
               eps_reg = eps_reg + error_tilde*theta[index+(*AR_lags)+1];
               for(z=1;z<(*MA_lags);z++)
               {
                   if(t+(*deb)-z-1>=0)
                   {
                       eps_reg = eps_reg + epsilon_Haas[((*deb)+t)-z-1 + taille_eps*i]*theta[index+(*AR_lags)+1+z];
                   }
               }
           }
//            mexPrintf("  error_tilde[%d,%d] = %f   ",t+(*deb),i,error_tilde);
//            mexPrintf("  transition[%d,%d] = %f  \n",t+(*deb),i,transition);
           epsilon_Haas[((*deb)+t) + taille_eps*i] = y[(*deb)+t]-eps_reg;
           if(transition!=0)
           {

                y_sq = epsilon_Haas[((*deb)+t) + taille_eps*i]*epsilon_Haas[((*deb)+t) + taille_eps*i];
                forward[t + taille_current*i] = (*temper)*(-0.5*log(2*PI*sig_current)-0.5*y_sq/sig_current) + log(transition);
                if(maxi<forward[t + taille_current*i])
                {
                   maxi = forward[t + taille_current*i];
                }
           }
           else
           {
               forward[t + taille_current*i] = 0;
           }

//           mexPrintf("  forward[%d,%d] = %f",t + taille_current*i,i,forward[t + taille_current*i]);
       }

       Q = 0;
       for(i=0;i<(*regime_sig);i++)
       {
           if(forward[t + taille_current*i]!=0)
           {
               forward[t + taille_current*i] = exp(forward[t + taille_current*i]-maxi);
               Q = Q + forward[t + taille_current*i];
           }
       }

       for(i=0;i<(*regime_sig);i++)
       {
           forward[t + taille_current*i] = forward[t + taille_current*i]/Q;

       }
//        mexPrintf("\n");
//        for(i=0;i<(*regime_sig);i++)
//        {
//            mexPrintf("  Epsilon[%d,%d] = %f",t + taille_current*i,i,epsilon_Haas[((*deb)+t) + (*fin)*i]);
//        }
//        mexPrintf("\n");
       log_like = log_like + log(Q) + maxi;

//        mexPrintf("log_dens = %f\n",log_like);
//         mexPrintf("maxi = %f\n",Q);
//             mexPrintf("theta[%d] = %f\n",1,theta[1]);
//             mexPrintf("theta[%d] = %f\n",2,theta[2]);
//             mexPrintf("theta[%d] = %f\n",3,theta[3]);
    }

    /*
      Procedure Backward
    */
    /*
     Creation d'une variable. On n'alloue la place sans en sp?cifier la valeur pour gagner du temps
    */

    double back_trans[(*regime_sig)];

//     mexPrintf("Debut Back");
    incr = 0;
    q_move = 0;
    q_stay = 0;
    prior_move = 0;
    prior_stay = 0;
    if((*fin)>=((*taille)-1))
    {

        unif = Unif_rand[taille_current-1-incr];
        Q = forward[taille_current-1];
        i = 0;
        while(Q<unif) {
            i = i+1;
            Q = Q + forward[taille_current-1 + taille_current*i];
        }
        sn_move[(*taille)-1] = i;
        q_move = log(forward[taille_current-1 + taille_current*i]);
        i = sn_current[(*taille)-1];
        q_stay = log(forward[taille_current-1 + taille_current*i]);
        incr = 1;
        (*fin) = (*taille)-2;

//         mexPrintf("fin = %d\n",fin);
//         mexPrintf("taille = %d\n",taille);
    }
//     for(i=0;i<(*regime_sig);i++)
//     {
//         for(j=0;j<(*regime_sig);j++)
//         {
//             mexPrintf("p[%d,%d] = %f\n",i,j,p_sig[i + j*(*regime_sig)]);
//         }
//     }

    //     mexPrintf("indice = %d\n",i);

    for(t=(*fin);t>=(*deb);t--) {
        maxi = -1e9;
        q = (int) sn_move[t+1];
//         mexPrintf("   sn_move[%d] = %f   sn_stay = %f   ",t+1,sn_move[t+1],sn_current[t+1]);
        Q = 0;
        for(z=0;z<(*regime_sig);z++)
        {
            back_trans[z] = forward[taille_current-1-incr + z*taille_current]*p_sig[z + q*(*regime_sig)];
            Q = Q + back_trans[z];

        }
//          mexPrintf("\n");
        for(i=0;i<(*regime_sig);i++)
        {
            back_trans[i] = back_trans[i]/Q;
//             mexPrintf("backward[%d,%d] = %f",taille_current-1-incr,i,back_trans[i]);
        }
//         mexPrintf("\n");
        unif = Unif_rand[taille_current-1-incr];
        Q = back_trans[0];
        i = 0;

        while(Q<unif && (i+1<(*regime_sig))) {
            i = i+1;
            Q = Q + back_trans[i];
        }

        sn_move[t] = i;

        if((*path_dep)==1)
        {
            q_move = q_move + log(back_trans[i]);
            prior_move = prior_move + log(p_sig[i + q*(*regime_sig)]);
            j = (int) sn_current[t+1];

            if(q==j)
            {
              i = (int) sn_current[t];
              q_stay = q_stay + log(back_trans[i]);
              prior_stay = prior_stay + log(p_sig[i + q*(*regime_sig)]);
            }
            else
            {
                Q = 0;
                for(z=0;z<(*regime_sig);z++) {

                    back_trans[z] = forward[taille_current-1-incr + z*taille_current]*p_sig[z + j*(*regime_sig)];
                    Q = Q + back_trans[z];
//                     mexPrintf("  forward[%d,%d] = %f",taille_current-1-incr,z,forward[taille_current-1-incr + z*taille_current]);
                }
//                 mexPrintf("\n");
                for(i=0;i<(*regime_sig);i++) {
                    back_trans[i] = back_trans[i]/Q;
//                     mexPrintf("  backward[%d,%d] = %f",taille_current-1-incr,i,back_trans[i]);
                }
//                  mexPrintf("\n");
                i = (int) sn_current[t];
                q_stay = q_stay + log(back_trans[i]);
                prior_stay = prior_stay + log(p_sig[i + j*(*regime_sig)]);
            }
//              mexPrintf("  q move[%d] = %f",taille_current-1-incr,q_move);
//              mexPrintf("  q stay[%d] = %f",taille_current-1-incr,q_stay);
//              mexPrintf("\n");
        }

        incr = incr +1;


    }


    if((*path_dep)==1)
    {
        if((*deb)!=0)
        {
           i = (int) sn_move[(*deb)-1];
           q = (int) sn_move[(*deb)];
           prior_move = prior_move + log(p_sig[i + q*(*regime_sig)]);
           i = (int) sn_current[(*deb)-1];
           j = (int) sn_current[(*deb)];
           prior_stay = prior_stay + log(p_sig[i + j*(*regime_sig)]);
        }

        log_q_prior_move[0] = q_move;
        log_q_prior_move[1] = q_stay;
        log_q_prior_move[2] = prior_move;
        log_q_prior_move[3] = prior_stay;
    }



    return;
}
