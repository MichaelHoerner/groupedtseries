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
 *[sn_sig log_like forward] = FB_sigma_c(eps_t,regime_sig,sigma,p_sig,Unif_rand)
 */
void FB_sigma(double *eps_t, int *regime_sig, double *sigma, double *p_sig, double *Unif_rand,
              double *temper, int *dependance, double *d_t, double *dist, double *p_lambda,
              int *taille, int *number_ts_segments, int *seg_start_index,
              int *seg_end_index, int *seg_length,
              double *sn_move, double *log_like_out, double *forward)
{


//declare variables
    double eps_sq ,maxi,Q,sig_current,log_like,transition,unif;
    double eps_reg;
    int  taille_current;
    int t,incr,i,z,q,j,taille_ARMA,index, seg;

    double PI = 3.141592;



    /*
     Creation des variables de sortie. On n'alloue la place sans en sp?cifier la valeur pour gagner du temps
    */


//     mexPrintf("Creation prob init");

    double prob_init[(*regime_sig)];


    for(i=0;i<(*regime_sig);i++)
    {
        prob_init[i] = 0;
    }
    prob_init[0] = 1;



    taille_current = (*taille);



//     mexPrintf("Debut Forward");
    log_like = 0;
    transition = 0;

for (seg=0; seg < (*number_ts_segments); seg++) {
    for(t=0;t<(seg_length[seg]);t++)
    {
       //mexPrintf("Iteration %d",t);
       maxi = -1e8;
       for(i=0;i<(*regime_sig);i++)
       {
           if(t>0)
           {
               transition = 0;
               for(z=0;z<(*regime_sig);z++)
               {
                   if((*dependance)==1 && d_t[seg_start_index[seg]-1+t]<(*dist))
                   {
                       transition = transition + p_lambda[z + i*(*regime_sig)]*forward[(seg_start_index[seg]-1+t-1) + z*taille_current];
                   }
                   else
                   {
                       transition = transition + p_sig[z + i*(*regime_sig)]*forward[(seg_start_index[seg]-1+t-1) + z*taille_current];
                   }
               }
           }
           else
           {
               transition = p_sig[0 + i*(*regime_sig)];
           }
           sig_current = sigma[i];
           eps_sq = eps_t[seg_start_index[seg]-1+t]*eps_t[seg_start_index[seg]-1+t];
           if(transition!=0)
           {
                forward[seg_start_index[seg]-1+t + taille_current*i] = (*temper)*(-0.5*log(2*PI*sig_current)-0.5*eps_sq/sig_current) + log(transition);
                if(maxi<forward[seg_start_index[seg]-1+t + taille_current*i])
                {
                   maxi = forward[seg_start_index[seg]-1+t + taille_current*i];
                }
           }
           else
           {
               forward[seg_start_index[seg]-1+t + taille_current*i] = 0;
           }


       }

       Q = 0;
       for(i=0;i<(*regime_sig);i++)
       {
           if(forward[seg_start_index[seg]-1+t + taille_current*i]!=0)
           {
               forward[seg_start_index[seg]-1+t + taille_current*i] = exp(forward[seg_start_index[seg]-1+t + taille_current*i]-maxi);
               Q = Q + forward[seg_start_index[seg]-1+t + taille_current*i];
           }
       }

       for(i=0;i<(*regime_sig);i++)
       {
           forward[seg_start_index[seg]-1+t + taille_current*i] = forward[seg_start_index[seg]-1+t + taille_current*i]/Q;
//             mexPrintf("  forward[%d,%d] = %f",t + taille_current*i,i,forward[t + taille_current*i]);
       }
//        mexPrintf("\n");
//        for(i=0;i<(*regime_sig);i++)
//        {
//            mexPrintf("  Epsilon[%d,%d] = %f",t + taille_current*i,i,epsilon_Haas[(deb+t) + fin*i]);
//        }
//        mexPrintf("\n");
       log_like = log_like + log(Q) + maxi;


//        mexPrintf("log_dens = %f\n",log_like);
//        mexPrintf("maxi = %f\n",Q);
//             mexPrintf("theta[%d] = %f\n",1,theta[1]);
//             mexPrintf("theta[%d] = %f\n",2,theta[2]);
//             mexPrintf("theta[%d] = %f\n",3,theta[3]);
    }
}
    (*log_like_out) = log_like;


    /*
      Procedure Backward
    */
    /*
     Creation d'une variable. On n'alloue la place sans en sp?cifier la valeur pour gagner du temps
    */

    double back_trans[(*regime_sig)];

//     mexPrintf("Debut Back");

for (seg=0; seg < (*number_ts_segments); seg++) {

    unif = Unif_rand[seg_end_index[seg]-1];
    Q = forward[seg_end_index[seg]-1];
    i = 0;

    while(Q<unif) {
        i = i+1;
        Q = Q + forward[seg_end_index[seg]-1 + i*taille_current];
    }
    sn_move[seg_end_index[seg]-1] = i;

    for(t=(seg_length[seg])-2;t>=0;t--) {
        maxi = -1e9;
        q = (int) sn_move[seg_start_index[seg]-1+t+1];
//         mexPrintf("   sn_move[%d] = %f   sn_stay = %f   ",t+1,sn_move[t+1],sn_current[t+1]);
        Q = 0;
        for(z=0;z<(*regime_sig);z++)
        {
            if((*dependance)==1 && d_t[seg_start_index[seg]-1+t]<(*dist))
            {
                back_trans[z] = forward[seg_start_index[seg]-1+t + z*taille_current]*p_lambda[z + q*(*regime_sig)];
            }
            else
            {
                back_trans[z] = forward[seg_start_index[seg]-1+t + z*taille_current]*p_sig[z + q*(*regime_sig)];
            }
            Q = Q + back_trans[z];

        }
//          mexPrintf("\n");
        for(i=0;i<(*regime_sig);i++)
        {
            back_trans[i] = back_trans[i]/Q;
//             mexPrintf("backward[%d,%d] = %f",taille_current-1-incr,i,back_trans[i]);
        }
//         mexPrintf("\n");
        unif = Unif_rand[seg_start_index[seg]-1+t];
        Q = back_trans[0];
        i = 0;
        while(Q<unif) {
            i = i+1;
            Q = Q + back_trans[i];
        }
        sn_move[seg_start_index[seg]-1+t] = i;


    }
}



    return;
}
