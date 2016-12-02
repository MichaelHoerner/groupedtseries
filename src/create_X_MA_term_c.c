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


/*  Fonction pour la densite d'un modele ARMA : 8 input : la serie y (T*1), les parametres GARCH, param?tres mu, cass , eta
    L'output : la densite et la variance conditionnelle
 *[log_dens sigma] = dens_GJR_GARCH_MCMC(y,theta_GARCH,theta_mu,cass,eta)
 *[log_dens] = dens_CP_ARMA_c++(y,X,AR_lags,MA_lags,theta,sigma,tau)
 */

void create_X_MA_term_c(double *y, double *X, int *AR_lags, int *MA_lags, double *theta, double *sn, int *taille, double *err_out, double *X_out)
{

//declare variables
    double *log_like,*err;//,*AR_lags,*MA_lags;
    int dimx, dimy, numdims;
    int i,q,z,iter,etat,index,ind_mu,taille_ARMA;
    double pi,eps,eps_sq;



//figure out dimensions
    dimy = (*taille);
    dimx = (*taille);



//associate outputs
    //log_like = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    //f_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

    taille_ARMA = 1+(*AR_lags)+(*MA_lags);
    //g_out_m = plhs[1] = mxCreateDoubleMatrix(dimy,taille_ARMA,mxREAL);
//      mexPrintf("taille_ARMA = %d\n",taille_ARMA);
//      mexPrintf("tau1 = %d\n",tau[0]);
//      mexPrintf("sigma1 = %f\n",sigma[0]);
//
//
//      mexPrintf("X 0 = %f\n",X[0]);
//      mexPrintf("X 1 = %f\n",X[1]);
//      mexPrintf("X taille+1 = %f\n",X[taille+1]);
//      mexPrintf("X taille+2 = %f\n",X[taille+2]);

//     mexPrintf("theta_all[%d] = %f\n",0,theta_all[0]);
//
//     mexPrintf("dimy = %d\n",dimy);
//     mexPrintf("dimx = %d\n",dimx);

//     mxDestroyArray(a_in_m);
//     mxDestroyArray(b_in_m);
//     mxDestroyArray(c_in_m);
//     mxDestroyArray(d_in_m);

//    log_like = mxGetPr(e_out_m);

    //err = mxGetPr(f_out_m);
    //X_out = mxGetPr(g_out_m);

    pi = 3.14159265;

    etat = sn[0]-1;
    index = (etat)*taille_ARMA;

        eps = 0;
        for(q=0;q<taille_ARMA;q++)
        {
          eps = eps + X[q*(*taille)+0]*theta[index+q];
          X_out[q*(*taille)+0] = X[q*(*taille)+0];
        }
        eps = y[0]-eps;
        err_out[0] = eps;
        if(MA_lags>0)
        {
            X[(1+(*AR_lags))*(*taille)+1] = eps;
            X_out[(1+(*AR_lags))*(*taille)+1] = eps;
        }

    //     mexPrintf("log_dens = %f\n",log_dens);
//             mexPrintf("theta[%d] = %f\n",0,theta[0]);
//             mexPrintf("theta[%d] = %f\n",1,theta[1]);
//             mexPrintf("theta[%d] = %f\n",2,theta[2]);
//             mexPrintf("theta[%d] = %f\n",3,theta[3]);
// //            mexPrintf("sigma_sq[%d] = %f\n",0,sigma_sq[0]);
//            mexPrintf("eps = %f\n",eps);
//            mexPrintf("y = %f\n",y[0]);
// //           mexPrintf("mu = %f\n",mu[ind-1]);
//            mexPrintf("log_dens = %f\n",log_dens);
//          system("pause");

    //do something
        for(i=1;i<(*taille);i++)
        {
            etat = sn[i]-1;
            index = (etat)*taille_ARMA;


            eps = 0;
            for(q=0;q<taille_ARMA;q++) {
                eps = eps + X[q*(*taille)+i]*theta[index+q];
                X_out[q*(*taille)+i] = X[q*(*taille)+i];
            }
            eps = y[i]-eps;
            if((*MA_lags)>0)
            {
                if(i+1<(*taille))
                {
                    X[(1+(*AR_lags))*(*taille)+i+1] = eps;
                    if((*MA_lags)>1)
                    {
                      if(i+1<(*MA_lags))
                      {
                          iter = i;
                      }
                      else
                      {
                          iter = (*MA_lags)-1;
                      }
                      for(z=0;z<iter;z++)
                      {
                          X[(1+(*AR_lags)+1+z)*(*taille)+i+1] = X[(1+(*AR_lags))*(*taille)+i-z];
                      }
                    }
                }
            }
            err_out[i] = eps;



//             if(sigma_sq[i]<0)
//             {
//                 mexPrintf("sigma[%d] = %f\n",i,sigma_sq[i]);
//             mexPrintf("log_dens[%d] = %f\n",i,log_dens);
//             mexPrintf("log_dens[%d] = %d\n",i,index);
//             mexPrintf("theta[%d] = %f\n",index+0,theta[index+0]);
//             mexPrintf("theta[%d] = %f\n",index+1,theta[index+1]);
//             mexPrintf("theta[%d] = %f\n",index+2,theta[index+2]);
//             mexPrintf("theta[%d] = %f\n",index+3,theta[index+3]);
//             }
        }


    return;
}
