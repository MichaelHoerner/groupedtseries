

#include <R.h>
#include <math.h>
#include <float.h>


void dens_CP_ARMA_c(double *y, double *X, int *AR_lags, int *MA_lags, double *theta, double *sigma, double *sn, double *sn_sig, int *deb, int *taille, double *log_dens, double *eps_out)
{

//declare variables
    //mxArray *e_out_m;
    //const dims;
    int dimx, dimy, numdims;
    int i,q,z,iter,etat,etat_sig,index,ind_mu,taille_ARMA;
    double pi,eps,eps_sq,temper;

//figure out dimensions
    /*dims = mxGetDimensions(X);*/
    dimy = (*taille); dimx = (*taille);

    // Create second
//associate outputs
    //log_like = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    /*e_out_m = mxCreateDoubleMatrix(1,1,mxREAL);*/
    //f_out_m = plhs[1] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

//associate pointers [log_dens] = log_dens_ARMA(y,X,AR_lags,MA_lags,theta,sigma,sn,deb)
    taille_ARMA = 1 + (*AR_lags) + (*MA_lags);
    (*deb) = (*deb)-1;
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

    pi = 3.14159265;

    etat = sn[(*deb)]-1;
    etat_sig = sn_sig[(*deb)]-1;
    index = (etat)*taille_ARMA;

        eps = 0;
        for(q=0;q<taille_ARMA;q++)
        {
          eps = eps + X[q*(*taille)+(*deb)]*theta[index+q];
        }
        eps = y[(*deb)]-eps;
        if((*MA_lags)>0)
        {X[(1+(*AR_lags))*(*taille)+1] = eps;}

        eps_out[(*deb)] = eps;
        eps_sq= eps*eps;
        (*log_dens) = -0.5*log(2*pi*sigma[etat_sig])-0.5*eps_sq/sigma[etat_sig];


        for(i=(*deb)+1;i<(*taille);i++)
        {

            etat = sn[i]-1;
            etat_sig = sn_sig[i]-1;
            index = (etat)*taille_ARMA;


            eps = 0;
            for (q=0;q<taille_ARMA;q++) {
                eps = eps + X[q*(*taille)+i]*theta[index+q];
            }
            eps = y[i]-eps;
            if((*MA_lags)>0)
            {
                if(i+1<(*taille))
                {
                    X[(1+(*AR_lags))*(*taille)+i+1] = eps;
                    //if((*MA_lags)>1)
                    //{
                    //  if(i+1<(*MA_lags))
                    //  {
                    //      iter = i;
                    //  }
                    //  else
                    //  {
                    //      iter = (*MA_lags)-1;
                    //  }
                    //  for(z=0;z<iter;z++)
                    //  {
                    //      X[(1+(*AR_lags)+1+z)*(*taille)+i+1] = X[(1+(*AR_lags))*(*taille)+i-z];
                    //  }
                    //}
                }
            }
            eps_out[i] = eps;
            eps_sq= eps*eps;
            (*log_dens) = (*log_dens) -0.5*log(2*pi*sigma[etat_sig])-0.5*eps_sq/sigma[etat_sig];

        }


        if((*log_dens)<-1e10 || (*log_dens)>1e10)
        {
            (*log_dens) = -1e10;
        }
        /* Pour le compilage sur visual c++ : include float.h et c'est la fonction _isnan
         if(_isnan(log_dens)==1)
         {
             log_dens = -1e10;
         }*/
        /* Pour le compilage sur SDK : include math.h et c'est la fonction isnan */
         if(isnan((*log_dens))==1)
         {
             (*log_dens) = -1e10;
         }


//          if(isnan(log_dens)==1 || isinf(log_dens)==1)
//          {
//              log_dens = -1e10;
//          }


    //log_like[0] = (*log_dens);

    return;
}
