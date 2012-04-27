#include "mex.h"
#include "matrix.h"

void encoder(double *u,double*Y,double*S,int n,int Ysize,double *y)
{
    int k,uidx,s = 0;
    
    for(k=0;k<n;k++)
    {        
        uidx=u[k];
        y[k<<1] = Y[uidx+(s<<1)];
        y[1+(k<<1)] = Y[uidx+(s<<1)+Ysize];
        s = (int)S[uidx+(s<<1)];
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *u,*S,*Y,*y;
    int n,Ysize;

    u = mxGetPr(prhs[0]);
    Y = mxGetPr(prhs[1]);
    S = mxGetPr(prhs[2]);
    
    n = mxGetN(prhs[0]);       /*Input size*/
    Ysize = mxGetN(prhs[1]);   /* A 2 x (2^nu) x 2 matrix is seen as a 2 x (2x2^nu) matrix in c
                                  Ysize is the size of a single 2 x (2^nu) matrix (that is 2x2^nu)*/
    
    plhs[0] = mxCreateDoubleMatrix(2,n, mxREAL);
    y = mxGetPr(plhs[0]);
    
    encoder(u,Y,S,n,Ysize,y);
}

