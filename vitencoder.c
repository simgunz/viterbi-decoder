#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *u,*S,*Y,*y;
    int s,n,k,uidx,mu;
    mu = 8;
    u = mxGetPr(prhs[0]);
    Y = mxGetPr(prhs[1]);
    S = mxGetPr(prhs[2]);
    n = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(2,n, mxREAL);
    y = mxGetPr(plhs[0]);

    s = 0;
    for(k=0;k<n;k++)
    {        
        uidx=u[k];
        y[k<<1] = Y[uidx+(s<<1)];
        y[1+(k<<1)] = Y[8+uidx+(s<<1)];
        s = (int)S[uidx+(s<<1)];
    }
}

