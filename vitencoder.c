#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *u,*S,*y,*valueY;
    mxArray *Y;
    int s,n,k;

    u = mxGetPr(prhs[0]);
    Y = mxDuplicateArray(prhs[1]);
    S = mxGetPr(prhs[2]);    
    n = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(2,n, mxREAL);
    y = mxGetPr(plhs[0]);

    s = 0;
    for(k=0;k<n;k++)
    {        
        valueY = mxGetPr(mxGetCell(Y,(int)u[k]+(s<<1)));         
        y[k<<1] = valueY[0];
        y[1+(k<<1)] = valueY[1];
        s = (int)S[(int)u[k]+(s<<1)];
    }
}

