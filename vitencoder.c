#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *uIn,*S,*y,*valueY;
    mxArray *Y,*N,*cellY;
    int u,s,n,k;

    uIn = mxGetPr(prhs[0]);
    S = mxGetPr(prhs[2]);
    Y = mxDuplicateArray(prhs[1]);
    n = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(2,n, mxREAL);
    y = mxGetPr(plhs[0]);

    s = 0;
    for(k=0;k<n;k++)
    {
        u = uIn[k];
        cellY = mxGetCell(Y,u+(s<<1));
        valueY = mxGetPr(cellY);
        y[k<<1] = valueY[0];
        y[1+(k<<1)] = valueY[1];
        s = (int)S[u+(s<<1)];
    }

}

