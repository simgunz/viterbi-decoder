#include "mex.h"
#include "matrix.h"
#include "math.h"

#define M(a) ((a) == (0) ? (a-1) : (a))

void decoder(double *r,double*Y,double*S,double*N,int n,int Ysize,int ns,double *u_out)
{
    int i,j,k,z;
    int bb,u1,u2,s1,s2,s,bestJ,tblen;
    double gamma[ns][2],tempgamma[2],maxgamma[2] = {0,0};

    tblen = 5*log2(ns);
    
    int survivors[ns][tblen][2];

    for(z=0;z<ns;z++)
    {
        gamma[z][1] = gamma[z][0] = -10000;
    }

    gamma[0][0] = 0;

    for(i=0;i<n;i++)
    {
        for(j=0;j<ns;j++)
        {
            u1 = (int)N[j];
            s1 = (int)N[j+Ysize];

            u2 = (int)N[j+ns];
            s2 = (int)N[j+ns+Ysize];

            tempgamma[0] = gamma[s1][i%2];
            tempgamma[1] = gamma[s2][i%2];

            for(z=0;z<2;z++)
            {
                tempgamma[0] += r[z+(i<<1)]*M((int)Y[u1+(s1<<1)+(z*Ysize)]);
                tempgamma[1] += r[z+(i<<1)]*M((int)Y[u2+(s2<<1)+(z*Ysize)]);
            }

            bb = (tempgamma[0] > tempgamma[1]) ? 0 : 1;

            gamma[j][(i+1)%2] = tempgamma[bb] - maxgamma[i%2];
            survivors[j][i%tblen][0] = (int)N[j+ns*bb];
            survivors[j][i%tblen][1] = (int)N[j+ns*bb+Ysize];

            if(gamma[j][(i+1)%2] > maxgamma[i%2])
            {
                maxgamma[(i+1)%2] = gamma[j][(i+1)%2];
                bestJ = j;
            }
        }
        if(i>=tblen-1)
        {
            s = bestJ;
            for(z=0;z<tblen-1;z++)          /* Performance bottleneck*/
            {
                s = survivors[s][(i-z)%tblen][1];
            }
            u_out[i-tblen+1] = survivors[s][(i+1)%tblen][0];
        }
    }

    s = 0;
    i--;
    for(z=0;z<log2(ns);z++)
    {
        s = survivors[s][(i-z)%tblen][1];
    }
    for(z=log2(ns);z<tblen;z++)
    {
        u_out[i-z] = survivors[s][(i-z)%tblen][0];
        s = survivors[s][(i-z)%tblen][1];
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *r,*Y,*S,*N,*u_out;
    int n,Ysize,ns;

    r = mxGetPr(prhs[0]);
    Y = mxGetPr(prhs[1]);
    S = mxGetPr(prhs[2]);
    N = mxGetPr(prhs[3]);

    n = mxGetN(prhs[0]);
    Ysize = mxGetN(prhs[1]);
    ns = mxGetN(prhs[1])/2;     /* Number of possible states*/

    plhs[0] = mxCreateDoubleMatrix(1,n-log2(ns), mxREAL);
    u_out = mxGetPr(plhs[0]);

    decoder(r,Y,S,N,n,Ysize,ns,u_out);
}