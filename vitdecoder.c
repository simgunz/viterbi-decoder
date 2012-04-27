#include "mex.h"
#include "matrix.h"
#include "math.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define M(a) ((a) == (0) ? (a-1) : (a))

void decoder(double *r,double*Y,double*S,double*N,int n,int Ysize,int ns,double *u_out)
{
    int i,j,k,z;
    int bb;
    double tblen,inf;
    double gamma[ns][2];
    int survivors[ns][n][2];
    
    tblen = 5*log2(ns);
    inf = -1000;
                    
    for(z=0;z<ns;z++)
    {
        gamma[z][0] = inf;
        gamma[z][1] = inf;
    }
        
    gamma[0][0] = 0;
 
    /*mxCalloc*/
                    
    for(i=0;i<n;i++)
    {
        int maxgamma[2] = {0,0};
        for(j=0;j<ns;j++)
        {                                          
            int u1 = (int)N[j];
            int s1 = (int)N[j+Ysize];

            int u2 = (int)N[j+ns];
            int s2 = (int)N[j+ns+Ysize];

            double tempgamma[2];
            tempgamma[0] = gamma[(int)N[j+Ysize]][i%2];
            tempgamma[1] = gamma[(int)N[j+ns+Ysize]][i%2];


            for(z=0;z<2;z++)
            {                    
                tempgamma[0] += r[z+(i<<1)]*M((int)Y[u1+(s1<<1)+(z*Ysize)]);
                tempgamma[1] += r[z+(i<<1)]*M((int)Y[u2+(s2<<1)+(z*Ysize)]);                    
            }
            
            bb = (tempgamma[0] > tempgamma[1]) ? 0 : 1;
               
            gamma[j][(i+1)%2] = tempgamma[bb] - maxgamma[i%2];
            survivors[j][i][0] = (int)N[j+ns*bb];
            survivors[j][i][1] = (int)N[j+ns*bb+Ysize];
            
            if(gamma[j][(i+1)%2] > maxgamma[i%2])
                maxgamma[(i+1)%2] = gamma[j][(i+1)%2];
        }        
    }

    int s = 0;
    int u[n];
    for(z=n-1;z>=0;z--)
    {
        u[z] = survivors[s][z][0];
        s = survivors[s][z][1];
        /*mexPrintf("%d",survivors[s][z][1]);*/
    }
    for(z=0;z<n-log2(ns);z++)
    {
        u_out[z] = u[z];
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