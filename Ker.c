#include "math.h"
#include "mex.h"

/*
 * timestwo.c - example found in API guide
 *
 * Computational function that takes a scalar and double it.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2000 The MathWorks, Inc.
 */
 
/* $Revision: 1.8 $ */

void Kernel(double alpha_hat00[],double alpha_hat1[],double T[],double T0[],int alpha_n,int tn,double h,int n, double m[],double m2[],int n2)
{
  int i,j,k,id1,id2;
  double u,K;
  for (id1=0;id1<n2;id1++)
    for (id2=0;id2<m2[id1];id2++)	
	{
	 K=0;
     for(i=0;i<n;i++)     
     for(j=0;j<m[i];j++) 
        {
		 if (T[i+n*j]>=0)
     		 u=(T[i+n*j]-T0[id1+n2*id2])/2/h; 
	     if (T[i+n*j]<0) 
		     u=1;
         if (u>=1) break;
         if(fabs(u)<1)
            {   
		     K+=(1-u*u);
             for(k=0;k<2*alpha_n;k++)
                 alpha_hat00[id1+n2*(tn*k+id2)]+=(1-u*u)*alpha_hat1[i+n*(tn*k+j)];
            }
        }
			 for(k=0;k<2*alpha_n;k++) 
	         if(K>0) 
		     alpha_hat00[id1+n2*(tn*k+id2)]/=K;
            }
    }

/***************************************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *alpha_hat1,*alpha_hat00,*T,*T0,h,*m,*m2;
 int  n,tn,alpha_n,n2;
  /* Check for proper number of arguments. */
  if(nrhs!=8) {
    mexErrMsgTxt("Seven pinput required.");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }
  
      /* Check data type of input argument */
  for(n=0;n<nrhs;n++)
  {
  if (!(mxIsDouble(prhs[n]))) 
  { 
      mexErrMsgTxt("Input array must be of type double.");
  } 
  }
    /************************************/
 n= mxGetM(prhs[0]);
  n2=mxGetM(prhs[2]);
alpha_hat1= mxGetPr(prhs[0]);
T= mxGetPr(prhs[1]);
T0= mxGetPr(prhs[2]);
alpha_n = (int)mxGetScalar(prhs[3]);
tn =(int) mxGetScalar(prhs[4]);
h =  mxGetScalar(prhs[5]); 
m= mxGetPr(prhs[6]);
m2= mxGetPr(prhs[7]);

 /* if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Input must be a noncomplex scalar double.");
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(n2,2*alpha_n*tn, mxREAL);
  /* Assign pointers to each input and output. */


alpha_hat00= mxGetPr(plhs[0]);
   /* Call the timestwo subroutine. */
Kernel(alpha_hat00,alpha_hat1,T,T0,alpha_n,tn,h,n,m,m2,n2); 
   /* Destroy array */

}
