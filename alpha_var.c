#include "math.h"
#include "mex.h"
#include <time.h> 
/*
 * timestwo.c - example found in API guide
 *
 * Computational function that takes a scalar and double it.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2000 The MathWorks, Inc.
 */
 
/* $Revision: 1.8 $ */
void sig_var(double Y[],double mu_hat[],double X[],double T[],double dN[],double h,double t,int gamma_n,int tn,double m[],double Sigma[],int n,double IExx[],double Exz[],double InvA,double temp00[],int  beta_n)
{
    int i,j,k,l;
    double u, *temp;
    mxArray *temp0[1];
   temp0[0]=mxCreateDoubleMatrix(gamma_n,1, mxREAL);
   temp=mxGetPr(temp0[0]);
    for(i=0;i<n;i++)
   { for(k=0;k<gamma_n;k++)temp[k]=0; 
       for(j=0;j<m[i];j++)
        { u=fabs(T[i+n*j]-t)/h;  
          if(u<1)
           for(k=0;k<gamma_n;k++)
              for(l=0;l<gamma_n;l++)
                  temp[k]+=IExx[k+l*gamma_n]*X[i+n*(tn*l+j)]*(1-u*u)*0.75/h*(Y[i+n*j]-mu_hat[i+n*j])*dN[i+j*n] ;
      }/*****************/
     for(k=0;k<gamma_n;k++) 
               for(l=0;l<gamma_n;l++)for(j=0;j<beta_n;j++)   
                     temp[k]-=IExx[k+l*gamma_n]*Exz[(l+gamma_n*j)]*temp00[i+n*j];
    for(k=0;k<gamma_n;k++)
              for(l=0;l<gamma_n;l++)
                   Sigma[k+gamma_n*l]+=temp[k]*temp[l];
   } 
   
   mxDestroyArray( temp0[0]);
}

/***************************************************************/
 
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *Y,*mu_hat,*X,*T,*dN,h,t,*Sigma,*m,*IExx,*Exz,InvA,*temp00;
 int tn , n,gamma_n,beta_n;
  /* Check for proper number of arguments. */
  if(nrhs!=14) {
    mexErrMsgTxt("Seven input required.");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }
  
      /* Check data type of input argument */
  for(tn=0;tn<nrhs;tn++)
  {
  if (!(mxIsDouble(prhs[tn]))) 
  { 
      mexErrMsgTxt("Input array must be of type double.");
  } 
  }
    /***********************************/


 n = mxGetM(prhs[0]);
 beta_n=mxGetN(prhs[13]);
Y= mxGetPr(prhs[0]);
mu_hat= mxGetPr(prhs[1]);
X= mxGetPr(prhs[2]);
T= mxGetPr(prhs[3]);
dN= mxGetPr(prhs[4]);
h=mxGetScalar(prhs[5]);
t=mxGetScalar(prhs[6]);
gamma_n=(int)mxGetScalar(prhs[7]);
tn=(int)mxGetScalar(prhs[8]);
m= mxGetPr(prhs[9]);
IExx= mxGetPr(prhs[10]);
Exz= mxGetPr(prhs[11]);
InvA=mxGetScalar(prhs[12]);
temp00= mxGetPr(prhs[13]);


 /* if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Input must be a noncomplex scalar double.");
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(gamma_n,gamma_n, mxREAL);
  /* Assign pointers to each input and output.
//[N,T,Y,Xwave,m]=poison(X,R,C,tao,alpha0,alpha1,gamma); */

Sigma= mxGetPr(plhs[0]);
   /* Call the timestwo subroutine. */
sig_var(Y,mu_hat,X,T,dN,h,t,gamma_n,tn,m,Sigma,n,IExx,Exz,InvA,temp00, beta_n);
   /* Destroy array
  //  mxDestroyArray( plhs[0]); */
}
