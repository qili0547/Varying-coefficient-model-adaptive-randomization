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
void Ua(double U[],double Jac[] ,double temp[],double Xw[],double T[],double t,double h ,int gamma_n,int n,int tn,double m[],double temp1[])
{
    int i,j,k,l;
    double u ;

    for(i=0;i<n;i++)
        for(j=0;j<m[i];j++)
    { u=fabs(T[i+n*j]-t)/h;  
      if(u<1)
       for(k=0;k<2*gamma_n;k++)
          {
            U[k]+=temp[i+n*j]*Xw[i+n*(tn*k+j)]*(1-u*u)*0.75/h;
            for(l=0;l<2*gamma_n;l++)
            {
               Jac[k+gamma_n*2*l]-=temp1[i+n*j]*Xw[i+n*(tn*k+j)]*Xw[i+n*(tn*l+j)]*(1-u*u)*0.75/h  ;
            }
          }
      
  }
}

/***************************************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *temp,*temp1,*Xw,*T,t,h,*U,*Jac,*m;
 int tn , n,gamma_n;
  /* Check for proper number of arguments. */
  if(nrhs!=9) {
    mexErrMsgTxt("Nine input required.");
  } else if(nlhs>2) {
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
temp= mxGetPr(prhs[0]);
Xw = mxGetPr(prhs[1]);
T= mxGetPr(prhs[2]);
t=mxGetScalar(prhs[3]);
h=mxGetScalar(prhs[4]);
gamma_n=(int)mxGetScalar(prhs[5]);
tn=(int)mxGetScalar(prhs[6]);
m= mxGetPr(prhs[7]);
temp1= mxGetPr(prhs[8]);
 
 /* if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Input must be a noncomplex scalar double.");
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(2*gamma_n,1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(2*gamma_n,2*gamma_n, mxREAL);
  /* Assign pointers to each input and output. */
/*/[N,T,Y,Xwave,m]=poison(X,R,C,tao,alpha0,alpha1,gamma);*/

U= mxGetPr(plhs[0]);
Jac= mxGetPr(plhs[1]);
   /* Call the timestwo subroutine. */
Ua(U,Jac,temp,Xw,T,t,h,gamma_n,n,tn,m,temp1);
   /* Destroy array */
  /*/  mxDestroyArray( plhs[0]);*/
}
