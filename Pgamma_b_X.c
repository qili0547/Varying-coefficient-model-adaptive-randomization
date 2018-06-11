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
double *MatrixOpp(double *a,int n);
double psy(double x, int link_id);
/* $Revision: 1.8 $ */

void Ub(double Pgx[],double gamma_hat1[],double X[],double V[],double Y[],double T[],double dN[],int tn,double m[],double h,double zeta00[],double dzeta00[],int n,int N,int gamma_n,double T0[],int ID,double IExx[],double Exz[],int link_id,int beta_n)
{
   int i,j,k,l,id1,id2 ;
   double u, t, mu_a=1,*Pgb,*temp0,*temp1,temp00,temp01;
   mxArray *temp[2];
   temp[0]=mxCreateDoubleMatrix(gamma_n,gamma_n, mxREAL);
   temp[1]=mxCreateDoubleMatrix(gamma_n,beta_n, mxREAL);
   temp0=mxGetPr(temp[0]);
   temp1=mxGetPr(temp[1]);
 
 if (ID==0) 
   {  
   for(id1=0;id1<n;id1++)
   for(id2=0;id2<m[id1];id2++)
     {  
	 t=T[id1+n*id2];
     for(k=0;k<gamma_n;k++)
        {    
        for(l=0;l<beta_n;l++)  temp1[k+gamma_n*l]=0;
        for(l=0;l<gamma_n;l++) temp0[k+gamma_n*l]=0;
        }
     for(i=0;i<n;i++)
     for(j=0;j<m[i];j++)
        { 
		u=fabs(T[i+n*j]-t)/h;  
        if(fabs(u)<1)
          {
   		   mu_a=0;
           for(k=0;k<gamma_n;k++) 
		       mu_a+=gamma_hat1[i+n*(tn*k+j)]*X[i+n*(tn*k+j)];
           temp01=mu_a; 
		   for(k=0;k<beta_n;k++) 
			    temp01+=zeta00[i+n*(tn*k+j)]*V[i+n*(tn*k+j)];
           mu_a=psy(temp01,link_id);				
           for(k=0;k<gamma_n;k++)
               {    
                for(l=0;l<beta_n;l++)
                   temp1[k+gamma_n*l]-=mu_a*X[i+n*(tn*k+j)]*dzeta00[i+n*(tn*l+j)]*V[i+n*(tn*l+j)]*(1-u*u)*0.75/h*dN[i+n*j];
                for(l=0;l<gamma_n;l++)
                   temp0[k+gamma_n*l]+=mu_a*X[i+n*(tn*k+j)]*X[i+n*(tn*l+j)]*(1-u*u)*0.75/h*dN[i+n*j];
                }
            }
        }    
    Pgb=MatrixOpp(temp0,gamma_n);
 
   for(j=0;j<beta_n;j++)  
   for(k=0;k<gamma_n;k++)  
   for(l=0;l<gamma_n;l++)
       Pgx[id1+n*(tn*j+id2)]+=temp1[k+gamma_n*j]*Pgb[k+gamma_n*l]*X[id1+n*(tn*l+id2)];
     }
   }

  else
  for(id1=0;id1<N;id1++)   
    {   
    t=T0[id1];
   for(k=0;k<gamma_n;k++)  
   for(l=0;l<gamma_n;l++)
       temp0[k+gamma_n*l]=0;    
   for(i=0;i<n;i++)
   for(j=0;j<m[i];j++)
        { u=fabs(T[i+n*j]-t)/h;           
          if(fabs(u)<1)
          {  mu_a=0; 
            for(k=0;k<gamma_n;k++) 
			    mu_a+=gamma_hat1[i+n*(tn*k+j)]*X[i+n*(tn*k+j)];
            temp01=mu_a; 
			for(k=0;k<beta_n;k++) 
			    temp01+=zeta00[i+n*(tn*k+j)]*V[i+n*(tn*k+j)];
            mu_a=psy(temp01,link_id);

            for(k=0;k<gamma_n;k++)
               {
  			   for(l=0;l<beta_n;l++)
                   Exz[id1+N*(k+gamma_n*l)]+=mu_a*X[i+n*(tn*k+j)]*dzeta00[i+n*(tn*l+j)]*V[i+n*(tn*l+j)]*(1-u*u)*0.75/h*dN[i+n*j]/n;  
                for(l=0;l<gamma_n;l++)
                   temp0[k+gamma_n*l]+=mu_a*X[i+n*(tn*k+j)]*X[i+n*(tn*l+j)]*(1-u*u)*0.75/h*dN[i+n*j]/n;
                } 
          }
        }
      Pgb=MatrixOpp(temp0,gamma_n);
 
     for(k=0;k<gamma_n;k++) 
             for(l=0;l<gamma_n;l++) 
             IExx[id1+N*(k+l*gamma_n)]=Pgb[k+l*gamma_n];
               
  }

    


mxDestroyArray(temp[0]);
mxDestroyArray(temp[1]);
free(Pgb);



}
/***************************************************************/
 
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *gamma_hat1,*Pgx,*X,*V,*Y,*T,*dN,*m,*zeta00,*dzeta00,h,*T0,*IExx,*Exz;
 int tn , n,gamma_n,N,ID, link_id,beta_n;
  /* Check for proper number of arguments. */
  if(nrhs!=14) {
    mexErrMsgTxt("twelve input required.");
  } else if(nlhs>3) {
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
    /************************************/

N = mxGetM(prhs[11]);
n=mxGetM(prhs[1]);
gamma_hat1= mxGetPr(prhs[0]);
X = mxGetPr(prhs[1]);
V= mxGetPr(prhs[2]);
Y= mxGetPr(prhs[3]);
T= mxGetPr(prhs[4]);
dN= mxGetPr(prhs[5]);
tn=(int)mxGetScalar(prhs[6]);
m= mxGetPr(prhs[7]);
h=mxGetScalar(prhs[8]);
zeta00=mxGetPr(prhs[9]);
dzeta00=mxGetPr(prhs[10]);
T0= mxGetPr(prhs[11]);
ID=(int)mxGetScalar(prhs[12]);
gamma_n=(int)(mxGetN(prhs[0])/2/tn);
beta_n=(int)(mxGetN(prhs[9])/tn);
 link_id=(int)mxGetScalar(prhs[13]);

 /* if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Input must be a noncomplex scalar double.");
  }
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(n,tn*beta_n, mxREAL);
/*  plhs[1] = mxCreateDoubleMatrix(N,gamma_n*gamma_n,mxREAL);
//  plhs[2] = mxCreateDoubleMatrix(N,gamma_n, mxREAL);*/
   plhs[1] = mxCreateDoubleMatrix(N,gamma_n*gamma_n,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(N,gamma_n*beta_n, mxREAL);
 /* plhs[2] = mxCreateDoubleMatrix(gamma_n,gamma_n, mxREAL);


  /* Assign pointers to each input and output. 
//[N,T,Y,Xwave,m]=poison(X,R,C,tao,alpha0,alpha1,gamma);*/


Pgx= mxGetPr(plhs[0]);
IExx= mxGetPr(plhs[1]);
Exz= mxGetPr(plhs[2]);
/*/Pgx1= mxGetPr(plhs[2]);
     Call the timestwo subroutine. */
Ub(Pgx,gamma_hat1,X,V,Y,T,dN,tn,m,h,zeta00,dzeta00,n,N,gamma_n,T0,ID,IExx,Exz, link_id,beta_n);
   /* Destroy array */
  /*/  mxDestroyArray( plhs[0]);*/
}



/*/%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/*/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double *MatrixOpp(double *A,int n)
{ 
int *is,*js,i,j,k,l,u,v;
    double d,p,*a;
    is=(int*)malloc(n*sizeof(int));
    js=(int*)malloc(n*sizeof(int));
    a=(double*)malloc(n*n*sizeof(double));
    for (i=0; i<=n-1; i++)
        for (j=0; j<=n-1; j++)
          a[i*n+j]= A[i*n+j];
    /*********************/
    for (k=0; k<=n-1; k++)
      { d=0.0;
        for (i=k; i<=n-1; i++)
        for (j=k; j<=n-1; j++)
          { l=i*n+j; p=fabs(a[l]);
            if (p>d) { d=p; is[k]=i; js[k]=j;}
          }
         if (d+1.0==1.0)
          { free(is); free(js); /*/printf("err**not inv\n");
           // return NULL;*/
           for (i=0; i<=n-1; i++)
           for (j=0; j<=n-1; j++)
               a[i*n+j]=0;
          return a;
          } 
        if (is[k]!=k)
          for (j=0; j<=n-1; j++)
            { u=k*n+j; v=is[k]*n+j;
              p=a[u]; a[u]=a[v]; a[v]=p;
            }
        if (js[k]!=k)
          for (i=0; i<=n-1; i++)
            { u=i*n+k; v=i*n+js[k];
              p=a[u]; a[u]=a[v]; a[v]=p;
            }
        l=k*n+k;
        a[l]=1.0/a[l];
        for (j=0; j<=n-1; j++)
          if (j!=k)
            { u=k*n+j; a[u]=a[u]*a[l];}
        for (i=0; i<=n-1; i++)
          if (i!=k)
            for (j=0; j<=n-1; j++)
              if (j!=k)
                { u=i*n+j;
                  a[u]=a[u]-a[i*n+k]*a[k*n+j];
                }
        for (i=0; i<=n-1; i++)
          if (i!=k)
            { u=i*n+k; a[u]=-a[u]*a[l];}
      }
    for (k=n-1; k>=0; k--)
      { if (js[k]!=k)
          for (j=0; j<=n-1; j++)
            { u=k*n+j; v=js[k]*n+j;
              p=a[u]; a[u]=a[v]; a[v]=p;
            }
        if (is[k]!=k)
          for (i=0; i<=n-1; i++)
            { u=i*n+k; v=i*n+is[k];
              p=a[u]; a[u]=a[v]; a[v]=p;
            }
      }
    free(is); free(js);
    return a;
}


double psy(double x, int link_id)
{/*calculate the dirivative of psy(x)*/
    double y;
    
 if( link_id==0 )  y=exp(x);
 if(link_id==1) y=1;
 if(link_id==2) y=exp(x)/((1+exp(x))*(1+exp(x)));
    return y;
}