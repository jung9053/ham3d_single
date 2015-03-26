/**
   Koren's differentiable limiter for
   3rd order accurate reconstruction
*/
#include <stdio.h>
#include <stdlib.h>

void muscld(double **f,
	    double **ql,
	    double **qr,
	    double **f2,
	    int is,
	    int ie,
	    double th,
	    double qt,
	    double eps,
	    int imax,
	    int nq)
{  
  int i,n,im;
  double thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt;
  
  im=imax-1;
  if (qt == 0.0)
    {
      for(n=0;n<nq;n++)
	{
	  for(i=is;i<=ie;i++)
	    {
	      ql[i][n]=f[i][n];
	      qr[i][n]=f[i][n];
	    }
	}
       return;
    }
  
  /*
   * do 3rd order otherwise
  */
  thm    = 1.0 - th;
  thp    = 1.0 + th;
  
  for(n=0;n<nq;n++)
    {
      for(i=0;i<im;i++)
	f2[i][n]=f[i+1][n]-f[i][n];
      for(i=is;i<=ie;i++)
	{
	  f2i    = f2[i][n];
	  f2i1   = f2[i-1][n];
	  a1     = 3.0*(f2i*f2i1+eps);
	  a2     = 2.0*(f2i-f2i1)*(f2i-f2i1) + a1;
	  f3i    =  a1/a2 ;
	  f3qt   = qt*f3i;
	  ql[i][n]= f[i][n]+f3qt*( thm*f2i1 + thp*f2i );
	  qr[i][n]= f[i][n]-f3qt*( thp*f2i1 + thm*f2i );
	}
    }
}

