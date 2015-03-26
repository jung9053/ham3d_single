// ##################################################################
//
// updateSoln.c
//
// Update the solution scheme (1st order EE)
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
void updateSoln(double *qsrc,double *qdest, double *res,double *sigma,
	double dt,double CFL,double coef,int ncells)
{
  int i,j,m;
  double dtfac;
  //
  // 1st order euler explicit for now
  //
  m=0;
  //
  for(i=0;i<ncells;i++)
    {
      for(j=0;j<NVAR;j++)
	{
     dtfac=coef*CFL/sigma[i];
	  qdest[m]=qsrc[m]+dtfac*res[m];
	  m++;
//	  tracef(dtfac);
	}
    }
}
     
  

