// ##################################################################
//
// lineGaussSeidel.c
//
// Diagonally dominant line Gauss seidel method
//
// Written by Dr. Jayanarayanan Sitaraman
// Modified by Yong Su Jung
// ##################################################################

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define NQ 5

void lineGaussSeidel(GRID *g,SOLN *s,double cflnum,double dt)
{
  //
  int i,j,k,l,m,f,n,ii;
  int mm1,f1,f2,iface,idv,chainSize,idiff;
  int iface1,jface;
  int jleftCell,jrightCell;
  double ds[3];
  double lmat[NQ][NQ];
  double rmat[NQ][NQ];
  double dsnorm,nynx,nx2ny;
  int node1,node2,leftCell,rightCell,icell;
  double x1,y1,x2,y2,r01,r02,a1,a2,b1,b2,pp,dtfac;
  double ref[3][3];
  double linearl2rho,l2rho,l2rhodq,linearlinfrho;
  double damping;
  double tval[NQ];
  int isweep,ntotal;
  FILE *fp1,*fp2,*fp3,*fp4,*fp5;
  //
/*   fp1=fopen("A.dat","w"); */
/*   fp2=fopen("B.dat","w"); */
/*   fp3=fopen("C.dat","w"); */
/*   fp4=fopen("F.dat","w"); */
/*   fp5=fopen("X.dat","w"); */
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  //outputr(g,s);
  ntotal = g->ncells*NVAR;

  for(i=0;i<g->ncells;i++)
    {
      dtfac=cflnum/s->sigma[i];
      for(m=0;m<NVAR;m++)
	{
         s->r[NVAR*i+m]*=dtfac;
	 s->r0[NVAR*i+m]=s->r[NVAR*i+m];
	 s->dq[NVAR*i+m]=0;
	}
    }
  //
  // find diagonals
  //
  findDiagonals(g,s,cflnum,dt);
  //
  for (isweep=0;isweep<g->msweep;isweep++)
    {
      for (i=0;i<ntotal;i++) s->ddq[i]=0;
	//for(ii=0;ii<67;ii++)
      for(ii=0;ii<2*g->nchains-1;ii++)
	{
	  //	  outputdq(g,s);
	  i=(ii<g->nchains) ? ii: (2*g->nchains-2-ii);
	  idiff=(ii < g->nchains) ? 1: -1;
	  //
	  // zero out block matrices
	  //
	  for(m=0;m<g->nmaxchain+5;m++)
	    for(k=0;k<NQ;k++)
	      {
		g->F[m][k]=0;
		for(j=0;j<NQ;j++)
		  g->A[m][k][j]=g->B[m][k][j]=g->C[m][k][j]=0;
	      }
	  //
	  // collect cells on the loop
	  // 
	  f1=g->faceStartPerChain[i];
	  f2=g->faceStartPerChain[i+1];
          //
          // here idv == 1 means a full chain 
          //
	  idv=(g->chainConn[f1]==g->chainConn[f2-1]);
	  m=0;
	  chainSize=(f2-idv-f1);
	  for(f=f1;f<f2-idv;f++)
	    {
	      iface=g->chainConn[f];
	      leftCell=g->faces[8*iface+4];
	      rightCell=g->faces[8*iface+6];
	      //
	      for(j=0;j<NQ;j++)
		for(k=0;k<NQ;k++)
		  {
		    lmat[j][k]=(g->ff[iface]).lmat[j][k];
		    rmat[j][k]=(g->ff[iface]).rmat[j][k];
		  }
	      //
	      if (idv==1) {
		mm1=(m==0)?chainSize-1:m-1;
		for(j=0;j<NQ;j++)
		  {
		    for(k=0;k<NQ;k++)
		      {
			g->A[m][j][k]+=(rmat[j][k]);
			g->C[mm1][j][k]-=(lmat[j][k]);
		      }
		    g->F[m][j]=s->r[NVAR*leftCell+j];
		  }
	      }
	      else
		{
		  mm1=(m-1);
		  if (rightCell < 0 && idv==0 && f==f2-1) m--;
		  for(j=0;j<NQ;j++)
		    {
		      for(k=0;k<NQ;k++)
			{
			  if (mm1 > -1 && rightCell > -1)
			    {
			      g->A[m][j][k]+=(rmat[j][k]);
			      g->C[mm1][j][k]-=(lmat[j][k]);
			    }
			}
		      g->F[m][j]=s->r[NVAR*leftCell+j];
		    }
		}  
	      m++;
	    }     
	  m=0;
	  chainSize=chainSize-(idv==0);
	  for(f=f1;f<f2-1;f++)
	    {
	      iface=g->chainConn[f];
	      leftCell=g->faces[8*iface+4];
	      dtfac=cflnum/s->sigma[leftCell];
	      for(j=0;j<NQ;j++)
		for(k=0;k<NQ;k++)
		  {
		    g->A[m][j][k]*=dtfac;
		    g->B[m][j][k]=s->D[leftCell][j][k];
		    g->C[m][j][k]*=dtfac;
		  }
	      iface1=g->chainConn[f+1];
	      icell=leftCell;
	      for (j=0;j<6;j++)
		{
		  jface=g->c2f[6*icell+j];
		  if (jface != iface  &&
		      jface != iface1)
		    {
		      
		      jleftCell=g->faces[8*jface+4];
		      jrightCell=g->faces[8*jface+6];
            
            //printf("icell:%d leftcell:%d rightcell:%d\n",icell,jleftCell,jrightCell);

		      if (icell==jleftCell && jrightCell > -1) 
			axb1((g->ff[jface]).rmat,&(s->ddq[NVAR*jrightCell]),
			       g->F[m],-dtfac,NVAR);
		      else if (icell==jrightCell)
			axb1((g->ff[jface]).lmat,&(s->ddq[NVAR*jleftCell]),
			     g->F[m],dtfac,NVAR);
		    }
		}
	      m++;
	    }
	  //
	  // invert using appropriate banded block solver
	  //
/* 	  if (i==0) { */
/* 	  for(m=0;m<chainSize;m++) */
/* 	    for(j=0;j<NQ;j++) */
/* 	      { */
/* 		for(k=0;k<NQ;k++) */
/* 		  { */
/* 		    fprintf(fp1,"%.16e ",g->A[m][j][k]); */
/* 		    fprintf(fp2,"%.16e ",g->B[m][j][k]); */
/* 		    fprintf(fp3,"%.16e ",g->C[m][j][k]); */
/* 		  } */

/* 		//g->F[m][j]=1; */
/* 		fprintf(fp4,"%.16e\n",g->F[m][j]); */
/* 		fprintf(fp1,"\n"); */
/* 		fprintf(fp2,"\n"); */
/* 		fprintf(fp3,"\n"); */
/* 	      } */
/* 	  fclose(fp1); */
/* 	  fclose(fp2); */
/* 	  fclose(fp3); */
/* 	  fclose(fp4); */
/* 	  } */
	  if (idv==1) blockTridagPeriodic(g->A,g->B,g->C,g->F,chainSize,NQ);
	  if (idv==0) blockTridag(g->A,g->B,g->C,g->F,chainSize,NQ);

/* 	  if (i==0) { */
/* 	  for(m=0;m<chainSize;m++) */
/* 	    for(j=0;j<NQ;j++) */
/* 		fprintf(fp5,"%.16e\n",g->F[m][j]); */
/* 	  fclose(fp5); */
/* 	  } */
	  //
	  // set ddq can be to be updated by linear RHS if called
          //
	  m=0;
	  for(f=f1;f<f2-1;f++)
	    {
	      //
	      iface=g->chainConn[f];
	      icell=g->faces[8*iface+4];
	      for(j=0;j<NQ;j++) s->ddq[NVAR*icell+j]=g->F[m][j];
	      m++;
	    }
	}
      for(i=0;i<ntotal;i++) s->dq[i]+=s->ddq[i];
      l2rho=0.;
      icell=0;
      for(i=0;i<g->ncells;i++)
	{
	  if ((l2rho) < fabs(s->ddq[5*i])) 
	  {
	    icell=i;
	    l2rho=fabs(s->ddq[5*i]);
	  }
	}
      //trace(icell);
      //tracef(l2rho);
      if (g->msweep > 1) {
       computeLinearRHS(g,s,cflnum,&linearl2rho,&linearlinfrho);
       tracef(linearl2rho);
       }
      //if (g->msweep > 1) tracef(linearl2rho);
    }

  l2rhodq=0.;
  for(i=0;i<g->ncells;i++)
    {
      if ((l2rhodq) < fabs(s->dq[5*i])) 
	{
	  icell=i;
	  l2rhodq=fabs(s->dq[5*i]);
	}
    }
  //trace(icell);
  //tracef(l2rhodq);

  //
  // update q
  //
  m=0;
  for(i=0;i<g->ncells;i++)
    for(j=0;j<NVAR;j++)
      {
	s->q[m]+=s->dq[m];
	m++;
      }
}
      
  
