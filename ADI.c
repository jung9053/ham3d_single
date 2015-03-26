// ##################################################################
//
// ADI.c
//
// Alternating Direction Implicit method
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

void ADI(GRID *g,SOLN *s,double cflnum,double dt)
{
  //
  int i,j,k,l,m,f,n;
  int mm1,f1,f2,iface,idv,chainSize;
  double ds[3];
  double lmat[NQ][NQ];
  double rmat[NQ][NQ];
  double dsnorm;
  int node1,node2,node3,node4,leftCell,rightCell,icell;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
  double xa,ya,za,xb,yb,zb;
  double r01,r02,r03,a1,a2,a3,b1,b2,b3,c1,c2,c3,pp,dtfac;
  double linearl2rho,linearlinfrho;
  int isweep,ntotal;
  double ref[3][3];

  ntotal = g->ncells*NVAR;   
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
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
  for(isweep=0;isweep<g->msweep;isweep++)
  {
   for(i=0;i<g->nchains;i++)
	{
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
	  idv=(g->chainConn[f1]==g->chainConn[f2-1]);
	  m=0;
	  chainSize=(f2-idv-f1);
	  
	  //
	  // make matrices A,B,C and vector F for
	  // inversion.
	  //
	  for(f=f1;f<f2-idv;f++)
	  {
	    iface=g->chainConn[f];
	    leftCell=g->faces[8*iface+4];
	    rightCell=g->faces[8*iface+6];
	    // 
	    // construct left and right matrices for this face
	    //
	   for(j=0;j<NQ;j++)
		for(k=0;k<NQ;k++)
		 {
		   lmat[j][k] = (g->ff[iface]).lmat[j][k];
		   rmat[j][k] = (g->ff[iface]).rmat[j][k];
		 }

	      if (idv==1) 
	      // close loop
	      {
		   mm1=(m==0)?chainSize-1:m-1;//mm1 = mm - 1
		                              // if mm=0, mm1 = nmax (periodic)
		   for(j=0;j<NQ;j++)
		     {
		       for(k=0;k<NQ;k++)
		       {
		   	   g->B[m][j][k]+=(lmat[j][k]);
		   	   g->B[mm1][j][k]-=(rmat[j][k]);
		   	   g->A[m][j][k]+=(rmat[j][k]);
		   	   g->C[mm1][j][k]-=(lmat[j][k]);
		       }
		       g->F[m][j]=s->r[NVAR*leftCell+j];
		     }
	      }
	      else
	      // open loop
		   {
		   mm1=(m-1);
		   if (rightCell==-2)
		   {

		      node1=g->faces[8*iface];
		      node2=g->faces[8*iface+1];
		      node3=g->faces[8*iface+2];
		      node4=g->faces[8*iface+3];

            x1=g->x[3*node1];
		      y1=g->x[3*node1+1];
		      z1=g->x[3*node1+2];

		      x2=g->x[3*node2];
		      y2=g->x[3*node2+1];
		      z2=g->x[3*node2+2];

            x3=g->x[3*node3];
		      y3=g->x[3*node3+1];
		      z3=g->x[3*node3+2];

		      x4=g->x[3*node4];
		      y4=g->x[3*node4+1];
		      z4=g->x[3*node4+2];

            // 3D face normal vector (direction?)
            xa = x1 - x3; xb = x2 - x4;
            ya = y1 - y3; yb = y2 - y4;
            za = z1 - z3; zb = z2 - z4;
    
            ds[0] = 0.5*(za*yb - ya*zb);
            ds[1] = 0.5*(xa*zb - za*xb);
            ds[2] = 0.5*(ya*xb - xa*yb);

            dsnorm=ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2];


            //make reflection matrix
            ref[0][0] = (-ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2])/dsnorm;
            ref[0][1] = -2.*ds[0]*ds[1]/dsnorm;
            ref[0][2] = -2.*ds[0]*ds[2]/dsnorm;
            ref[1][1] = (-ds[1]*ds[1]+ds[0]*ds[0]+ds[2]*ds[2])/dsnorm;
            ref[1][0] = -2.*ds[1]*ds[0]/dsnorm;
            ref[1][2] = -2.*ds[1]*ds[2]/dsnorm;
            ref[2][2] = (-ds[2]*ds[2]+ds[0]*ds[0]+ds[1]*ds[1])/dsnorm;
            ref[2][0] = -2.*ds[2]*ds[0]/dsnorm;
            ref[2][1] = -2.*ds[2]*ds[1]/dsnorm;

            a1 = ref[0][0]; 
            b1 = ref[0][1];
            c1 = ref[0][2];
            a2 = ref[1][0];
            b2 = ref[1][1];
            c2 = ref[1][2];
            a3 = ref[2][0];
            b3 = ref[2][1];
            c3 = ref[2][2];

            
		     for (j=0;j<NQ;j++)
			  {
			    r01=rmat[j][1];
			    r02=rmat[j][2];
			    r03=rmat[j][3];

			    rmat[j][1]=a1*r01+b1*r02+c1*r03;
			    rmat[j][2]=a2*r01+b2*r02+c2*r03;
			    rmat[j][3]=a3*r01+b3*r02+c3*r03;
			  }
		   }
		  //
		  if (rightCell < 0 && idv==0 && f==f2-1) m--;

		  //
		  for(j=0;j<NQ;j++)
		  {
		     for(k=0;k<NQ;k++)
			  {
			      g->B[m][j][k]+=(lmat[j][k]);
			      if (mm1 > -1 && rightCell > -1)
			      {
			        g->A[m][j][k]+=(rmat[j][k]);
			        g->B[mm1][j][k]-=(rmat[j][k]);
			        g->C[mm1][j][k]-=(lmat[j][k]);
			      }
			      else
			      {
			        if (rightCell==-2) g->B[m][j][k]+=(rmat[j][k]);
			      }
			  }
		     g->F[m][j]=s->r[NVAR*leftCell+j];
		  }
		 } // idv	  
	    m++;
	  } //f = f1 ~ f2-1 
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
		    g->B[m][j][k]*=dtfac;
		    g->C[m][j][k]*=dtfac;
		  }
	      for(j=0;j<NQ;j++)
		   g->B[m][j][j]+=1.0;
	      m++;
	    }
	  //
	  // invert using appropriate banded block solver
	  //
     
	  if (idv==1) blockTridagPeriodic(g->A,g->B,g->C,g->F,chainSize,NQ);
	  if (idv==0) blockTridag(g->A,g->B,g->C,g->F,chainSize,NQ);
        
	  //
	  // reassign values back at the unknown locations
	  //
	  m=0;
	  for(f=f1;f<f2-1;f++)
	    {
	      iface=g->chainConn[f];
	      leftCell=g->faces[8*iface+4];
	      for(j=0;j<NQ;j++)
		   s->r[NVAR*leftCell+j]=g->F[m][j];
	      m++;
	    }      
	}

   //for(i=0;i<g->ncells*NVAR;i++) s->dq[i]+=s->r[i];
   for(i=0;i<ntotal;i++) s->dq[i]+=s->r[i];
   computeLinearRHS(g,s,cflnum,&linearl2rho,&linearlinfrho);

   tracef(linearl2rho);
   //tracef(s->dq[0]);
  }
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
