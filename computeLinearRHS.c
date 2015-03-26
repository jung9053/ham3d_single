#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"


void computeLinearRHS(GRID *g,SOLN *s,double cflnum,double *l2rho,double *linfrho)
{
  //
  int i,j,k,m,f,n;
  int f1,f2;
  int is,ie;
  int iface;
  int idv;
  int chainSize,ntotal;
  double ds[3];
  double leftState[NVAR];
  double rightState[NVAR];
  double dleftState[NVAR];
  double drightState[NVAR];
  double consVar[NVAR];
  double dqvar[NVAR];
  double flux[NVAR];
  double gm1=gamm-1.0;
  double gamma1=gamm;
  double specRadius;
  double faceVel=0.;
  double dsnorm;
  double dtfac;
  int node1,node2,node3,node4,leftCell,rightCell,icell,iflag,nbase;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;  
  double xa,ya,za,xb,yb,zb;
  double pp;
  double th,qt,eps;
  double dscheck[2];
  double ref[3][3];

  int nghost,order;

  order = g->order;
               nghost = 2;
  if(order==5) nghost = 3;



  nbase = (g->nfaces-(g->ncells)/(g->nstrand-1)*g->nstrand)/(g->nstrand-1);
  //
  // add diagonal term to the linear residual
  //
  ntotal = g->ncells*NVAR;
  for(i=0;i<ntotal;i++) s->r[i]=(s->r0[i]-s->dq[i]);
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  for(i=0;i<g->nchains;i++)
    {
      iflag = 0;
      f1=g->faceStartPerChain[i];
      f2=g->faceStartPerChain[i+1];
      m=nghost;
     for(f=f1;f<f2;f++)
	  {
	    iface=g->chainConn[f];
	    g->cindx[m]=g->faces[8*iface+4];
	    m++;
	  }
      //
      // add buffer cells to the chain
      //
      if (g->chainConn[f1]==g->chainConn[f2-1])
	  {
	  //
	  // this is a closed chain
	  // make it periodic
	  //
	  iflag         =  0;
	  f             =  f1+1;
	  iface         =  g->chainConn[f];
	  g->cindx[m]   =  g->faces[8*iface+4];
	  m++;
	  chainSize     =  m;
	  m             =  0;
	  for(f=f2-nghost-1;f<f2-1;f++)
	    {
	      iface=g->chainConn[f];
	      g->cindx[m]=g->faces[8*iface+4];
	      m++;
	    }
	}
      else
	{
	  //
	  // this is a open chain
	  // -ve index indicates necessity to create
	  // ghost cells
	  //
	  iflag = 1;
	  
	  // solid bc 
     if(g->test != 1)
     {
        if(order==5) 
        {
          m--;
          g->cindx[m] = -g->cindx[m];
          m++;
          g->cindx[m] = -g->cindx[m-3];
          m++;
          g->cindx[m] = -g->cindx[m-5];
          chainSize = m+1;
          m = 0;
          g->cindx[m] = -g->cindx[m+5];
          m = 1;
          g->cindx[m] = -g->cindx[m+3];
          m = 2;
          g->cindx[m] = -g->cindx[m+1];
        }
        else
        {
	       m--;
	       g->cindx[m]=-g->cindx[m];
	       m++;
	       g->cindx[m]=-g->cindx[m-3];
	       chainSize=m+1;
	       m=0;
	       g->cindx[m]=-g->cindx[m+3];
	       m=1;
	       g->cindx[m]=-g->cindx[m+1];
        }

	          
       // prriodic bc at only for strand grid 
       if(g->test ==2){
         iface = g->chainConn[f1];
         
         if(iface>nbase*(g->nstrand-1)-1)
         {
           if(order==5)
           {
             m = chainSize-1;
             g->cindx[2]   = g->cindx[g->nstrand];
             g->cindx[1]   = g->cindx[g->nstrand-1];
             g->cindx[0]   = g->cindx[g->nstrand-2];
             g->cindx[m-2] = g->cindx[3];
             g->cindx[m-1] = g->cindx[4]; 
             g->cindx[m]   = g->cindx[5]; 
           }
           else
           {
             m = chainSize-1;
             g->cindx[1]   = g->cindx[g->nstrand];
             g->cindx[0]   = g->cindx[g->nstrand-1];
             g->cindx[m-1] = g->cindx[2]; 
             g->cindx[m]   = g->cindx[3];
           }
         }
       }
	  }
     if(g->test==1)
     {
       apply_periodic(&g[0],f1,f2,m);
       chainSize=m+1;
       if(order==5) chainSize = m+2;
     }
	}
//============ print for verification===========	
    //if(iflag==1) //open chain:1, close :0
    //{
    //printf("i:%d Stopping code\n",i);
    //printf("chain size=%d\n",chainSize);

    //for(k=0;k<=chainSize-1;k++) 
    //{
    //printf("cell index:%d\n",g->cindx[k]);  
    //}
    ////exit(1); 
    //}
//============================================

	//======================================
   for(j=0;j<chainSize;j++)
	{
	  icell=g->cindx[j];
	  if (icell >=0 ||(icell==0&&j==nghost)||(icell==0&&j==chainSize-nghost-1)) 
	  //if (icell >=0)
	  {
	   m=NVAR*icell;
	   for(k=0;k<NVAR;k++)
		{
		  consVar[k]=s->q[m];
		  dqvar[k]=s->dq[m];
        m++;
		}
	     // collect primitive variables in 
	     g->f[j][0]=consVar[0];
	     g->f[j][1]=consVar[1];
	     g->f[j][2]=consVar[2];
	     g->f[j][3]=consVar[3];
	     g->f[j][4]=consVar[4];

	     // 
	     g->df[j][0]=dqvar[0];
	     g->df[j][1]=dqvar[1];
	     g->df[j][2]=dqvar[2];
	     g->df[j][3]=dqvar[3];
	     g->df[j][4]=dqvar[4];
	  }
	  //else
     if(icell<0||(icell==0&&j==nghost-1)||(icell==0&&j==chainSize-nghost))

	  {
	  //
	  // do ghost cells
	  // based on whether they are on the solid boundary on that
	   if (j < nghost) 
		{
		  iface=g->chainConn[f1];
		}

	   else 
		{
		iface=g->chainConn[f2-1];
		}

	   rightCell=g->faces[8*iface+6];
	      
	   if (rightCell==-2)  /* this is a face on solid wall */
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
		  
        if(icell>0){printf("Stopping! negative icell at i: %d\n",i);exit(1);}
		  icell=-icell;
		  //printf("icell:%d\n",icell);
		  m=NVAR*icell;
		  for(k=0;k<NVAR;k++)
        {
		    consVar[k]=s->q[m];
          dqvar[k]=s->dq[m];
          m++;
        }
		  dsnorm=ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2];
        //make reflection matrix
        //ref[0][0] = 1.-2.*ds[0]*ds[0]/dsnorm;
        //ref[0][1] = -2.*ds[0]*ds[1]/dsnorm;
        //ref[0][2] = -2.*ds[0]*ds[2]/dsnorm;
        //ref[1][1] = 1.-2.*ds[1]*ds[1]/dsnorm;
        //ref[1][0] = -2.*ds[1]*ds[0]/dsnorm;
        //ref[1][2] = -2.*ds[1]*ds[2]/dsnorm;
        //ref[2][2] = 1.-2.*ds[2]*ds[2]/dsnorm;
        //ref[2][0] = -2.*ds[2]*ds[0]/dsnorm;
        //ref[2][1] = -2.*ds[2]*ds[1]/dsnorm;

        ref[0][0] = (-ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2])/dsnorm;
        ref[0][1] = -2.*ds[0]*ds[1]/dsnorm;
        ref[0][2] = -2.*ds[0]*ds[2]/dsnorm;
        ref[1][1] = (-ds[1]*ds[1]+ds[0]*ds[0]+ds[2]*ds[2])/dsnorm;
        ref[1][0] = -2.*ds[1]*ds[0]/dsnorm;
        ref[1][2] = -2.*ds[1]*ds[2]/dsnorm;
        ref[2][2] = (-ds[2]*ds[2]+ds[0]*ds[0]+ds[1]*ds[1])/dsnorm;
        ref[2][0] = -2.*ds[2]*ds[0]/dsnorm;
        ref[2][1] = -2.*ds[2]*ds[1]/dsnorm;

        //
		  //=======================================================
		  g->f[j][0]  =  consVar[0];
		  g->f[j][1]  =  (consVar[1]*ref[0][0]+consVar[2]*ref[0][1]+consVar[3]*ref[0][2]);
		  g->f[j][2]  =  (consVar[1]*ref[1][0]+consVar[2]*ref[1][1]+consVar[3]*ref[1][2]);
		  g->f[j][3]  =  (consVar[1]*ref[2][0]+consVar[2]*ref[2][1]+consVar[3]*ref[2][2]);
		  g->f[j][4]  =  consVar[4];		  	

		  g->df[j][0]  =  dqvar[0];
		  g->df[j][1]  =  (dqvar[1]*ref[0][0]+dqvar[2]*ref[0][1]+dqvar[3]*ref[0][2]);
		  g->df[j][2]  =  (dqvar[1]*ref[1][0]+dqvar[2]*ref[1][1]+dqvar[3]*ref[1][2]);
		  g->df[j][3]  =  (dqvar[1]*ref[2][0]+dqvar[2]*ref[2][1]+dqvar[3]*ref[2][2]);
		  g->df[j][4]  =  dqvar[4];		  	
        //========================================================

		}
	   else //this is for far field bc. 
		{
		  if(g->test==1) 
        {
          printf("Periodic bc has a problem!\n");
          exit(1);
        } 
		  g->f[j][0]  =  rinf;
		  g->f[j][1]  =  rinf*s->uinf;
		  g->f[j][2]  =  rinf*s->vinf;
        g->f[j][3]  =  rinf*s->winf;
		  g->f[j][4]  =  pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf+s->winf*s->winf);
 
		  g->df[j][0]=0;
		  g->df[j][1]=0;
		  g->df[j][2]=0;
		  g->df[j][3]=0;
		  g->df[j][4]=0;
		 }
	
	  }

	}
      
      is=nghost-1;
      ie=chainSize-1;
      th=1./3;
      qt=0.25;
      if (g->order==1) qt=0.0;
      eps=1e-10;

      if(order==1 || order==3)
      {  
      muscld_deriv(g->f,g->dql,g->dqr,g->f2,g->df,is,ie,th,qt,eps,chainSize,NVAR);
      }
    
      if(order==5) weno_deriv(g->f,g->dql,g->dqr,g->df,is,ie,eps,chainSize,NVAR); 

      n=is;
      idv=(g->chainConn[f1]==g->chainConn[f2-1]);

      for(f=f1;f<f2-idv;f++)
	   {
	   iface=g->chainConn[f];
	   leftCell=g->faces[8*iface+4];
	   rightCell=g->faces[8*iface+6];

	   for(m=0;m<NVAR;m++)
	   {
	    if (f==f2-idv-1 && idv==0) 
		 {
		  dleftState[m]=g->dql[n][m];
		  drightState[m]=g->dqr[n+1][m];
		 }
	    else
		 {
		  dleftState[m]=g->dqr[n+1][m];
		  drightState[m]=g->dql[n][m];
		 }
	   }
	  
	  for(j=0;j<NVAR;j++)
	  {
	    flux[j]=0; 
		 for(k=0;k<NVAR;k++)
       {
		  flux[j]+=(((g->ff[iface]).lmat[j][k]*dleftState[k])+
			    ((g->ff[iface]).rmat[j][k]*drightState[k]));
	    }
	  }
	  //
	  m=NVAR*leftCell;
	  dtfac=cflnum/s->sigma[leftCell];
	  for(j=0;j<NVAR;j++)
	    {
	      s->r[m]-=(flux[j]*dtfac);
	      m++;
	    }
	  if (rightCell > -1) 
	    {
	      m=NVAR*rightCell;
	      dtfac=cflnum/s->sigma[rightCell];
	      for(j=0;j<NVAR;j++)
         {
		     s->r[m]+=(flux[j]*dtfac);
           m++;
         }
	    } 
	  n++;
	  }
    }

  *linfrho = 0.;
  *l2rho   = 0.;
  for(i=0;i<g->ncells;i++)
    {
      if((*linfrho) < fabs(s->r[5*i]))
      {
       icell=i;
       *linfrho=fabs(s->r[5*i]);
      }
      *l2rho = *l2rho + (s->r[5*i])*(s->r[5*i]);
    }
    *l2rho = sqrt(*l2rho/(g->ncells));
}
      
