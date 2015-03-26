#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void computeRHSkv(GRID *g,SOLN *s,double *l2rho, double *linfrho)
{
  int i,j,k,m,f,n,kv;
  int f1,f2;
  int is,ie;
  int iface;
  int idv;
  int chainSize;
  int itype;
  int imode=1;
  double ds[3];
  double leftState[NVAR];
  double rightState[NVAR];
  double leftStatev[NVAR];
  double rightStatev[NVAR];
  double leftState0[NVAR];
  double rightState0[NVAR];
  double lmat[NVAR][NVAR];
  double rmat[NVAR][NVAR];
  double lmatv[NVAR][NVAR];
  double rmatv[NVAR][NVAR];
  double consVar[NVAR];
  double flux[NVAR];
  double fluxv[NVAR];
  double gm1=gamm-1.0;
  double gamma1=gamm;
  double specRadius;
  double faceVel=0.;
  double dsnorm,nynx,nx2ny;
  int node1,node2,node3,node4,leftCell,rightCell,icell;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;  
  double xa,ya,za,xb,yb,zb;
  double pp;
  double th,qt,eps;
  double ref[3][3];
  double dscheck[2];
  int    iflag,nbase;
  int    nghost,order;

  order = g->order;
               nghost = 2;
  if(order==5) nghost = 3;


  nbase = (g->nfaces-(g->ncells)/(g->nstrand-1)*g->nstrand)/(g->nstrand-1);
  //
  // zero out residual and spectral radii
  //
  dscheck[0]=dscheck[1]=0;
  for(i=0;i<NVAR*g->ncells;i++) s->r[i]=0.0;
  for(i=0;i<g->ncells;i++) s->sigma[i]=0.0;
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  for(i=0;i<g->nchains;i++)
  {
    f1=g->faceStartPerChain[i];
    f2=g->faceStartPerChain[i+1];
    m=nghost;
    for(f=f1;f<f2;f++)
	{
	  iface=g->chainConn[f];
	  g->cindx[m]=g->faces[8*iface+4];
	  g->ctype[m]=1;
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
	  f            = f1+1;
	  iface        = g->chainConn[f];
	  g->cindx[m]  = g->faces[8*iface+4];
	  g->ctype[m]  = 1;
	  m++;
	  chainSize    = m;
	  m=0;
	  for(f=f2-nghost-1;f<f2-1;f++)
	    {
	      iface       = g->chainConn[f];
	      g->cindx[m] = g->faces[8*iface+4];
	      g->ctype[m] = 1;
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
	  //solid bc
	  if(g->test!=1)
	  {
        if(order==5) 
        {
          m--;
          g->cindx[m] = -g->cindx[m];
          g->ctype[m]=-1;
          m++;
          g->cindx[m] = -g->cindx[m-3];
          g->ctype[m]=-1;
          m++;
          g->cindx[m] = -g->cindx[m-5];
          g->ctype[m]=-1;
          chainSize = m+1;
          m = 0;
          g->cindx[m] = -g->cindx[m+5];
          g->ctype[m]=-1;
          m = 1;
          g->cindx[m] = -g->cindx[m+3];
          g->ctype[m]=-1;
          m = 2;
          g->cindx[m] = -g->cindx[m+1];
          g->ctype[m]=-1;

        }
        else
        {
	       m--;
	       g->cindx[m]=-g->cindx[m];
	       g->ctype[m]=-1;
	       m++;
	       g->cindx[m]=-g->cindx[m-3];
	       g->ctype[m]=-1;
	       chainSize=m+1;
	       m=0;
	       g->cindx[m]=-g->cindx[m+3];
	       g->ctype[m]=-1;
	       m=1;
	       g->cindx[m]=-g->cindx[m+1];
	       g->ctype[m]=-1;
        }
     
     // periodic bc at only strand grid
     if(g->test==2){
       iface = g->chainConn[f1];
       //strand grid
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
           g->ctype[2]   = 1;
           g->ctype[1]   = 1;
           g->ctype[0]   = 1;
           g->ctype[m-2] = 1;
           g->ctype[m-1] = 1;
           g->ctype[m]   = 1;
         }
         else
         {
           m = chainSize-1;
           g->cindx[1]   = g->cindx[g->nstrand];
           g->cindx[0]   = g->cindx[g->nstrand-1];
           g->cindx[m-1] = g->cindx[2]; 
           g->cindx[m]   = g->cindx[3]; 
           g->ctype[1]   = 1;
           g->ctype[0]   = 1;
           g->ctype[m-1] = 1;
           g->ctype[m]   = 1;
         }
       }
     }
     }
	  //periodic bc
	  if(g->test==1){
     apply_periodic_LHS(&g[0],f1,f2,m);
     chainSize=m+1;
     if(order==5) chainSize = m+2;
	  }
	}
   //end

  
   for(j=0;j<chainSize;j++)
	{
	  icell=g->cindx[j];
	  itype=g->ctype[j];
	  if (itype >=0) 
	  {
	     m=NVAR*icell;
	     for(k=0;k<NVAR;k++)
		  {
		    consVar[k]=s->q[m];
		    m++;
		  }
	     g->f[j][0]=consVar[0];
	     g->f[j][1]=consVar[1];
	     g->f[j][2]=consVar[2];
	     g->f[j][3]=consVar[3];	
	     g->f[j][4]=consVar[4];	

	     g->fv[j][0]=consVar[0];
	     g->fv[j][1]=consVar[1];
	     g->fv[j][2]=consVar[2];
	     g->fv[j][4]=consVar[4];	

	  }
	  else
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

        // for check!
        if(icell>0){printf("Stopping! negative icell at i: %d\n",i);exit(1);}
		  icell=-icell;
		  //printf("icell:%d\n",icell);
		  //printf("=======\n");

        m=NVAR*icell;
		  for(k=0;k<NVAR;k++)
		  {
		    consVar[k]=s->q[m];
		    m++;
		  }

		  dsnorm=ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2];
        //make reflection matrix
        ref[0][0] = 1.-2.*ds[0]*ds[0]/dsnorm;
        ref[0][1] = -2.*ds[0]*ds[1]/dsnorm;
        ref[0][2] = -2.*ds[0]*ds[2]/dsnorm;
        ref[1][1] = 1.-2.*ds[1]*ds[1]/dsnorm;
        ref[1][0] = -2.*ds[1]*ds[0]/dsnorm;
        ref[1][2] = -2.*ds[1]*ds[2]/dsnorm;
        ref[2][2] = 1.-2.*ds[2]*ds[2]/dsnorm;
        ref[2][0] = -2.*ds[2]*ds[0]/dsnorm;
        ref[2][1] = -2.*ds[2]*ds[1]/dsnorm;

        // calculate ghost cell conserve variables
		  g->f[j][0]  =  consVar[0];
		  g->f[j][1]  =  (consVar[1]*ref[0][0]+consVar[2]*ref[0][1]+consVar[3]*ref[0][2]);
		  g->f[j][2]  =  (consVar[1]*ref[1][0]+consVar[2]*ref[1][1]+consVar[3]*ref[1][2]);
		  g->f[j][3]  =  (consVar[1]*ref[2][0]+consVar[2]*ref[2][1]+consVar[3]*ref[2][2]);
		  g->f[j][4]  =  consVar[4];		  	
        
        g->fv[j][0] =  consVar[0];
		  g->fv[j][1] =  (consVar[1]*ref[0][0]+consVar[2]*ref[0][1]+consVar[3]*ref[0][2]);
		  g->fv[j][2] =  (consVar[1]*ref[1][0]+consVar[2]*ref[1][1]+consVar[3]*ref[1][2]);
		  g->fv[j][3] =  (consVar[1]*ref[2][0]+consVar[2]*ref[2][1]+consVar[3]*ref[2][2]);
		  g->fv[j][4] =  consVar[4];		  	

        //===============================================================
		}
	   else if (rightCell==-3) 
		{
		  icell=-icell;
		  m=NVAR*icell;
		  for(k=0;k<NVAR;k++)
		  {
		    consVar[k]=s->q[m];
		    m++;
		  }
		  
		  g->f[j][0]=consVar[0];
		  g->f[j][1]=0; //-consVar[1];
		  g->f[j][2]=0; //-consVar[2];		  
		  g->f[j][3]=0; //-consVar[3];		  
		  g->f[j][4]=consVar[4]-(consVar[1]*consVar[1]+consVar[2]*consVar[2]
		     +consVar[3]*consVar[3])*0.5/consVar[0];

		  g->fv[j][0]=consVar[0];
		  g->fv[j][1]=-consVar[1];
		  g->fv[j][2]=-consVar[2];		  
		  g->fv[j][3]=-consVar[3];
		  g->fv[j][4]=consVar[4];
		}
	   else
		{
		  g->f[j][0]=rinf;
		  g->f[j][1]=rinf*s->uinf;
		  g->f[j][2]=rinf*s->vinf;
		  g->f[j][3]=rinf*s->winf;
		  g->f[j][4]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf+s->winf*s->winf);

		  g->fv[j][0]=rinf;
		  g->fv[j][1]=rinf*s->uinf;
		  g->fv[j][2]=rinf*s->vinf;
		  g->fv[j][3]=rinf*s->winf;
		  g->fv[j][4]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf+s->winf*s->winf);
		}
	 } //icell
	} //chainSize
      
      is=nghost-1;
      ie=chainSize-1;
      th=1./3;
      qt=0.25;
      if (g->order==1) qt=0.0;
      eps=1e-10;

     //reconstruction for each chain
     if(order==1 || order==3) 
     {
       muscld(g->f,g->ql,g->qr,g->f2,is,ie,th,qt,eps,chainSize,NVAR);
     }
     if(order==5) 
     {
       weno(g->f,g->ql,g->qr,is,ie,eps,chainSize,NVAR); 
     }

      n=is;
      idv=(g->chainConn[f1]==g->chainConn[f2-1]);

   for(f=f1;f<f2-idv;f++)
	{
	  iface=g->chainConn[f];
	  node1=g->faces[8*iface];
	  node2=g->faces[8*iface+1];
	  node3=g->faces[8*iface+2];
	  node4=g->faces[8*iface+3];

	  leftCell=g->faces[8*iface+4];
	  rightCell=g->faces[8*iface+6];

     //3D case
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
     xa = x3 - x1; xb = x2 - x4;
     ya = y3 - y1; yb = y2 - y4;
     za = z3 - z1; zb = z2 - z4;
     ds[0] = 0.5*(za*yb - ya*zb);
     ds[1] = 0.5*(xa*zb - za*xb);
     ds[2] = 0.5*(ya*xb - xa*yb);

	  
	  for(m=0;m<NVAR;m++)
	  {
	   if (f==f2-idv-1 && idv==0) 
		{
		  leftState[m]=g->ql[n][m];
		  rightState[m]=g->qr[n+1][m];
		  leftStatev[m]=g->f[n][m];
		  rightStatev[m]=g->f[n+1][m];
		}
	   else
		{
		  leftState[m]=g->qr[n+1][m];
		  rightState[m]=g->ql[n][m];
		  leftStatev[m]=g->f[n+1][m];
		  rightStatev[m]=g->f[n][m];
		}
	      
	   leftState0[m]=s->q[NVAR*leftCell+m]; //g->ql[j][m];	            
	   if (rightCell > -1) 
		{
		  rightState0[m]=s->q[NVAR*rightCell+m]; //g->qr[j+1][m];
		}
	  }
	  
     if (rightCell==-1 && g->test!=1) 
     {
	    rightState0[0]=rinf;
	    rightState0[1]=rinf*s->uinf;
       rightState0[2]=rinf*s->vinf;
       rightState0[3]=rinf*s->winf;
       rightState0[4]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+
       s->vinf*s->vinf+s->winf*s->winf);
     }
	  else if (rightCell==-2 && g->test!=1) 
	  {
	    dsnorm=ds[0]*ds[0]+ds[1]*ds[1];
       ref[0][0] = (-ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2])/dsnorm;
       ref[0][1] = -2.*ds[0]*ds[1]/dsnorm;
       ref[0][2] = -2.*ds[0]*ds[2]/dsnorm;
       ref[1][1] = (-ds[1]*ds[1]+ds[0]*ds[0]+ds[2]*ds[2])/dsnorm;
       ref[1][0] = -2.*ds[1]*ds[0]/dsnorm;
       ref[1][2] = -2.*ds[1]*ds[2]/dsnorm;
       ref[2][2] = (-ds[2]*ds[2]+ds[0]*ds[0]+ds[1]*ds[1])/dsnorm;
       ref[2][0] = -2.*ds[2]*ds[0]/dsnorm;
       ref[2][1] = -2.*ds[2]*ds[1]/dsnorm;

	    // which fomular is in ?....is it variables? flux?
		 rightState0[0]  =  leftState0[0];
		 rightState0[1]  =  (leftState0[1]*ref[0][0]+leftState0[2]*ref[0][1]+leftState0[3]*ref[0][2]);
		 rightState0[2]  =  (leftState0[1]*ref[1][0]+leftState0[2]*ref[1][1]+leftState0[3]*ref[1][2]);
		 rightState0[3]  =  (leftState0[1]*ref[2][0]+leftState0[2]*ref[2][1]+leftState0[3]*ref[2][2]);
		 rightState0[4]  =  leftState0[4];		  		  
	  }
	  
	  if (rightState[1]!=rightState0[1] && 0) 
	  {
	    trace(leftCell);
	    trace(rightCell);
            for(k=0;k<chainSize;k++)
              printf("%d %d %.16e\n",k,g->cindx[k],g->f[k][1]);
	    trace(n);
	    tracef(leftState[1]);
	    tracef(leftState0[1]);
	    tracef(rightState[1]);
	    tracef(rightState0[1]);
	    tracef(s->q[leftCell*NVAR+1])
	    trace(i);
       trace(n);
	    trace(f);
	    exit(0);
	  }
	  //

	  if (rightCell < -1) 
	  {
	    wallflux_(ds,leftState,flux,&specRadius,&gamma1);
	    wallfluxjacobian_(ds,leftState,lmat,rmat,&gamma1);
	  }
	  else
	  {
	    flux_roe3d_(ds,leftState,rightState,flux,&specRadius,&gamma1);
	    jac_roe_(ds,leftState,rightState,lmat,rmat,&gamma1,&imode);
	  }
	  //
	  
	  //flux_visc_2d_(ds,&(g->vol[leftCell]),leftStatev,rightStatev,fluxv,
	  //	       &(s->rey),&(s->pr),&(s->prtr),&gamma1,&(s->c2b),
	  //	       &(s->rgas));
	  jac_visc_3d_(ds,&(g->vol[leftCell]),leftStatev,rightStatev,fluxv,lmatv,rmatv,
	  	       &(s->rey),&(s->pr),&(s->prtr),&gamma1,&(s->c2b),
	  	       &(s->rgas));
     

	  for(j=0;j<NVAR;j++)
     for(k=0;k<NVAR;k++)
     {
	    (g->ff[iface]).lmat[j][k]=(lmat[k][j]-lmatv[k][j]);
	    (g->ff[iface]).rmat[j][k]=(rmat[k][j]-rmatv[k][j]);
     }
	  //
	  m=NVAR*leftCell;
	  for(j=0;j<NVAR;j++)
	  {
	    s->r[m]-=flux[j];
	    s->r[m]+=fluxv[j];
	    m++;
	  }
	  s->sigma[leftCell]+=specRadius;
	  if (rightCell > -1) 
	  {
	    m=NVAR*rightCell;
	    for(j=0;j<NVAR;j++)
       {
		   s->r[m]+=flux[j];
		   s->r[m]-=fluxv[j];
         m++;
       }
	    s->sigma[rightCell]+=specRadius;
	  } 
	  n++;
	}
  }

  *l2rho=0.;
  *linfrho=0.;
  for(i=0;i<g->ncells;i++)
  {
    if ((*linfrho) < fabs(s->r[5*i])) 
	 {
	   icell=i;
	   *linfrho=fabs(s->r[5*i]);
	 }
    *l2rho+=(s->r[5*i]*s->r[5*i]);
  }
  *l2rho/=g->ncells;
  *l2rho=sqrt(*l2rho);
  //tracef(s->sigma[0]);
  //trace(icell);
}
      
      
  


