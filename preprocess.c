// ##################################################################
//
// preprocess.c
//
// Preprocessing code
// Written by Dr. Jayanarayanan Sitaraman
// Modified by Yong Su Jung
// Note: When accessing elements of 'x', the access is random.
//       Efforts can be focussed towards a more linear accessing
//       of data to better utilize cache and increase efficiency
//       of parallelization
// ##################################################################

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define nearBodyRadius 1.1

void preprocess(GRID *g)
{
  int i,m,i1,i2,i3,i4,kk;
  int leftCell,rightCell;
  int leftFaceIndx,rightFaceIndx;
  int f,f1,f2;
  int icell,iface;
  double r1,r2,r3,r4,dis,denom;
  double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,dvedge;
  FILE *fp, *fp1, *fp2, *fp3;
  double volmin,volmax;
  double r_mid[2],s_vec[2];
  double xa,xb,ya,yb,za,zb;
  
  int j,k,ff1,ff2,nbase,nbloop,tmp;
  int nsidefaces, nbases;

  int neigb[10000];

  // use periodic bc condition  
  if(g->test==1) periodic_bc(&g[0]); 
  // use surface boundary condition
  if(g->test==0) g->nbfaces = g->ncells/(g->nstrand-1);
  // use infinite wing case 
  if(g->test==2) 
  {
     nbases = g->ncells/(g->nstrand-1); // face on one surface layer
     // nsidefaces = the number of faces only side faces
     nsidefaces = g->nfaces - nbases*g->nstrand;

     g->nbfaces=0;
     for (i=0;i<g->nfaces;i++)
     {
       if(g->faces[8*i+6] == -1 && i < nsidefaces)
       {
         i1=g->faces[8*i];
         i2=g->faces[8*i+1];
         i3=g->faces[8*i+2];
         i4=g->faces[8*i+3];
         //
         x1=g->x[3*i1];
         y1=g->x[3*i1+1];

         x2=g->x[3*i2];
         y2=g->x[3*i2+1];

         x3=g->x[3*i3];
         y3=g->x[3*i3+1];

         x4=g->x[3*i4];
         y4=g->x[3*i4+1];

         r_mid[0] = (x1+x2+x3+x4)/4.;
         r_mid[1] = (y1+y2+y3+y4)/4.;

         dis = sqrt(r_mid[0]*r_mid[0]+r_mid[1]*r_mid[1]);
         if(dis<nearBodyRadius) g->nbfaces++;
        } //if
     } //nbface
  } //test=2

  //
  if(g->test!=1) g->bfaces=(int *) malloc(sizeof(int)*g->nbfaces);
  g->vol=(double *) malloc(sizeof(double)*g->ncells);
  g->cx=(double *) malloc(sizeof(double)*g->nbfaces);
  g->cy=(double *) malloc(sizeof(double)*g->nbfaces);
  g->cz=(double *) malloc(sizeof(double)*g->nbfaces);
  
  g->ncx=(double *) malloc(sizeof(double)*g->nbfaces);
  g->ncy=(double *) malloc(sizeof(double)*g->nbfaces);
  g->ncz=(double *) malloc(sizeof(double)*g->nbfaces);
  
  g->neig=(int *) malloc(6*sizeof(int)*g->ncells);
  g->c2f=(int *) malloc(6*sizeof(int)*g->ncells);
  g->c2chain=(int *) malloc(3*sizeof(int)*g->ncells);  
  //
  for(i=0;i<g->ncells;i++) g->vol[i]=0.;
  //
  m=0;
  for(i=0;i<g->nfaces;i++)
  {
      leftCell      = g->faces[8*i+4];
      rightCell     = g->faces[8*i+6];
      leftFaceIndx  = g->faces[8*i+5];
      rightFaceIndx = g->faces[8*i+7];
        
      g->c2f[6*leftCell+leftFaceIndx]=i;
      if (rightCell > -1) 
      {
	      g->c2f[6*rightCell+rightFaceIndx]=i;
      }

      //
      // collect neighbors
      //
      g->neig[6*leftCell+leftFaceIndx]=rightCell;
      if (rightCell > -1) 
      {
	     g->neig[6*rightCell+rightFaceIndx]=leftCell;
      }

      //
      // indices of face nodes
      //
      i1=g->faces[8*i];
      i2=g->faces[8*i+1];
      i3=g->faces[8*i+2];
      i4=g->faces[8*i+3];

      //
      x1=g->x[3*i1];
      y1=g->x[3*i1+1];
      z1=g->x[3*i1+2];

      x2=g->x[3*i2];
      y2=g->x[3*i2+1];
      z2=g->x[3*i2+2];

      x3=g->x[3*i3];
      y3=g->x[3*i3+1];
      z3=g->x[3*i3+2];

      x4=g->x[3*i4];
      y4=g->x[3*i4+1];
      z4=g->x[3*i4+2];


      r_mid[0] = (x1+x2+x3+x4)/4.;
      r_mid[1] = (y1+y2+y3+y4)/4.;
      r_mid[2] = (z1+z2+z3+z4)/4.;

      xa = x3-x1;
      xb = x2-x4;
      ya = y3-y1;
      yb = y2-y4;
      za = z3-z1;
      zb = z2-z4;

      s_vec[0]=0.5*(za*yb-ya*zb);
      s_vec[1]=0.5*(xa*zb-za*xb);
      s_vec[2]=0.5*(ya*xb-xa*yb);
     

      dvedge=1./3.*(r_mid[0]*s_vec[0]+r_mid[1]*s_vec[1]+r_mid[2]*s_vec[2]);

      // 
      // compute cell volumes
      // hexahedron volume by formula of Gauss
      //
      g->vol[leftCell]+=dvedge;
      
      if (rightCell > -1) g->vol[rightCell]-=dvedge;

      if(g->test==0) //Sphere
      {
        if(g->faces[8*i+6]==-2)
        {
          g->bfaces[m]    = i;
          g->faces[8*i+6] = -(g->visc+2);
          icell           = g->faces[8*m+4];       // left cell indx
          iface           = g->faces[8*m+5];       // left face 
          g->neig[6*icell+iface] = -(g->visc+2);   //
  
          // face normal vector(x,y,z)
          g->cx[m] = s_vec[0]; 
          g->cy[m] = s_vec[1];
          g->cz[m] = s_vec[2];
          
          denom = sqrt(pow(s_vec[0],2)+pow(s_vec[1],2)+pow(s_vec[2],2)); 
          g->ncx[m] = s_vec[0]/denom; 
          g->ncy[m] = s_vec[1]/denom;
          g->ncz[m] = s_vec[2]/denom;
          m++;
        }
      }

      if(g->test==2) //infinite wing case
      { // if a boundary face
        if(g->faces[8*i+6] == -1 && i < nsidefaces)
        {
         dis = sqrt(r_mid[0]*r_mid[0]+r_mid[1]*r_mid[1]);
         if (dis<nearBodyRadius)
          {
            g->bfaces[m]     = i;
            // face normal vector(x,y,z)
            g->cx[m] = s_vec[0]; 
            g->cy[m] = s_vec[1];
            g->cz[m] = s_vec[2];
            
            m++;
            g->faces[8*i+6] = -(g->visc+2);
            icell           = g->faces[8*i+4];       // left cell indx
            iface           = g->faces[8*i+5];       // left face 
            g->neig[6*icell+iface] = -(g->visc+2);   //??
          }
        }
      }
    }//nface
      
    //Initial value for volmin and volmax
    volmin=1E15;
    volmax=0.;
  //find maximum and minimum v(lume
  for(i=0;i<g->ncells;i++)
  {
   volmin=(volmin < g->vol[i]) ? volmin : g->vol[i]; 
   volmax=(volmax > g->vol[i]) ? volmax : g->vol[i];
  }
  tracef(volmin);
  tracef(volmax);
  tracef(g->vol[0]);
  
  //=============================================================
  // find cell to chain connectivity
  //============================================================
  for(i=0;i<g->ncells;i++)
  { 
    g->c2chain[3*i]=-1;
    g->c2chain[3*i+1]=-1;
    g->c2chain[3*i+2]=-1;
  }
  for(i=0;i<g->nchains;i++)
  {
    f1=g->faceStartPerChain[i];
    f2=g->faceStartPerChain[i+1];

    for(f=f1;f<f2-1;f++)
	 {
	   iface=g->chainConn[f]; //iface = face index that consist chain
	   leftCell=g->faces[8*iface+4]; // left cell

	   if (g->c2chain[3*leftCell] > -1) 
	   {
          if (g->c2chain[3*leftCell+1] >-1) 
          { 
	       g->c2chain[3*leftCell+2]=i;
          }
          else
          {
             g->c2chain[3*leftCell+1]=i;
          }
	   }
	   else
	   {
	       g->c2chain[3*leftCell]=i;
	   }
	 }
  }

   // WRITE TO VERIFY (Solid wall face)
   fp = fopen("./output/bface.dat","w");
   for (i=0;i<g->nbfaces;i++)
   {
     iface = g->bfaces[i];
     i1=g->faces[8*iface];
     i2=g->faces[8*iface+1];
     i3=g->faces[8*iface+2];
     i4=g->faces[8*iface+3];
     fprintf(fp," %f %f %f\n",g->x[3*i1],g->x[3*i1+1],g->x[3*i1+2]);
     fprintf(fp," %f %f %f\n",g->x[3*i2],g->x[3*i2+1],g->x[3*i2+2]);
     fprintf(fp," %f %f %f\n",g->x[3*i3],g->x[3*i3+1],g->x[3*i3+2]);
     fprintf(fp," %f %f %f\n",g->x[3*i4],g->x[3*i4+1],g->x[3*i4+2]);
   }
   fclose(fp);
}
		

