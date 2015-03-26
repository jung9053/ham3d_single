// ##################################################################
//
// computeForce.c
//
// calculate lift and drag(pressure + skin friction)
//
// Written by Dr. Jayanarayanan Sitaraman
// Modified by Yong Su Jung
// ##################################################################

#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void computeForce(GRID *g,SOLN *s)
{
  int i,n;
  double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
  double x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8;
  double x0,xp,y0,yp,z0,zp;
  double dis,vel,tau;
  double dudx,dvdy,dwdz,dudy,dvdx,dwdx,dudz,dvdz,dwdy;
  double txx,tyy,tzz,txy,txz,tyz;
  double tracx,tracy,tracz;
  double temp,mu;
  double rgas=1./1.4, c2b=0.3678;


  int iface,node1,node2,node3,node4,icell;
  int node5,node6,node7,node8;
  
  double rho,rhou,rhov,rhow,e,p;
  double sref,fac;
  double cs,ss;

  sref = 3.141592/4.0; //Sphere
  s->fx=s->fy=s->fz=s->cl=s->cd=0;

  for (i=0;i<g->nbfaces;i++)
    {
      iface = g->bfaces[i];
      node1 = g->faces[8*iface];
      node2 = g->faces[8*iface+1];
      node3 = g->faces[8*iface+2];
      node4 = g->faces[8*iface+3];       
       
      x1    = g->x[3*node1];
      y1    = g->x[3*node1+1];
      z1    = g->x[3*node1+2];

      x2    = g->x[3*node2];
      y2    = g->x[3*node2+1];
      z2    = g->x[3*node2+2];

      x3    = g->x[3*node3];
      y3    = g->x[3*node3+1];
      z3    = g->x[3*node3+2];

      x4    = g->x[3*node4];
      y4    = g->x[3*node4+1];
      z4    = g->x[3*node4+2];

      icell = g->faces[8*iface+4]; //left cell
      //
      rho   = s->q[NVAR*icell];
      rhou  = s->q[NVAR*icell+1];
      rhov  = s->q[NVAR*icell+2];
      rhow  = s->q[NVAR*icell+3];
      e     = s->q[NVAR*icell+4];
      p     = (gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
      
      s->fx += p*g->cx[i]; 
      s->fy += p*g->cy[i]; 
      s->fz += p*g->cz[i]; 

      if(g->visc)
      {
         //find avg. node points
         //node5 = g->conn[8*icell+1]-1;
         //node6 = g->conn[8*icell+2]-1;
         //node7 = g->conn[8*icell+5]-1;
         //node8 = g->conn[8*icell+6]-1;

         //x5    = g->x[3*node5];
         //y5    = g->x[3*node5+1];
         //z5    = g->x[3*node5+2];

         //x6    = g->x[3*node6];
         //y6    = g->x[3*node6+1];
         //z6    = g->x[3*node6+2];

         //x7    = g->x[3*node7];
         //y7    = g->x[3*node7+1];
         //z7    = g->x[3*node7+2];

         //x8    = g->x[3*node8];
         //y8    = g->x[3*node8+1];
         //z8    = g->x[3*node8+2];

         //x0 = 0.125*(x1+x2+x3+x4+x5+x6+x7+x8);
         //y0 = 0.125*(y1+y2+y3+y4+y5+y6+y7+y8);
         //z0 = 0.125*(z1+z2+z3+z4+z5+z6+z7+z8);

         //xp = 0.25*(x1+x2+x3+x4);
         //yp = 0.25*(y1+y2+y3+y4);
         //zp = 0.25*(z1+z2+z3+z4);

         //dis = sqrt((xp-x0)*(xp-x0)+(yp-y0)*(yp-y0)+(zp-z0)*(zp-z0)); 

         //printf("icell:%d: %d %d %d %d %d %d %d %d\n",icell,node1,node2,node3,node4,node5,node6,node7,node8);
         //printf("distance:%f\n",dis);
         dis = 0.005;
         temp = p/rho/rgas;
         mu = (c2b+1.)*temp*sqrt(temp)/(c2b+temp);
         mu = 1./(s->rey)*mu;
         //
         vel = sqrt(pow((rhou/rho),2)+pow((rhov/rho),2)+pow((rhow/rho),2));
         tau = mu*vel/(dis/g->vol[icell]);
         s->fx += tau*g->cx[i];
         s->fy += tau*g->cy[i];
         s->fz += tau*g->cz[i];
      }
    }

  fac   = 0.5*rinf*s->mach*s->mach*sref;
  s->fx= s->fx/=fac;
  s->fy= s->fy/=fac;
  s->fz= s->fz/=fac;

  cs    = cos(s->alpha*deg2rad);
  ss    = sin(s->alpha*deg2rad);

  if(g->test==2) //3d wing
  {
   s->cl = s->fy*cs-s->fx*ss;
   s->cd = s->fy*ss+s->fx*cs;      
  }
  else //general (sphere or ROBIN)
  {
   s->cl = s->fz*cs-s->fx*ss;
   s->cd = s->fz*ss+s->fx*cs;
  }

}


