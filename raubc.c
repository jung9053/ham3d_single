// ##################################################################
//
// raubc.c
//
// ##################################################################
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <string.h>
#include <stdio.h>


void raubc(GRID *g, SOLN *s)
{
  int i,k,ibtype;
  int rightCell,leftCell,node1,node2,node3,node4;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xa,ya,za,xb,yb,zb,dsnorm;
  double ds[3];
  
  
  ibtype = 1; // characteristic inflow/outflow

  if(ibtype==1) 
  {
    
    for (i=0;i<g->nfaces;i++)
    {
       rightCell = g->faces[8*i+6];
       leftCell  = g->faces[8*i+4];

       if(rightCell == -1)
       {

         for(k=0;k<NVAR;k++)
         {
           s->qb[k] = s->q[NVAR*leftCell+k]; 
         }

         s->qb[5] = rinf;
         s->qb[6] = s->uinf;
         s->qb[7] = s->vinf;
         s->qb[8] = s->winf;
         s->qb[9] = (gamm-1.0)*(s->einf-0.5*(s->uinf*rinf*s->uinf*rinf+s->vinf*rinf*s->vinf*rinf+s->winf*rinf*s->winf*rinf)/rinf);

		  node1=g->faces[8*i];
		  node2=g->faces[8*i+1];
        node3=g->faces[8*i+2];
        node4=g->faces[8*i+3];
		  
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

        s->qb[10] = ds[0]/sqrt(dsnorm);
        s->qb[11] = ds[1]/sqrt(dsnorm);
        s->qb[12] = ds[2]/sqrt(dsnorm);


        // using riemann invariant
        rie1d(g,s);

         // set again
         for(k=0;k<NVAR;k++)
         {
           s->q[NVAR*leftCell+k] = s->qb[k]; 
         }
       }      
    }
  }  
  else // other boundary condition.. 
  {
    printf("Need to be implemented\n");
  }




}
// ####################################
// END OF FILE
// ###################################

