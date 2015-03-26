// ##################################################################
//
// initflow.c
//
// Initialize the flow 
// Written by Dr. Jayanarayanan Sitaraman
// Modified by Yong Su Jung
// ##################################################################
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NQ 5 //  for 3D
//
void initflow(GRID *g,SOLN *s, int irest)
{
  int i,j,m,k; 
  FILE *fp;
 //
  double xc,yc,zc;
  int n1,n2,n3,n4,n5,n6,n7,n8;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
         x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;
  double vstr,xv0,yv0,zv0;
  double pi,c1,c2,xx,yy,zz,r2,rr;
  double r,p,du,dv,dw,tp,dtp,tp0;
  double gm1;
  double xv,yv,zv;

  //
  fp = fopen("input.ham3d","r");
  fscanf(fp,"Mach=%lf\n",&s->mach);
  fscanf(fp,"alpha=%lf\n",&s->alpha);
  fscanf(fp,"beta=%lf\n",&s->beta);
  fscanf(fp,"rey=%lf\n",&s->rey);
  s->rey=s->rey/s->mach;
  fclose(fp);
  tracef(s->mach);
  tracef(s->alpha);
  tracef(s->beta);
  //
  if(g->test==0 || g->test==1)
  {
    s->uinf = s->mach*cos(s->alpha*deg2rad)*cos(s->beta*deg2rad);
    s->vinf = s->mach*cos(s->alpha*deg2rad)*sin(s->beta*deg2rad);
    s->winf = s->mach*sin(s->alpha*deg2rad);
  }
  else if(g->test==2)
  {
    s->uinf = s->mach*cos(s->alpha*deg2rad)*cos(s->beta*deg2rad);
    s->winf = s->mach*cos(s->alpha*deg2rad)*sin(s->beta*deg2rad);
    s->vinf = s->mach*sin(s->alpha*deg2rad);
  }

  s->einf = pinf/(gamm-1)+0.5*rinf*(s->uinf*s->uinf+
                                    s->vinf*s->vinf+s->winf*s->winf);
  tracef(s->uinf);
  tracef(s->vinf);
  tracef(s->winf);
  s->gm1=gamm-1;
  s->c2b=0.3678;
  s->rgas=1./gamm;
  s->pr=0.72;
  s->prtr=0.3333;
  //
  // allocate all the solution arrays
  //
  s->q=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->qt=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->qtt=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->r=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->dq=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->ddq=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->ddqb=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->ddqf=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->dtac=(double *) malloc(sizeof(double)*g->ncells);
  s->r0=(double *) malloc(sizeof(double)*NVAR*g->ncells);
  s->D=(double ***) malloc(sizeof(double **)*g->ncells);
  s->itag=(int *) malloc(sizeof(int)*g->ncells);
  g->ff=(faceMat *) malloc(sizeof(faceMat) *g->nfaces);

  //
  for(i=0;i<g->ncells;i++)
    {
      s->D[i]=(double **) malloc(sizeof(double *)*NQ);
      for(j=0;j<NQ;j++)
	s->D[i][j]=(double *) malloc(sizeof(double)*NQ);
    }
  s->sigma=(double *) malloc(sizeof(double)*g->ncells);
  //

  m=0;

  // far-field bc
  if(g->test != 1) 
  {
    //read restart files
    if(irest==1) 
    {
      fp = fopen("QuadData/urest.tc1","r");
      for(i=0;i<g->ncells;i++)
      {
        fscanf(fp,"%lf %lf %lf %lf %lf\n",&(s->q[NVAR*i]),&(s->q[NVAR*i+1]),
        &(s->q[NVAR*i+2]),&(s->q[NVAR*i+3]),&(s->q[NVAR*i+4]));
      
        s->qt[NVAR*i]   = s->q[NVAR*i];
        s->qt[NVAR*i+1] = s->q[NVAR*i+1];
        s->qt[NVAR*i+2] = s->q[NVAR*i+2];
        s->qt[NVAR*i+3] = s->q[NVAR*i+3];
        s->qt[NVAR*i+4] = s->q[NVAR*i+4];
      }
      fscanf(fp,"%d\n",&(s->nt));
      fclose(fp);
      printf("#ham3d: Reading restart files\n");

    }
    else
    {
      for (i=0;i<g->ncells;i++)
      {
        //
        // initialize freestream conservative variables
        //
        s->q[m]=rinf;
        s->qt[m]=s->q[m];
        m++;
        s->q[m]=rinf*s->uinf;
        s->qt[m]=s->q[m];
        m++;
        s->q[m]=rinf*s->vinf;
        s->qt[m]=s->q[m];
        m++;
        s->q[m]=rinf*s->winf;
        s->qt[m]=s->q[m];
        m++;
        s->q[m] = s->einf;
        s->qt[m]=s->q[m];
        m++;
      }
    }
  }

    //periodic boundary case
    if(g->test == 1)     
    {
      if(irest==1)
      {
        fp = fopen("QuadData/urest.tc1","r");
        for(i=0;i<g->ncells;i++)
        {
          fscanf(fp,"%lf %lf %lf %lf %lf\n",&(s->q[NVAR*i]),&(s->q[NVAR*i+1]),
          &(s->q[NVAR*i+2]),&(s->q[NVAR*i+3]),&(s->q[NVAR*i+4]));
        
          s->qt[NVAR*i]   = s->q[NVAR*i];
          s->qt[NVAR*i+1] = s->q[NVAR*i+1];
          s->qt[NVAR*i+2] = s->q[NVAR*i+2];
          s->qt[NVAR*i+3] = s->q[NVAR*i+3];
          s->qt[NVAR*i+4] = s->q[NVAR*i+4];
        }
        fscanf(fp,"%d\n",&(s->nt));
        fclose(fp);
        printf("#ham3d: Reading restart files\n");
      }
      else
      {
        for(i=0;i<g->ncells;i++)
        {
          n1 = g->conn[8*i]-1;
          n2 = g->conn[8*i+1]-1;
          n3 = g->conn[8*i+2]-1;
          n4 = g->conn[8*i+3]-1;
          n5 = g->conn[8*i+4]-1;
          n6 = g->conn[8*i+5]-1;
          n7 = g->conn[8*i+6]-1;
          n8 = g->conn[8*i+7]-1;

          xc = (g->x[3*n1]+g->x[3*n2]+g->x[3*n3]+g->x[3*n4]
                +g->x[3*n5]+g->x[3*n6]+g->x[3*n7]+g->x[3*n8])/8.;
          yc = (g->x[3*n1+1]+g->x[3*n2+1]+g->x[3*n3+1]+g->x[3*n4+1]
                +g->x[3*n5+1]+g->x[3*n6+1]+g->x[3*n7+1]+g->x[3*n8+1])/8.;
          zc = (g->x[3*n1+2]+g->x[3*n2+2]+g->x[3*n3+2]+g->x[3*n4+2]
                +g->x[3*n5+2]+g->x[3*n6+2]+g->x[3*n7+2]+g->x[3*n8+2])/8.;

          pi  = deg2rad*180.;
          gm1 = gamm-1.;
          tp0 = 1.0/gamm;
    
          s->uinf = 0.5;
          s->vinf = 0.0;
          s->winf = 0.0;
          vstr    = 5.0;

          xv0     = 5.;
          yv0     = 5.0;
          zv0     = 0.;
    
          c1 = 0.5*vstr/pi;
          c2 = 0.125*gm1/(gamm*gamm)*pow(vstr,2)/pow(pi,2);//gamm^2?
    
          xv = 0;
          yv = 0;
          zv = 0;
    
          xx = xc - (xv0+xv);
          yy = yc - (yv0+yv);
          zz = zc - (zv0+zv);
    
          r2 = xx*xx + yy*yy;
          //r2 = xx*xx + zz*zz;

          rr = sqrt(r2);
    
          du = c1*exp(0.5*(1.-r2))*(-yy);
          dv = c1*exp(0.5*(1.-r2))*(xx);
          dw = 0;   
          //du = c1*exp(0.5*(1.-r2))*(-zz);
          //dv = 0.;
          //dw = c1*exp(0.5*(1.-r2))*(xx);
 

          dtp = -c2*exp(1.-r2);    
          tp  = tp0 + dtp;
    
          r   = pow(gamm*tp,1.0/gm1);
    
          s->uinf = s->uinf + du;
          s->vinf = s->vinf + dv;
          s->winf = s->winf + dw;

          p   = pow(r,gamm);

          s->q[m] = r;
          s->qt[m] = s->q[m];
          m++;
    
          s->q[m] = r*s->uinf;
          s->qt[m] = s->q[m];
          m++;
    
          s->q[m] = r*s->vinf;
          s->qt[m] = s->q[m];
          m++;

          s->q[m] = r*s->winf;
          s->qt[m] = s->q[m];
          m++;

          s->q[m] = p/(gamm-1)+0.5*r*(s->uinf*s->uinf+s->vinf*s->vinf
                                      +s->winf*s->winf);
          s->qt[m]=s->q[m]; //what is this?
          m++;
          }
       } //irest
    } //if(test==1)

}
