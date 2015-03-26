#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define nearBodyRadius 3.0

void periodic_bc(GRID *g)
{
  int i,m,i1,i2,i3,i4,kk;
  int leftCell,rightCell;
  int leftFaceIndx,rightFaceIndx;
  int f,f1,f2;
  int icell,iface;
  double r1,r2,r3,r4;
  double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,dvedge;
  FILE *fp;
  double volmin,volmax;
  double r_mid[2],s_vec[2];
  double xa,xb,ya,yb,za,zb;
  
  int j,k,ff1,ff2,nbase,nbloop,tmp;
  double bf1n[10000],bf2n[10000],bf3n[10000],bf4n[10000];

//////////////////////////////////////////////////////////////////
  // Find far-boundary face (just for 1st layers) 
  //
  //       3
  //   |-------|
  //   |       |
  // 4 |       |2
  //   ---------
  //       1
  //
  //   top surface is 5,
  //   bottom surface is 6. 
  // strand grid is only alined with z-axis 
  g->m1=0;
  g->m2=0;
  g->m3=0;
  g->m4=0;
  g->m5=0;
  g->m6=0;
  
  // number of cells in each layer
  // nbase = number of side face for each layer
  // nbloop = nloops - strand
  nbase = (g->nfaces-(g->ncells)/(g->nstrand-1)*g->nstrand)/(g->nstrand-1);
  nbloop = g->nchains-(g->ncells)/(g->nstrand-1);

  for(i=0;i<nbase;i++)
  {
    if(g->faces[8*i+6] == -1)
    {
     i1 = g->faces[8*i];
     i2 = g->faces[8*i+2];

     if(g->x[3*i1+1]==g->x[3*i2+1]) // x parallel
     {
        if(g->x[3*i1+1] > 5.0) //3
        {
           g->bf3[g->m3] = i;
           bf3n[g->m3] = 0.5*(g->x[3*i1]+g->x[3*i2]);
            
           g->m3 = g->m3+1;
        }
        else //1
        {
           g->bf1[g->m1] = i;
           bf1n[g->m1] = 0.5*(g->x[3*i1]+g->x[3*i2]);
           g->m1 = g->m1+1;
        } 
     }
     if(g->x[3*i1] == g->x[3*i2]) // y parallel
     {
        if(g->x[3*i1]>5.0) //2
        {
           g->bf2[g->m2] = i;
           bf2n[g->m2] = 0.5*(g->x[3*i1+1]+g->x[3*i2+1]);
           g->m2 = g->m2+1;
        }
        else //4
        {
           g->bf4[g->m4] = i;
           bf4n[g->m4] = 0.5*(g->x[3*i1+1]+g->x[3*i2+1]);
           g->m4 = g->m4+1;
        }
     } //y parallel
    } //face = -1

  } // nface

    // m=1
    for(i=0;i<g->m1;i++)
    {
     for(j=i+1;j<g->m1;j++)
     {
      if(bf1n[i]>bf1n[j])
      {
       tmp = g->bf1[i];
       g->bf1[i] = g->bf1[j];
       g->bf1[j] = tmp;
      }
     }
    }
    
    // m=2
    for(i=0;i<g->m2;i++)
    {
     for(j=i+1;j<g->m2;j++)
     {
      if(bf2n[i]>bf2n[j])
      {
       tmp = g->bf2[i];
       g->bf2[i] = g->bf2[j];
       g->bf2[j] = tmp;
      }
     }
    }

    // m=3
    for(i=0;i<g->m3;i++)
    {
     for(j=i+1;j<g->m3;j++)
     {
      if(bf3n[i]>bf3n[j])
      {
       tmp = g->bf3[i];
       g->bf3[i] = g->bf3[j];
       g->bf3[j] = tmp;
      }
     }
    }

    // m=4
    for(i=0;i<g->m4;i++)
    {
     for(j=i+1;j<g->m4;j++)
     {
      if(bf4n[i]>bf4n[j])
      {
       tmp = g->bf4[i];
       g->bf4[i] = g->bf4[j];
       g->bf4[j] = tmp;
      }
     }
    }

    // m1
    for (j=0;j<g->m1;j++)
    {
      for(k=0;k<nbloop/(g->nstrand-1);k++)
      {
         f1 = g->faceStartPerChain[k];
         f2 = g->faceStartPerChain[k+1];
         ff1=g->chainConn[f1];
         ff2=g->chainConn[f2-1];

         if(g->bf1[j] == ff1) 
         {
          g->bf1_neigf[j][0] = g->chainConn[f1+1]; 
          g->bf1_neigf[j][1] = g->chainConn[f1+2]; 
          g->bf1_neigf[j][2] = g->chainConn[f1+3]; 
          //start
          for (kk=1;kk<g->nstrand-1;kk++)
          {
            f1  = g->faceStartPerChain[k+kk*nbloop/(g->nstrand-1)];
            g->bf1[j+(g->m1)*kk]          = g->chainConn[f1];
            g->bf1_neigf[j+(g->m1)*kk][0] = g->chainConn[f1+1];
            g->bf1_neigf[j+(g->m1)*kk][1] = g->chainConn[f1+2];
            g->bf1_neigf[j+(g->m1)*kk][2] = g->chainConn[f1+3];  

//    printf("face index = %d, n_face index = %d %d %d\n",g->bf1[j+(g->m1)*kk],g->bf1_neigf[j+g->m1*kk][0],g->bf1_neigf[j+g->m1*kk][1],g->bf1_neigf[j+g->m1*kk][2]); 
          }
          //end
          break;
         }
         if(g->bf1[j] ==ff2) 
         {
          g->bf1_neigf[j][0] = g->chainConn[f2-2]; 
          g->bf1_neigf[j][1] = g->chainConn[f2-3]; 
          g->bf1_neigf[j][2] = g->chainConn[f2-4]; 
          //start
          for (kk=1;kk<g->nstrand-1;kk++)
          {            
            f1  = g->faceStartPerChain[(k+1)+kk*nbloop/(g->nstrand-1)];
            g->bf1[j+(g->m1)*kk]          = g->chainConn[f1-1];
            g->bf1_neigf[j+(g->m1)*kk][0] = g->chainConn[f1-2];
            g->bf1_neigf[j+(g->m1)*kk][1] = g->chainConn[f1-3];
            g->bf1_neigf[j+(g->m1)*kk][2] = g->chainConn[f1-4];
//    printf("face index = %d, n_face index = %d %d %d\n",g->bf1[j+(g->m1)*kk],g->bf1_neigf[j+g->m1*kk][0],g->bf1_neigf[j+g->m1*kk][1],g->bf1_neigf[j+g->m1*kk][2]);
 
          }
          //end
          break;
         }   
      }

//    printf("face index = %d, n_face index = %d %d %d\n",g->bf1[j+g->m1],g->bf1_neigf[j+g->m1][0],g->bf1_neigf[j+g->m1][1],g->bf1_neigf[j+g->m1][2]);
      
      }
 
//exit(1);

    // m2
    for (j=0;j<g->m2;j++)
    {
    for(k=0;k<nbloop/(g->nstrand-1);k++)
    {
       f1 = g->faceStartPerChain[k];
       f2 = g->faceStartPerChain[k+1];
       ff1=g->chainConn[f1];
       ff2=g->chainConn[f2-1];
       if(g->bf2[j] == ff1) 
       {
        g->bf2_neigf[j][0] = g->chainConn[f1+1]; 
        g->bf2_neigf[j][1] = g->chainConn[f1+2]; 
        g->bf2_neigf[j][2] = g->chainConn[f1+3]; 
          //start
          for (kk=1;kk<g->nstrand-1;kk++)
          {
            f1  = g->faceStartPerChain[k+kk*nbloop/(g->nstrand-1)];
            g->bf2[j+(g->m2)*kk]          = g->chainConn[f1];
            g->bf2_neigf[j+(g->m2)*kk][0] = g->chainConn[f1+1];
            g->bf2_neigf[j+(g->m2)*kk][1] = g->chainConn[f1+2];
            g->bf2_neigf[j+(g->m2)*kk][2] = g->chainConn[f1+3];  
          }
        break;
       } 
       if(g->bf2[j] ==ff2) 
       {
        g->bf2_neigf[j][0] = g->chainConn[f2-2]; 
        g->bf2_neigf[j][1] = g->chainConn[f2-3]; 
        g->bf2_neigf[j][2] = g->chainConn[f2-4];
          //start
          for (kk=1;kk<g->nstrand-1;kk++)
          {
            f1  = g->faceStartPerChain[(k+1)+kk*nbloop/(g->nstrand-1)];
            g->bf2[j+(g->m2)*kk]          = g->chainConn[f1-1];
            g->bf2_neigf[j+(g->m2)*kk][0] = g->chainConn[f1-2];
            g->bf2_neigf[j+(g->m2)*kk][1] = g->chainConn[f1-3];
            g->bf2_neigf[j+(g->m2)*kk][2] = g->chainConn[f1-4];
          }
          //end 
        break;
       } 
    }
//    printf("face index = %d, n_face index = %d %d %d\n",g->bf2[j],g->bf2_neigf[j][0],g->bf2_neigf[j][1],g->bf2_neigf[j][2]);
    }
//exit(1);


    // m3
    for (j=0;j<g->m3;j++)
    {
    for(k=0;k<nbloop/(g->nstrand-1);k++)
    {
       f1 = g->faceStartPerChain[k];
       f2 = g->faceStartPerChain[k+1];
       ff1=g->chainConn[f1];
       ff2=g->chainConn[f2-1];
       if(g->bf3[j] == ff1) 
       {
        g->bf3_neigf[j][0] = g->chainConn[f1+1]; 
        g->bf3_neigf[j][1] = g->chainConn[f1+2]; 
        g->bf3_neigf[j][2] = g->chainConn[f1+3]; 
          for (kk=1;kk<g->nstrand-1;kk++)
          {
            f1  = g->faceStartPerChain[k+kk*nbloop/(g->nstrand-1)];
            g->bf3[j+(g->m3)*kk]          = g->chainConn[f1];
            g->bf3_neigf[j+(g->m3)*kk][0] = g->chainConn[f1+1];
            g->bf3_neigf[j+(g->m3)*kk][1] = g->chainConn[f1+2];
            g->bf3_neigf[j+(g->m3)*kk][2] = g->chainConn[f1+3];  
          }
        break;
       } 
       if(g->bf3[j] ==ff2) 
       {
        g->bf3_neigf[j][0] = g->chainConn[f2-2]; 
        g->bf3_neigf[j][1] = g->chainConn[f2-3]; 
        g->bf3_neigf[j][2] = g->chainConn[f2-4]; 
          for (kk=1;kk<g->nstrand-1;kk++)
          {
            f1  = g->faceStartPerChain[(k+1)+kk*nbloop/(g->nstrand-1)];
            g->bf3[j+(g->m3)*kk]          = g->chainConn[f1-1];
            g->bf3_neigf[j+(g->m3)*kk][0] = g->chainConn[f1-2];
            g->bf3_neigf[j+(g->m3)*kk][1] = g->chainConn[f1-3];
            g->bf3_neigf[j+(g->m3)*kk][2] = g->chainConn[f1-4];
          }
          //end 

        break;
       } 
    }
//    printf("face index = %d, n_face index = %d %d %d\n",g->bf3[j],g->bf3_neigf[j][0],g->bf3_neigf[j][1],g->bf3_neigf[j][2]);
    }
//exit(1);


    // m4
    for (j=0;j<g->m4;j++)
    {
    for(k=0;k<nbloop/(g->nstrand-1);k++)
    {
       f1 = g->faceStartPerChain[k];
       f2 = g->faceStartPerChain[k+1];
       ff1=g->chainConn[f1];
       ff2=g->chainConn[f2-1];
       if(g->bf4[j] == ff1) 
       {
        g->bf4_neigf[j][0] = g->chainConn[f1+1]; 
        g->bf4_neigf[j][1] = g->chainConn[f1+2]; 
        g->bf4_neigf[j][2] = g->chainConn[f1+3]; 
          for (kk=1;kk<g->nstrand-1;kk++)
          {
            f1  = g->faceStartPerChain[k+kk*nbloop/(g->nstrand-1)];
            g->bf4[j+(g->m4)*kk]          = g->chainConn[f1];
            g->bf4_neigf[j+(g->m4)*kk][0] = g->chainConn[f1+1];
            g->bf4_neigf[j+(g->m4)*kk][1] = g->chainConn[f1+2];
            g->bf4_neigf[j+(g->m4)*kk][2] = g->chainConn[f1+3];  
          }

        break;
       } 
       if(g->bf4[j] ==ff2) 
       {
        g->bf4_neigf[j][0] = g->chainConn[f2-2]; 
        g->bf4_neigf[j][1] = g->chainConn[f2-3]; 
        g->bf4_neigf[j][2] = g->chainConn[f2-4]; 
          for (kk=1;kk<g->nstrand-1;kk++)
          {
            f1  = g->faceStartPerChain[(k+1)+kk*nbloop/(g->nstrand-1)];
            g->bf4[j+(g->m4)*kk]          = g->chainConn[f1-1];
            g->bf4_neigf[j+(g->m4)*kk][0] = g->chainConn[f1-2];
            g->bf4_neigf[j+(g->m4)*kk][1] = g->chainConn[f1-3];
            g->bf4_neigf[j+(g->m4)*kk][2] = g->chainConn[f1-4];
          }

        break;
       } 
    }
//    printf("face index = %d, n_face index = %d %d %d\n",g->bf4[j],g->bf4_neigf[j][0],g->bf4_neigf[j][1],g->bf4_neigf[j][2]);
    }


//    for(i=0;i<g->m1*(g->nstrand-1);i++)
//    {
//      printf("n:%d face index:%d, %d %d %d\n",i+1,g->bf3[i],g->bf3_neigf[i][0],g->bf3_neigf[i][1],g->bf3_neigf[i][2]);
//    }
//    exit(1);


    // m5
    g->m5=g->ncells/(g->nstrand-1); 
    g->m6=g->ncells/(g->nstrand-1); 
    for(j=0;j<g->ncells/(g->nstrand-1);j++)
    {
      k =  nbloop + j;
      f1 = g->faceStartPerChain[k];
      f2 = g->faceStartPerChain[k+1];
     
      g->bf5_neigf[j][0] = g->chainConn[f1+1];
      g->bf5_neigf[j][1] = g->chainConn[f1+2];
      g->bf5_neigf[j][2] = g->chainConn[f1+3];

      g->bf6_neigf[j][0] = g->chainConn[f2-2];
      g->bf6_neigf[j][1] = g->chainConn[f2-3];
      g->bf6_neigf[j][2] = g->chainConn[f2-4];
//    printf("strand: %d %d %d\n",g->bf6_neigf[j][0],g->bf6_neigf[j][1],g->bf6_neigf[j][2]);
    }
//subroutine end
}

