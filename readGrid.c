// ##################################################################
//
// readGrid.c
//
// Function that reads the grid based on the outputs of the Matlab
// mesh generation code
//
// Data files read from:
//  - coord.dat   (data coordinate points)
//  - conn.dat    (node connectivity information for the triangles)
//  - qedges.dat  (writes the QEdge matrix, 6 cols, refer documentation)
//  - ncolors.dat (number of loops of each colour)
//  - iqloops.dat (index of inner and outer loops of each node)
//  - qloops.dat  (Inner and outer loops  for each node)
//
// Written by Dr. Jayanarayanan Sitaraman
// Modified by Yong Su Jung
// ##################################################################

#include <stdio.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <stdlib.h>
#define NQ 5

void readGrid(GRID *g)
{
  FILE *fp;
  char c;
  int i,l1,l2,l3,j,ii;
  int workSize;
  double fdummy;

  //
  // read coordinates
  //
  fp=fopen("QuadData/coord.dat","r");
  fscanf(fp,"%d",&(g->nnodes));
  g->x=(double *) malloc(sizeof(double)*3*g->nnodes);

  for(i=0;i<g->nnodes;i++)
    fscanf(fp,"%lf %lf %lf\n",&(g->x[3*i]),&(g->x[3*i+1]),&(g->x[3*i+2]));
  fclose(fp);

 
  //
  // read connectivity
  //
  fp=fopen("QuadData/conn.dat","r");
  fscanf(fp,"%d",&(g->ncells));

  g->conn=(int *) malloc(sizeof(int)*8*g->ncells);

  for(i=0;i<g->ncells;i++)    
    fscanf(fp,"%d %d %d %d %d %d %d %d\n",
    &(g->conn[8*i]),&(g->conn[8*i+1]),
    &(g->conn[8*i+2]),&(g->conn[8*i+3]),
    &(g->conn[8*i+4]),&(g->conn[8*i+5]),
    &(g->conn[8*i+6]),&(g->conn[8*i+7]));
  fclose(fp);
 
  //
  // read faces
  //
  fp=fopen("QuadData/ofaces.dat","r");
  fscanf(fp,"%d",&(g->nfaces));
  g->faces=(int *) malloc(sizeof(int)*8*g->nfaces);
  for(i=0;i<g->nfaces;i++)
   {
    fscanf(fp,"%d %d %d %d %d %d %d %d\n",&(g->faces[8*i]),//n1
	   &(g->faces[8*i+1]),//n2
	   &(g->faces[8*i+2]),//n3
	   &(g->faces[8*i+3]),//n4
	   &(g->faces[8*i+4]),//c1
	   &(g->faces[8*i+6]),//c2
	   &(g->faces[8*i+5]),//e1
	   &(g->faces[8*i+7]));//e2

      // boundary cell should be -1 at upperlayers
      // this is due to just add the total number of faces at each layer
      // to the previous layer face index
      if(g->faces[8*i+7]==-1){
         g->faces[8*i+6]=-1;
      }

     //
     // swap for boundary cells(-1: far boundary, -2: solid wall)  
     //
     if (g->faces[8*i+4]==-1 || g->faces[8*i+4]==-2) 
	  {
         swap(g->faces[8*i]   , g->faces[8*i+1]);//n2,n4 
         swap(g->faces[8*i+2] , g->faces[8*i+3]);//n2,n4 

	      swap(g->faces[8*i+4] , g->faces[8*i+6]);//cell index
         swap(g->faces[8*i+5] , g->faces[8*i+7]);//element index  
	  }
    }
  fclose(fp);

  //
  // read colors
  //
  fp=fopen("QuadData/ncolors.dat","r");
  fscanf(fp,"%d",&(g->ncolors));
  g->chainsPerColor=(int *) malloc(sizeof(int)*g->ncolors);
  for(i=0;i<g->ncolors;i++)
    fscanf(fp,"%d\n",&(g->chainsPerColor[i]));
  fclose(fp);


  //
  // read chain information (chain size)
  //
  fp=fopen("QuadData/iqloops.dat","r");
  fscanf(fp,"%d\n",&(g->nchains));
  g->nchains--;
  g->faceStartPerChain=(int *) malloc(sizeof(int)*((g->nchains+1)));
  g->nmaxchain=0;

  ii=0;
  for(i=0;i<g->nchains+1;i++)
    {      
      fscanf(fp,"%d\n",&(g->faceStartPerChain[ii]));
        
      // this routine is for delete the duplicated one
      // at every layers
      if (i > 0) 
	   {
        if(g->faceStartPerChain[ii]==g->faceStartPerChain[ii-1])
        {    
          ii--;
        }
	  g->nmaxchain=max(g->nmaxchain,g->faceStartPerChain[ii]-g->faceStartPerChain[ii-1]);
	  }
    ii++;
   }
  fclose(fp);
  g->nchains =ii-1;

  //
  trace(g->nmaxchain);
  //
  workSize=g->nmaxchain+5;
  g->ql=(double **)malloc(sizeof(double *)*(workSize));
  g->qr=(double **)malloc(sizeof(double *)*(workSize));
  g->dql=(double **)malloc(sizeof(double *)*(workSize));
  g->dqr=(double **)malloc(sizeof(double *)*(workSize));
  g->f=(double **) malloc(sizeof(double *)*(workSize));
  g->fv=(double **) malloc(sizeof(double *)*(workSize));
  g->df=(double **) malloc(sizeof(double *)*(workSize));
  g->f2=(double **) malloc(sizeof(double *)*(workSize));
  g->cindx=(int *)malloc(sizeof(int)*(workSize));
  g->ctype=(int *)malloc(sizeof(int)*(workSize));
  g->A=(double ***) malloc(sizeof(double **)*(workSize));
  g->B=(double ***) malloc(sizeof(double **)*(workSize));
  g->C=(double ***) malloc(sizeof(double **)*(workSize));  
  g->F=(double **) malloc(sizeof(double *)*(workSize));
  g->Q=(double **) malloc(sizeof(double *)*(workSize));
  //
  for(i=0;i<workSize;i++)
  {
    g->ql[i]=(double *)malloc(sizeof(double)*NVAR);
    g->qr[i]=(double *)malloc(sizeof(double)*NVAR);
    g->dql[i]=(double *)malloc(sizeof(double)*NVAR);
    g->dqr[i]=(double *)malloc(sizeof(double)*NVAR);
    g->f[i]=(double *)malloc(sizeof(double)*NVAR);
    g->fv[i]=(double *)malloc(sizeof(double)*NVAR);
    g->df[i]=(double *)malloc(sizeof(double)*NVAR);
    g->f2[i]=(double *)malloc(sizeof(double)*NVAR);
    g->A[i]=(double **) malloc(sizeof(double *)*NQ);
    g->B[i]=(double **) malloc(sizeof(double *)*NQ);
    g->C[i]=(double **) malloc(sizeof(double *)*NQ);
    g->F[i]=(double *) malloc(sizeof(double)*NQ);
    g->Q[i]=(double *) malloc(sizeof(double)*NQ);
    for(j=0;j<NQ;j++)
    {
      g->A[i][j]=(double *)malloc(sizeof(double)*NQ);
      g->B[i][j]=(double *)malloc(sizeof(double)*NQ);
      g->C[i][j]=(double *)malloc(sizeof(double)*NQ);
    }
  }
  //
  // read chain information (chain connecitivity)
  //
  fp=fopen("QuadData/qloops.dat","r");
  fscanf(fp,"%d",&(g->nchainFaces));
  g->chainConn=(int *)malloc(sizeof(int) * g->nchainFaces);
  for(i=0;i<g->nchainFaces;i++)
    fscanf(fp,"%d",&(g->chainConn[i]));
  fclose(fp);

  fp=fopen("QuadData/nstrands.dat","r");
  fscanf(fp,"%d",&(g->nstrand));
  fclose(fp);

  printf("#ham3d: Finished reading files\n");
  trace(g->nnodes);
  trace(g->ncells);
  trace(g->nfaces);
  trace(g->ncolors);
  trace(g->nchains);
  trace(g->nchainFaces);
  trace(g->nstrand);
}

