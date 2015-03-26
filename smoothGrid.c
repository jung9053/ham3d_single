#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void smoothGrid(GRID *g, int msweep)
{
  int **indx;
  int *nodeCount;
  int *iflag;
  int i,m,j;
  int node1,node2;
  double *x1;
  double norm;
  //


  indx=(int **) malloc(sizeof(int *)*g->nnodes);
  iflag=(int *) malloc(sizeof(int)*g->nnodes);
  nodeCount=(int *) malloc(sizeof(int)*g->nnodes);
  x1=(double *) malloc(sizeof(double)*2*g->nnodes);
  //
  for(i=0;i<g->nnodes;i++) iflag[i]=nodeCount[i]=0;
  //
  for(i=0;i<g->nfaces;i++)
    {
      node1=g->faces[6*i];
      node2=g->faces[6*i+1];
      nodeCount[node1]++;
      nodeCount[node2]++;      
      if (g->faces[6*i+4] == -1) 
	{
	  iflag[node1]=1;
	  iflag[node2]=1;
	}
    }


  for(i=0;i<g->nnodes;i++)
    indx[i]=(int *)malloc(sizeof(int)*nodeCount[i]);

  // why nodecount set 0  again?
  //
  for(i=0;i<g->nnodes;i++) nodeCount[i]=0;

  //
  for(i=0;i<g->nfaces;i++)
    {
      node1=g->faces[6*i];
      node2=g->faces[6*i+1];
      indx[node1][nodeCount[node1]]=node2;
      nodeCount[node1]++;
      indx[node2][nodeCount[node2]]=node1;
      nodeCount[node2]++;
    }

  // except surface node, inner nodes location relocated 
  // by averaging neighbor nodes locations 
  //
  for(m=0;m<msweep;m++)
    {      
      for(i=0;i<g->nnodes;i++)
	{
	  if (iflag[i]==0)
	    {
	      x1[2*i]=x1[2*i+1]=0;
	      for(j=0;j<nodeCount[i];j++)                 
		{
		  x1[2*i]+=g->x[2*indx[i][j]];
		  x1[2*i+1]+=g->x[2*indx[i][j]+1];
		}
	      x1[2*i]/=nodeCount[i];
	      x1[2*i+1]/=nodeCount[i];
	    }
	  else
	    {
	      x1[2*i]=g->x[2*i];
	      x1[2*i+1]=g->x[2*i+1];
	    }
	}
      norm=0.0;

      for(i=0;i<2*g->nnodes;i++)
	{
	  norm+=(g->x[i]-x1[i])*(g->x[i]-x1[i]);
	 
          // array containing the original data points (g->x) are
          // re-written based on the "smoothed" grid positions x1
          // Therefore, the output.plt may vary as compared to the
          // input data points
          g->x[i]=x1[i]; 
	}
      norm=sqrt(norm)/2/g->nnodes;
      //tracef(norm);
    }

  free(nodeCount);
  free(x1);
  for(i=0;i<g->nnodes;i++)
    free(indx[i]);
  free(indx);
  free(iflag);
}
