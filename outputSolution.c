#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void outputSolution(GRID *g,SOLN *s,int nn)
{
  int i,j,n;
  int iface,node1,node2,node3,node4;
  int icell;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,rho,rhou,rhov,rhow,e,pp,cp;
  FILE *fp,*fp1;
  char fname[80];

  //
//  fp=fopen("output/output.plt","w");
  if        (nn < 10   ) {sprintf(fname,"./output/vol00%d.dat",nn);}
  else if   (nn < 100  ) {sprintf(fname,"./output/vol0%d.dat",nn);}
  else if   (nn < 1000 ) {sprintf(fname,"./output/vol%d.dat",nn);}
  fp = fopen(fname,"w");

  //default format
  if(s->outform==0 || g->test==2) //3D case only can this NOW
  {
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"RHOU\",\"RHOV\",\"RHOW\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEBRICK N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i+1]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i+2]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      {
       fprintf(fp,"%f\n",s->q[5*i+n]);
      }
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d %d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],g->conn[8*i+2],g->conn[8*i+3],g->conn[8*i+4],g->conn[8*i+5],g->conn[8*i+6],g->conn[8*i+7]);

  fclose(fp);
  }


  if(s->outform==1 && g->test!=2) //3D cannot this now!
  {
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"U\",\"V\",\"W\",\"CP\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEBRICK N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i+1]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i+2]);
  fprintf(fp,"\n");

  //rho
  for(i=0;i<g->ncells;i++)
     {
      rho = s->q[5*i];
      fprintf(fp,"%f\n",rho);
     }
  // u
  for(i=0;i<g->ncells;i++)
     {
      rho = s->q[5*i];
      rhou = s->q[5*i+1];
      fprintf(fp,"%f\n",rhou/rho);
     }
  // v
  for(i=0;i<g->ncells;i++)
     {
      rho = s->q[5*i];
      rhov = s->q[5*i+2];
      fprintf(fp,"%f\n",rhov/rho);
     }
  // w
  for(i=0;i<g->ncells;i++)
     {
      rho = s->q[5*i];
      rhow = s->q[5*i+3];
      fprintf(fp,"%f\n",rhow/rho);
     }
  // cp
  for(i=0;i<g->ncells;i++)
      {
       rho = s->q[5*i];
       rhou = s->q[5*i+1];
       rhov = s->q[5*i+2];
       rhow = s->q[5*i+3];
       e = s->q[5*i+4];
       pp=(gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
       cp=(pp-pinf)/(0.5*s->mach*s->mach);	
       fprintf(fp,"%f\n",cp);
      }

  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d %d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],g->conn[8*i+2],g->conn[8*i+3],g->conn[8*i+4],g->conn[8*i+5],g->conn[8*i+6],g->conn[8*i+7]);

  fclose(fp);
  }
 
//======================================================================
  // output for surface boundary
  if        (nn < 10   ) {sprintf(fname,"./output/surf00%d.dat",nn);}
  else if   (nn < 100  ) {sprintf(fname,"./output/surf0%d.dat",nn);}
  else if   (nn < 1000 ) {sprintf(fname,"./output/surf%d.dat",nn);}
  fp = fopen(fname,"w");
 
  //default
  if(s->outform==0 && g->test!=2)
  {
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"RHOU\",\"RHOV\",\"RHOW\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes/g->nstrand,g->nbfaces); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i+1]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i+2]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->nbfaces;i++)
      {
       fprintf(fp,"%f\n",s->q[5*i+n]);
      }
  fprintf(fp,"\n");
  for(i=0;i<g->nbfaces;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],g->conn[8*i+2],g->conn[8*i+3]);

  fclose(fp);
  }
  
  //new
  if(s->outform==1 && g->test!=2)
  {
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"U\",\"V\",\"W\",\"CP\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes/g->nstrand,g->nbfaces); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i+1]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i+2]);
  fprintf(fp,"\n");
  // rho
  for(i=0;i<g->nbfaces;i++)
    {
     rho = s->q[5*i];
     fprintf(fp,"%f\n",rho);
    }
  // u
  for(i=0;i<g->nbfaces;i++)
    {
     rho = s->q[5*i];
     rhou = s->q[5*i+1];
     fprintf(fp,"%f\n",rhou/rho);
    }
  // v
  for(i=0;i<g->nbfaces;i++)
    {
     rho = s->q[5*i];
     rhov = s->q[5*i+2];
     fprintf(fp,"%f\n",rhov/rho);
    }
  // w
  for(i=0;i<g->nbfaces;i++)
    {
     rho = s->q[5*i];
     rhow = s->q[5*i+3];
     fprintf(fp,"%f\n",rhow/rho);
    }
  // cp
  for(i=0;i<g->nbfaces;i++)
      {
       rho = s->q[5*i];
       rhou = s->q[5*i+1];
       rhov = s->q[5*i+2];
       rhow = s->q[5*i+3];
       e = s->q[5*i+4];
       pp=(gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
       cp=(pp-pinf)/(0.5*s->mach*s->mach);	
       fprintf(fp,"%f\n",cp);
      }
  fprintf(fp,"\n");
  for(i=0;i<g->nbfaces;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],g->conn[8*i+2],g->conn[8*i+3]);

  fclose(fp);
  }


  //3D wing
  if(g->test == 2)
  {
fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"RHOU\",\"RHOV\",\"RHOW\",\"E\"\n");
fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nbfaces*4,g->nbfaces); 
fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
  
  for(i=0;i<g->nbfaces;i++) // x 
  {
     iface = g->bfaces[i];
     for(j=0;j<4;j++)
     {
       n = g->faces[8*iface+j];
       fprintf(fp,"%f\n",g->x[3*n]);
     }
  }
       fprintf(fp,"\n");
  for(i=0;i<g->nbfaces;i++)  // y
  {
     iface = g->bfaces[i];
     for(j=0;j<4;j++)
     {
       n = g->faces[8*iface+j];
       fprintf(fp,"%f\n",g->x[3*n+1]);
     }
  }
       fprintf(fp,"\n");
  for(i=0;i<g->nbfaces;i++) // z
  {
     iface = g->bfaces[i];
     for(j=0;j<4;j++)
     {
       n = g->faces[8*iface+j];
       fprintf(fp,"%f\n",g->x[3*n+2]);
     }
  }
       fprintf(fp,"\n");
  //flow variables
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->nbfaces;i++)
      {
       iface = g->bfaces[i];
       icell = g->faces[8*iface+4];
       fprintf(fp,"%f\n",s->q[5*icell+n]);
      }
  fprintf(fp,"\n");

 for(i=0;i<g->nbfaces;i++) //connectivity
 {   
    fprintf(fp,"%d %d %d %d\n",4*i+1,4*i+2,4*i+3,4*i+4);
 }
  fclose(fp);
 }








//======================================================================
//  fp=fopen("./output/cp.dat","w");
//  for(i=0;i<g->nbfaces;i++)
//   {
//    iface=g->bfaces[i];
//    node1=g->faces[8*iface];
//    node2=g->faces[8*iface+1];
//    node3=g->faces[8*iface+2];
//    node4=g->faces[8*iface+3];
//
//    x1=g->x[3*node1];
//    y1=g->x[3*node1+1];
//    z1=g->x[3*node1+2];
//    x2=g->x[3*node2];
//    y2=g->x[3*node2+1];
//    z2=g->x[3*node2+2];
//    x3=g->x[3*node3];
//    y3=g->x[3*node3+1];
//    z3=g->x[3*node3+2];
//    x4=g->x[3*node4];
//    y4=g->x[3*node4+1];
//    z4=g->x[3*node4+2];
//
//    icell = g->faces[8*iface+4];
//    rho   = s->q[NVAR*icell];
//    rhou  = s->q[NVAR*icell+1];
//    rhov  = s->q[NVAR*icell+2];
//    rhow  = s->q[NVAR*icell+3];
//    e     = s->q[NVAR*icell+4];
//    pp=(gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
//    cp=(pp-pinf)/(0.5*s->mach*s->mach);	
//    fprintf(fp,"%.16e %.16e %.16e\n",(x1+x2+x3+x4)*0.25,(y1+y2+y3+y4)*0.25,cp);
//  }
//  fclose(fp);
}





void outputdq(GRID *g,SOLN *s)
{
  int i,n;
  int iface,node1,node2;
  int icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE *fp;
  char fname[80];
  static int istep0=0;
  //
  sprintf(fname,"dq%d.plt",istep0);
  fp=fopen(fname,"w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%f\n",s->ddq[4*i+n]);
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);
  istep0++;
}

void outputr(GRID *g,SOLN *s)
{
  int i,n;
  int iface,node1,node2;
  int icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE *fp;
  char fname[80];
  static int istep0=0;
  //
  sprintf(fname,"r%d.plt",istep0);
  fp=fopen(fname,"w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%f\n",s->r[4*i+n]);
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);
  istep0++;
}


