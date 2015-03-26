// #######################################################
//
// ham3d.c
//
// Written dy Yong Su Jung
// University of Marylnad
// ########################################################
// HAM3D V6
// 1. strand grid is set
// 2. Inviscid solid surface boundary condition is set
// 3. Only explicit time integration can be used
// 4. It does not need nearbody Radius
// 5. periodic boundary condition is set
// 6. Add output option
// 7. integrate 3D wing case 2014.12.17
// HAM3D V8
// 1. Add restart option
// 2. Laminar Viscous term add (explicit, implicit)
//======================================================
//bug fix report
//
//1. computeRHS.c : line 357 (2014.11.19)
// else if (rightCell==-2) 
// -> else if(rightCell==-2 && g->test==0)
//2. computeRHS.c: line339
// if(rightCell>-1||g->test==1)->if(rightCell>-1)
//
//3. same as computeRHSk.c : line 310, line:325, line: 396
//4. outputSolution.c fix
//4. apply_periodic_LHS.c bug fixed : 2014.11.21
//5. jac_roe.f90 : add the muliply at last line
//6. error fixed in apply_periodic_LHS.c : 2014.12.17
//======================================================
//############################################################

#include <stdio.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <time.h>

int main()
{
  GRID *g;
  SOLN *s;
  int ngrids;
  int nsteps, ntotal;
  int nwrite;
  int i,n,nn;
  int irest,nrest; //restart
  double dt;
  double CFL;
  double l2rho,linfrho;
  double sref;
  int msweep;
  char c;
  char fname[20];
  char scheme[20];
  char timeInteg[20];
  int order,timeacc,visc,test,outform;
  FILE *fp;
  clock_t start, end;
  double cpu_time_used;
  // 
  fp=fopen("input.ham3d","r");
  while((c=fgetc(fp))!='\n');
  while((c=fgetc(fp))!='\n');
  while((c=fgetc(fp))!='\n');
  while((c=fgetc(fp))!='\n');
  fscanf(fp,"scheme=%s\n",scheme);
  fscanf(fp,"time integration=%s\n",timeInteg);
  fscanf(fp,"order=%d\n",&order);
  fscanf(fp,"timeacc=%d\n",&timeacc);
  fscanf(fp,"nsteps=%d\n",&nsteps);
  fscanf(fp,"nwrite=%d\n",&nwrite);
  fscanf(fp,"dt=%lf\n",&dt);
  fscanf(fp,"CFL=%lf\n",&CFL);
  fscanf(fp,"msweep=%d\n",&msweep);
  fscanf(fp,"visc=%d\n",&visc);
  fscanf(fp,"testcase=%d\n",&test);
  fscanf(fp,"irest=%d\n",&irest);
  fscanf(fp,"nrest=%d\n",&nrest);
  fscanf(fp,"output=%d\n",&outform);
  fclose(fp);
  trace(nsteps);
  tracef(dt);

  //
  ngrids=1;
  g=(GRID *) malloc(sizeof(GRID)*ngrids);
  s=(SOLN *) malloc(sizeof(SOLN)*ngrids);
  //
  // preprocess grids
  // code is written with an overset
  // framework in mind (but not implemented yet)
  // (ngrid==1) for now
  //

  for(i=0;i<ngrids;i++) 
    {
      
      g[i].visc=visc;
      readGrid(&g[i]); //need to strand grid reading and make connectivity
      
      g[i].test=test;
      s[i].outform=outform;
      
      preprocess(&g[i]);
      initflow(&g[i],&s[i],irest);

      g[i].order=order;
      g[i].timeacc=timeacc;
      g[i].CFL=CFL;  
      s[i].cflnum=CFL;

      g[i].msweep=msweep;

      if (strcmp(timeInteg,"bdf1")==0) g[i].timeInteg=BDF1;
      if (strcmp(timeInteg,"bdf2")==0) g[i].timeInteg=BDF2;
    }
    tracef(CFL)
 
  // now run chosen number of time steps
  //
  if(irest==0) 
  {
    fp = fopen("./output/sol_his.dat","w");
    fclose(fp);
  }

  printf("#ham3d : using %s scheme for inversion\n",scheme);
  cpu_time_used=0;

  // main iteration
  ntotal = 0;
  if(irest==1) ntotal = s->nt;
  for (n=ntotal;n<ntotal+nsteps;n++) 
  { 
     fp = fopen("./output/sol_his.dat","a+");
     for(i=0;i<ngrids;i++) 
     { 
	    start=clock();
	    
	    stepSolution(scheme,&g[i],&s[i],dt,&l2rho,&linfrho);
	    
	    computeForce(&g[i],&s[i]);
	    
	    end = clock();
	    
	    cpu_time_used+=(((double) (end - start)) / CLOCKS_PER_SEC);
       
       printf("%d %e %e %2.4f %2.4f %2.4f\n",n,l2rho,linfrho,s[i].cl
               ,s[i].cd,cpu_time_used);
       fprintf(fp,"%d %e %e %2.4f %2.4f %2.4f\n",n,l2rho,linfrho
               ,s[i].cl,s[i].cd,cpu_time_used);
       fclose(fp);
	  }
     //write solution
     if((n+1)%nwrite==0||n==ntotal)
     { 
       nn = (n+1)/nwrite;
       outputSolution(&g[0],&s[0],nn); 
     }

     //write restart
     if((n+1)%nrest==0)
     {
       nn = (n+1)/nrest;
       wrest(&g[0],&s[0],n,nn); 
     }
  } 
}

//######################################################################
// END OF PROGRAM
//######################################################################

