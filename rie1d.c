// ##################################################################
//
// rie1d.c
//
// Determine far-field boundary data using quasi 1-d 
// characteristic relations.
// Written by Yong Su Jung
// ##################################################################

#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void rie1d(GRID *g, SOLN *s)
{
  int i,j;
  double rin,uin,vin,win,pin;
  double rfr,ufr,vfr,wfr,pfr;
  double dx,dy,dz;
  double ain,afr,enin,enfr,ubin,ubfr;
  double rplus,rminus,ubar,abar;
  double ub,vb,wb,enb,uu;
  double u,v,w,rho,p,gm1,xgm1;
  
  gm1 = gamm-1;
  xgm1 = 1.0/gm1;

  // calculate R+ R-
  rin = s->qb[0];
  uin = s->qb[1]/rin;
  vin = s->qb[2]/rin;
  win = s->qb[3]/rin;
  pin = gm1*(s->qb[4]-0.5*(s->qb[1]*s->qb[1]+s->qb[2]*s->qb[2]+s->qb[3]*s->qb[3])/rin);
  rfr = s->qb[5];
  ufr = s->qb[6];
  vfr = s->qb[7];
  wfr = s->qb[8];
  pfr = s->qb[9];
  dx  = s->qb[10];
  dy  = s->qb[11];
  dz  = s->qb[12];
  ain = sqrt(gamm*pin/rin);
  afr = sqrt(gamm*pfr/rfr);
  enin = pin*pow(rin,-gamm);
  enfr = pfr*pow(rfr,-gamm);
  ubin = dx*uin + dy*vin + dz*win;
  ubfr = dx*ufr + dy*vfr + dz*wfr;

  rplus = ubin + 2.0*ain*xgm1;
  rminus = ubfr - 2.0*afr*xgm1;
  ubar = 0.5*(rplus+rminus);
  abar = 0.25*gm1*(rplus-rminus);


  // if unorm > 0 this is outflow: take variables from inside
  // if unorm > 0 this is outflow: take variables from inside
  
  if(ubar>0.0) 
  {
    ub = uin;
    vb = vin;
    wb = win;
    enb = enin;
    uu = ubin;
  }
  else
  {
    ub = ufr;
    vb = vfr;
    wb = wfr;
    enb = enfr;
    uu = ubfr;
  }


  rho = pow((abar*abar/(gamm*enb)),xgm1);
  p   = rho*abar*abar/gamm;
  u = ub + dx*(ubar - uu);
  v = vb + dy*(ubar - uu);
  w = wb + dz*(ubar - uu);

  s->qb[0] = rho;
  s->qb[1] = u*rho;
  s->qb[2] = v*rho;
  s->qb[3] = w*rho; 
  s->qb[4] = p/gm1 + 0.5*rho*(u*u+v*v+w*w);
}
