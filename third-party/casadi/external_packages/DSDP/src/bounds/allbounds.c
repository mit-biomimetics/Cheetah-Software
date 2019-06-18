#include "dsdpcone_impl.h"
#include "dsdpsys.h"
#include "dsdp5.h"

/*!
  \file allbounds.c
  \brief Bound all the variables y in (D) and implement DSDPCone operations.
 */

#define COMPUTELBS(a,b,c)  ( (a)+(b)-(c))
#define COMPUTEUBS(a,b,c)  (-(a)-(b)-(c))

struct LUBounds_C{
  double r,muscale,minx;
  int invisible;
  int keyid;
  int setup,iter;
  double lbound,ubound;
  double maxratio;
  DSDPVec YD,YP,DYD;
  DSDP dsdp;
  DSDPTruth skipit;
  double pobj,pax,sumx,xdots;
};
typedef struct LUBounds_C* LUBounds;

#define LUKEY  5432

#define LUConeValid(a) {if (!(a)||((a)->keyid!=LUKEY)){ DSDPSETERR(101,"DSDPERROR: Invalid LUCone object\n");}}

static int LUBoundsS(void* dcone,DSDPVec Y, DSDPDualFactorMatrix flag, 
		     DSDPTruth *psdefinite);

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsSetUp"
static int LUBoundsSetup(void *dcone,DSDPVec y){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsSetUp2"
static int LUBoundsSetup2(void *dcone, DSDPVec Y, DSDPSchurMat M){
  int info;
  LUBounds lucone=(LUBounds)dcone;
  DSDPFunctionBegin;
  LUConeValid(lucone);
  if (lucone->setup==0){
    info=DSDPVecDuplicate(Y,&lucone->DYD);DSDPCHKERR(info);
    info=DSDPVecDuplicate(Y,&lucone->YD);DSDPCHKERR(info);
    info=DSDPVecDuplicate(Y,&lucone->YP);DSDPCHKERR(info);
    info=DSDPVecSet(lucone->lbound,lucone->YD);DSDPCHKERR(info);
    info=DSDPVecSetR(lucone->YD,-1e30);DSDPCHKERR(info);
    info=DSDPVecSetC(lucone->YD,-1e30);DSDPCHKERR(info);
    info=DSDPVecPointwiseMax(lucone->YD,Y,Y);DSDPCHKERR(info);
    info=DSDPVecSet(lucone->ubound,lucone->YD);DSDPCHKERR(info);
    info=DSDPVecSetR(lucone->YD,1e30);DSDPCHKERR(info);
    info=DSDPVecSetC(lucone->YD,1e30);DSDPCHKERR(info);
    info=DSDPVecPointwiseMin(lucone->YD,Y,Y);DSDPCHKERR(info);
    lucone->setup=1;
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsDestroy"
static int LUBoundsDestroy(void *dcone){
  int info;
  LUBounds lucone=(LUBounds)dcone;
  DSDPFunctionBegin;
  LUConeValid(lucone);
  info=DSDPVecDestroy(&lucone->DYD);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lucone->YD);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lucone->YP);DSDPCHKERR(info);
  DSDPFREE(&dcone,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsSize"
static int LUBoundsSize(void *dcone, double *n){
  int nn,info;
  LUBounds lucone=(LUBounds)dcone;
  DSDPFunctionBegin;
  LUConeValid(lucone);
  *n=0;
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  info=DSDPVecGetSize(lucone->YD,&nn);DSDPCHKERR(info);
  *n=(2*(nn-2)*lucone->muscale);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LUBoundsHessian"
static int LUBoundsHessian(void* dcone, double mu, DSDPSchurMat M,
			   DSDPVec vrhs1, DSDPVec vrhs2){
  int info,i,m;
  LUBounds lucone=(LUBounds)dcone;
  double assa,as,asrs=0,sl,su;
  double dd,yy,rs=0,srrs=0;
  double cc,rr,r,r0=lucone->r;
  double lbound, ubound;
  DSDPVec DScale=lucone->YP,Y=lucone->YD;

  DSDPFunctionBegin;
  LUConeValid(lucone);

  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  mu*=lucone->muscale;
  info=DSDPSchurMatDiagonalScaling(M,DScale);DSDPCHKERR(info);
  info=DSDPVecGetSize(DScale,&m);DSDPCHKERR(info);
  info=DSDPVecGetC(Y,&cc);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&rr);DSDPCHKERR(info);
  lbound=cc*lucone->lbound;
  ubound=cc*lucone->ubound;
  r=rr*lucone->r;
  info=DSDPVecSetC(DScale,0);DSDPCHKERR(info);
  info=DSDPVecSetR(DScale,0);DSDPCHKERR(info);
  for (i=1;i<m-1;i++){
    info=DSDPVecGetElement(DScale,i,&dd);DSDPCHKERR(info);
    info=DSDPVecSetElement(DScale,i,0);DSDPCHKERR(info);
    info=DSDPVecGetElement(Y,i,&yy);DSDPCHKERR(info);
    sl=1.0/COMPUTELBS(lbound,yy,r);
    su=1.0/COMPUTEUBS(ubound,yy,r);
    assa=mu*(su*su + sl*sl);
    as=mu*(su-sl);
    if (r){
      asrs=mu*r0*(su*su - sl*sl);
      rs+=(su+sl);
      srrs+=(su*su + sl*sl);
    }
    if (dd==0) continue;
    info=DSDPVecAddElement(vrhs2,i,dd*as);DSDPCHKERR(info);
    /*   info=DSDPVecAddElement(vrhs3,i,dd*asrs);DSDPCHKERR(info); */
    info=DSDPVecSetElement(DScale,i,dd*assa);DSDPCHKERR(info);
  }
  info=DSDPSchurMatAddDiagonal(M,DScale);DSDPCHKERR(info);
  info=DSDPVecAddR(vrhs2,r0*mu*rs);DSDPCHKERR(info);
  /*
  info=DSDPVecAddR(vrhs3,r0*mu*srrs);DSDPCHKERR(info);
  */
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsMultiply"
static int LUBoundsMultiply(void* dcone,  double mu, DSDPVec vrow, DSDPVec vin, DSDPVec vout){
  int info,i,m;
  double vv,ww,yy,sl,su,assa,rr;
  LUBounds lucone=(LUBounds)dcone;
  double cc, r=lucone->r;
  double lbound=lucone->lbound, ubound=lucone->ubound;
  DSDPVec Y=lucone->YD;

  DSDPFunctionBegin;
  LUConeValid(lucone);
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  mu*=lucone->muscale;
  info=DSDPVecGetR(Y,&rr);DSDPCHKERR(info);
  info=DSDPVecGetC(Y,&cc);DSDPCHKERR(info);
  r=r*rr;
  lbound*=cc; ubound*=cc;
  info=DSDPVecGetSize(vin,&m);DSDPCHKERR(info);
  for (i=1;i<m-1;i++){
    info=DSDPVecGetElement(vrow,i,&ww);DSDPCHKERR(info);
    info=DSDPVecGetElement(vin,i,&vv);DSDPCHKERR(info);
    info=DSDPVecGetElement(Y,i,&yy);DSDPCHKERR(info);
    if (vv==0 || ww==0) continue; 
    sl=(1.0)/COMPUTELBS(lbound,yy,r);
    su=(1.0)/COMPUTEUBS(ubound,yy,r);
    assa=mu*ww*vv*(su*su + sl*sl);
    info=DSDPVecAddElement(vout,i,assa);DSDPCHKERR(info);
  }

  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BoundYConeAddX"
int BoundYConeAddX(LUBounds lucone, double mu, DSDPVec Y, DSDPVec DY, DSDPVec XLU, double *tracexs){
  int i,m,info;
  double sl,su,dsl,dsu,ux,lx;
  double dy,yy,xdots=0,xsum1=0,xsum2=0;
  double r,dr,rr,drr,cr;
  double lbound, ubound;

  DSDPFunctionBegin;
  LUConeValid(lucone);
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  info=DSDPVecGetSize(Y,&m);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&rr);DSDPCHKERR(info);
  info=DSDPVecGetC(Y,&cr);DSDPCHKERR(info);
  info=DSDPVecGetR(DY,&drr);DSDPCHKERR(info);
  r=rr*lucone->r; dr=drr*lucone->r;
  lbound=cr*lucone->lbound; ubound=cr*lucone->ubound;
  mu*=lucone->muscale;
  for (i=1;i<m-1;i++){
    info=DSDPVecGetElement(DY,i,&dy);DSDPCHKERR(info);
    info=DSDPVecGetElement(Y,i,&yy);DSDPCHKERR(info);
    sl=1.0/COMPUTELBS(lbound,yy,r);
    su=1.0/COMPUTEUBS(ubound,yy,r);
    dsl=COMPUTELBS(0,dy,dr);
    dsu=COMPUTEUBS(0,dy,dr);
    lx=mu*(sl-sl*dsl*sl);
    ux=mu*(su-su*dsu*su);
    info=DSDPVecAddElement(XLU,i,ux-lx);DSDPCHKERR(info);
    /*   printf("%d) %4.4f %4.4f %4.4f ",i,ux,lx,ux-lx); */
    xsum1+=lx;
    xsum2+=ux;
    xdots+=lx/sl + ux/su;
  }
  info=DSDPVecAddC(XLU,ubound*xsum1-lbound*xsum2);DSDPCHKERR(info);
  info=DSDPVecAddR(XLU,xsum1+xsum2);DSDPCHKERR(info);
  *tracexs+=xdots;
  /*
  DSDPLogInfo(0,3,"YBounds: Minimum x: %4.4e.  Maximum x: %4.4e, Objective: %4.4e, Trace(xs): %4.4e, Max Ratio of x and s: %4.4e\n",minx,maxx,ppobj,xdots,ratioxs);
  */
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BoundYConeAddS"
int BoundYConeAddS(LUBounds lucone, DSDPVec SL, DSDPVec SU){
  DSDPFunctionBegin;
  LUConeValid(lucone);
  DSDPFunctionReturn(1);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsS"
static int LUBoundsS(void* dcone,DSDPVec Y, DSDPDualFactorMatrix flag, 
		     DSDPTruth *psdefinite){
  int i,m,info;
  LUBounds lucone=(LUBounds)dcone;
  double cr,rr,r;
  double yy,sl,su;
  double lbound, ubound;
  DSDPSchurMat M={0,0,0};

  DSDPFunctionBegin;
  LUConeValid(lucone);
  *psdefinite=DSDP_TRUE;
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  if (lucone->setup==0){
    info=LUBoundsSetup2(dcone,Y,M);DSDPCHKERR(info);
  }
  info=DSDPVecGetC(Y,&cr);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&rr);DSDPCHKERR(info);
  lbound=cr*lucone->lbound; ubound=cr*lucone->ubound;
  r=rr*lucone->r;
  *psdefinite=DSDP_TRUE;
  if (flag==DUAL_FACTOR){
    info=DSDPVecCopy(Y,lucone->YD);DSDPCHKERR(info);
  } else {
    info=DSDPVecCopy(Y,lucone->YP);DSDPCHKERR(info);
  }

  info=DSDPVecGetSize(Y,&m);DSDPCHKERR(info);
  for (i=1;i<m-1;i++){
    info=DSDPVecGetElement(Y,i,&yy);DSDPCHKERR(info);
    sl=COMPUTELBS(lbound,yy,r);
    su=COMPUTEUBS(ubound,yy,r);
    if (sl<=0 || su<=0){ *psdefinite=DSDP_FALSE; break;}
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUInvertS"
static int LUInvertS(void* dcone){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsSetX"
static int LUBoundsSetX(void* dcone,double mu, DSDPVec Y,DSDPVec DY){
  LUBounds lucone=(LUBounds)dcone;
  int info;
  DSDPFunctionBegin;
  LUConeValid(lucone);
  info=DSDPVecCopy(Y,lucone->YP);DSDPCHKERR(info);
  info=DSDPVecCopy(DY,lucone->DYD);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsX"
static int LUBoundsX(void* dcone,double mu, DSDPVec Y,DSDPVec DY, DSDPVec AX , double *tracexs){
  LUBounds lucone=(LUBounds)dcone;
  int info,invisible=lucone->invisible;

  DSDPFunctionBegin;
  LUConeValid(lucone);
  info=LUBoundsSetX(dcone,mu,Y,DY);DSDPCHKERR(info);
  if (!invisible){
    info=BoundYConeAddX(lucone,mu,Y,DY,AX,tracexs);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "LUBoundsComputeMaxStepLength"
static int LUBoundsComputeMaxStepLength(void* dcone, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
  int i,m,info;
  double mstep=1.0e200;
  LUBounds lucone=(LUBounds)dcone;
  double lbound=lucone->lbound, ubound=lucone->ubound;
  double cc,rr,ddr,r,dsl,dsu,sl,su,yy,dy,dr;
  DSDPVec Y;
  DSDPFunctionBegin;

  LUConeValid(lucone);
  *maxsteplength=mstep;
  if (flag==PRIMAL_FACTOR){
    info=DSDPVecCopy(DY,lucone->DYD);DSDPCHKERR(info);
  }
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  info=DSDPVecGetR(DY,&ddr);DSDPCHKERR(info);
  dr=ddr*lucone->r;
  if (flag==DUAL_FACTOR){
    Y=lucone->YD;
  } else {
    Y=lucone->YP;
  }
  info=DSDPVecGetC(Y,&cc);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&rr);DSDPCHKERR(info);
  lbound*=cc; ubound*=cc;
  r=rr*lucone->r;

  info=DSDPVecGetSize(Y,&m);DSDPCHKERR(info);
  for (i=1;i<m-1;i++){
    info=DSDPVecGetElement(DY,i,&dy);DSDPCHKERR(info);
    info=DSDPVecGetElement(Y,i,&yy);DSDPCHKERR(info);
    sl=COMPUTELBS(lbound,yy,r);
    su=COMPUTEUBS(ubound,yy,r);
    dsl=COMPUTELBS(0,dy,dr);
    dsu=COMPUTEUBS(0,dy,dr);

    if (dsl<0){mstep=DSDPMin(mstep,-sl/dsl);}
    if (dsu<0){mstep=DSDPMin(mstep,-su/dsu);}
  }
  *maxsteplength=mstep;
  DSDPLogInfo(0,8,"YBounds: max step: %4.4e\n",mstep);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LUBoundsPotential"
static int LUBoundsPotential(void* dcone, double *logobj, double *logdet){
  LUBounds lucone=(LUBounds)dcone;
  int i,m,info;
  double sl,su,yy,sumlog=0;
  double cc,rr,r;
  double lbound=lucone->lbound, ubound=lucone->ubound;
  DSDPVec Y=lucone->YD;
  DSDPFunctionBegin;

  LUConeValid(lucone);
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  info=DSDPVecGetC(Y,&cc);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&rr);DSDPCHKERR(info);
  lbound*=cc; ubound*=cc;
  r=rr*lucone->r;

  info=DSDPVecGetSize(Y,&m);DSDPCHKERR(info); 
  for (i=1;i<m-1;i++){
    info=DSDPVecGetElement(Y,i,&yy);DSDPCHKERR(info);
    sl=COMPUTELBS(lbound,yy,r);
    su=COMPUTEUBS(ubound,yy,r);
    sumlog+= log(su*sl);
  }
  *logdet=lucone->muscale*sumlog;
  *logobj=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsSparsity"
static int LUBoundsSparsity(void *dcone,int row, int *tnnz, int rnnz[], int m){
  LUBounds lucone=(LUBounds)dcone;
  DSDPFunctionBegin;
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  *tnnz=1;rnnz[row]++;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LPANorm2"
static int LPANorm2( void *dcone, DSDPVec ANorm){
  LUBounds lucone=(LUBounds)dcone;
  double yy,cnorm2;
  int i,m,info;
  DSDPFunctionBegin;
  LUConeValid(lucone);
  if (lucone->invisible){DSDPFunctionReturn(0);}
  info=DSDPVecGetSize(ANorm,&m);DSDPCHKERR(info);
  for (i=1;i<m-1;i++){
    yy=2.0;
    info=DSDPVecAddElement(ANorm,i,yy);
  }
  cnorm2=m*lucone->lbound*lucone->lbound + m*lucone->ubound*lucone->ubound;
  cnorm2=m+1.0;
  info=DSDPVecAddC(ANorm,cnorm2);DSDPCHKERR(info);
  info=DSDPVecAddR(ANorm,2*lucone->r);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LUBoundsView"
int LUBoundsView(LUBounds lucone){
  double lbound=lucone->lbound, ubound=lucone->ubound;
  DSDPFunctionBegin;
  LUConeValid(lucone);
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}
  printf("Lower Bounds for all y variables: %4.8e\n",lbound);
  printf("Upper Bounds for all y variables: %4.8e\n",ubound);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LUBoundsRHS"
static int LUBoundsRHS( void *dcone, double mu, DSDPVec vrow, DSDPVec vrhs1, DSDPVec vrhs2){
  int info,i,m;
  LUBounds lucone=(LUBounds)dcone;
  double as,sl,su;
  double dd,yy,rs=0;
  double cc,rr,r,r0=lucone->r;
  double lbound, ubound;
  DSDPVec Y=lucone->YD;

  DSDPFunctionBegin;
  if (lucone->skipit==DSDP_TRUE){DSDPFunctionReturn(0);}

  LUConeValid(lucone);
  mu*=lucone->muscale;
  info=DSDPVecGetSize(vrow,&m);DSDPCHKERR(info);
  info=DSDPVecGetC(Y,&cc);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&rr);DSDPCHKERR(info);
  lbound=cc*lucone->lbound;
  ubound=cc*lucone->ubound;
  r=rr*lucone->r;
  for (i=1;i<m-1;i++){
    info=DSDPVecGetElement(vrow,i,&dd);DSDPCHKERR(info);
    info=DSDPVecGetElement(Y,i,&yy);DSDPCHKERR(info);
    sl=1.0/COMPUTELBS(lbound,yy,r);
    su=1.0/COMPUTEUBS(ubound,yy,r);
    as=mu*(su-sl);
    if (r){
      rs+=(su+sl);
    }
    if (dd==0) continue;
    info=DSDPVecAddElement(vrhs2,i,dd*as);DSDPCHKERR(info);
  }
  info=DSDPVecAddR(vrhs2,r0*mu*rs);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LUBoundsMonitor"
static int LUBoundsMonitor(void* dcone,int tag){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


static struct  DSDPCone_Ops kops;
static const char *luconename="Bound Y Cone";

#undef __FUNCT__
#define __FUNCT__ "LUBoundsOperationsInitialize"
static int LUBoundsOperationsInitialize(struct  DSDPCone_Ops* coneops){
  int info;
  if (coneops==NULL) return 0;
  info=DSDPConeOpsInitialize(coneops); DSDPCHKERR(info);
  coneops->conehessian=LUBoundsHessian;
  coneops->conesetup=LUBoundsSetup;
  coneops->conesetup2=LUBoundsSetup2;
  coneops->conedestroy=LUBoundsDestroy;
  coneops->conemonitor=LUBoundsMonitor;
  coneops->conecomputes=LUBoundsS;
  coneops->coneinverts=LUInvertS;
  coneops->conecomputex=LUBoundsX;
  coneops->conesetxmaker=LUBoundsSetX;
  coneops->conemaxsteplength=LUBoundsComputeMaxStepLength;
  coneops->conerhs=LUBoundsRHS;
  coneops->conelogpotential=LUBoundsPotential;
  coneops->conesize=LUBoundsSize;
  coneops->conesparsity=LUBoundsSparsity;
  coneops->conehmultiplyadd=LUBoundsMultiply;
  coneops->coneanorm2=LPANorm2;
  coneops->id=12;
  coneops->name=luconename;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BoundYConeSetBounds"
/*!
\fn int BoundYConeSetBounds(LUBounds lucone, double lb, double ub);
\brief Set bounds on the variables
\param lucone cone of bounds.
\param lb lower bound on variables.
\param ub upper bound
*/
int BoundYConeSetBounds(LUBounds lucone, double lb, double ub){
  DSDPFunctionBegin;
  LUConeValid(lucone);
  lucone->lbound=lb;
  lucone->ubound=ub;
  if (lb==0 && ub==0){lucone->skipit=DSDP_TRUE;}
  else { lucone->skipit=DSDP_FALSE;}
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "BoundYConeGetBounds"
/*!
\fn int BoundYConeGetBounds(LUBounds lucone, double *lb, double *ub);
\brief Get bounds on the variables
\param lucone cone of bounds.
\param lb lower bound on variables.
\param ub upper bound
*/
int BoundYConeGetBounds(LUBounds lucone, double *lb, double *ub){
  DSDPFunctionBegin;
  LUConeValid(lucone);
  *lb=lucone->lbound;
  *ub=lucone->ubound;
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPAddLUBounds"
/*!
\fn int DSDPAddLUBounds(DSDP dsdp,LUBounds lucone);
\brief Set the constraints to the solver.
\param lucone cone of bounds.
\param dsdp the solver
*/
int DSDPAddLUBounds(DSDP dsdp,LUBounds lucone){
  int info;
  DSDPFunctionBegin;
  LUConeValid(lucone);
  info=LUBoundsOperationsInitialize(&kops); DSDPCHKERR(info);
  info=DSDPAddCone(dsdp,&kops,(void*)lucone); DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateLUBoundsCone"
/*!
\fn int DSDPCreateLUBoundsCone(DSDP dsdp, LUBounds *dspcone);
\brief Create bounds cone.
\param dsdp the solver
\param dspcone cone of bounds.
*/
int DSDPCreateLUBoundsCone(DSDP dsdp, LUBounds *dspcone){
  int m,info;
  struct LUBounds_C *lucone;
  DSDPFunctionBegin;
  if (!dsdp){DSDPFunctionReturn(1);}
  DSDPCALLOC1(&lucone,struct LUBounds_C,&info);DSDPCHKERR(info);
  *dspcone=lucone;
  lucone->keyid=LUKEY;
  info=DSDPAddLUBounds(dsdp,lucone);DSDPCHKERR(info);
  info=DSDPGetNumberOfVariables(dsdp,&m);DSDPCHKERR(info);
  lucone->muscale=1.0;
  lucone->r=0.0;
  lucone->skipit=DSDP_FALSE;
  lucone->pobj=0; lucone->pax=0; lucone->sumx=0; lucone->xdots=0;
  info=BoundYConeSetBounds(lucone,-1000000,1000000);DSDPCHKERR(info);
  lucone->invisible=1;
  lucone->setup=0;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LUBoundsScaleBarrier"
int LUBoundsScaleBarrier(LUBounds lucone,double muscale){
  DSDPFunctionBegin;
  LUConeValid(lucone);
  if (muscale>0){
    lucone->muscale=muscale;
  }
  DSDPFunctionReturn(0); 
}

