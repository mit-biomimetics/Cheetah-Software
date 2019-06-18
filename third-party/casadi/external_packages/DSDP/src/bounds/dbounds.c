#include "dsdpcone_impl.h"
#include "dsdpsys.h"
#include "dsdp5.h"

#define COMPUTEUBS(a,b,c)  (-(a)-(b)-(c))
/*!
\file dbounds.c
\brief Individually bound variables y.
 */
/*!
struct BCone_C

Internal Structure of bounds on variables.
 */
struct BCone_C{
  int keyid;
  int nn,nnmax;
  int *ib;
  double *u,*au,*us,*uss,*ux,*uds;
  double r,muscale;
  int m;
  double *xuout;
  DSDPVec WY,WY2;
};

#define BKEY  5432
#define BConeValid(a) {if (!(a)||((a)->keyid!=BKEY)){ DSDPSETERR(101,"DSDPERROR: Invalid Bcone object\n");}}

#undef __FUNCT__  
#define __FUNCT__ "BConeSetUp"
static int BConeSetup(void *dcone,DSDPVec y){
  BCone bcone=(BCone)dcone;
  int nn=bcone->nn;
  int info;

  DSDPFunctionBegin;
  if (bcone->nn<1) return 0;
  DSDPCALLOC2(&bcone->us,double,nn,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&bcone->uss,double,nn,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&bcone->ux,double,nn,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&bcone->uds,double,nn,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeSetUp2"
static int BConeSetup2(void *dcone, DSDPVec Y, DSDPSchurMat M){
  int info;
  BCone bcone=(BCone)dcone;
  DSDPFunctionBegin;
  info=DSDPVecDuplicate(Y,&bcone->WY);DSDPCHKERR(info);
  info=DSDPVecDuplicate(Y,&bcone->WY2);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeDestroy"
static int BConeDestroy(void *dcone){
  int info;
  BCone bcone=(BCone)dcone;
  DSDPFunctionBegin;
  DSDPFREE(&bcone->ib,&info);DSDPCHKERR(info);
  DSDPFREE(&bcone->u,&info);DSDPCHKERR(info);
  DSDPFREE(&bcone->au,&info);DSDPCHKERR(info);
  DSDPFREE(&bcone->us,&info);DSDPCHKERR(info);
  DSDPFREE(&bcone->uss,&info);DSDPCHKERR(info);
  DSDPFREE(&bcone->uds,&info);DSDPCHKERR(info);
  DSDPFREE(&bcone->ux,&info);DSDPCHKERR(info);

  info=DSDPVecDestroy(&bcone->WY);DSDPCHKERR(info);
  info=DSDPVecDestroy(&bcone->WY2);DSDPCHKERR(info);

  DSDPFREE(&dcone,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeSize"
static int BConeSize(void *dcone, double *n){
  BCone bcone=(BCone)dcone;
  DSDPFunctionBegin;
  *n=(double)(bcone->nn);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "BConeComputeS"
static int BConeComputeS(BCone bcone,DSDPVec Y,double ss[], int n){
  int i,ii,info;
  int *ib=bcone->ib, nn=bcone->nn;
  double cr,rr,r,yy,*au=bcone->au,*u=bcone->u;

  DSDPFunctionBegin;
  info=DSDPVecGetC(Y,&cr);
  info=DSDPVecGetR(Y,&rr);
  r=rr*bcone->r;
  for (i=0;i<nn;i++){
    ii=ib[i];
    info=DSDPVecGetElement(Y,ii,&yy);
    ss[i]=-cr*u[i]-au[i]*yy-r;
  }  
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeHessian"
static int BConeHessian(void* dcone, double mu, DSDPSchurMat M,
			DSDPVec vrhs1, DSDPVec vrhs2){
  int info,i,ii;
  BCone bcone=(BCone)dcone;
  int *ib=bcone->ib, nn=bcone->nn;
  double *us=bcone->us;
  double *au=bcone->au,*u=bcone->u;
  double dd,cc,rr,cs,as,rs;
  double r=bcone->r;
  DSDPVec DD=bcone->WY,MScale=bcone->WY2;
  
  DSDPFunctionBegin;
  if (bcone->nn<1) return 0;
  mu*=bcone->muscale;
  info=DSDPVecZero(DD);DSDPCHKERR(info);

  info=DSDPSchurMatDiagonalScaling(M,MScale);DSDPCHKERR(info);
  info=DSDPVecGetC(MScale,&cc);DSDPCHKERR(info);
  info=DSDPVecGetR(MScale,&rr);DSDPCHKERR(info);
  for (i=0;i<nn;i++){
    ii=ib[i];
    info=DSDPVecGetElement(MScale,ii,&dd);DSDPCHKERR(info);

    cs=cc*u[i]/us[i];
    as=dd*au[i]/us[i];
    rs=rr*r/us[i];

    if (cs){
      info=DSDPVecAddC(vrhs2,mu*cs);DSDPCHKERR(info);
      info=DSDPVecAddC(DD,mu*cs*cs);DSDPCHKERR(info);
      info=DSDPSchurMatAddC(M,ii,mu*as*cs);DSDPCHKERR(info);
      info=DSDPSchurMatAddR(M,0,mu*cs*rs);DSDPCHKERR(info);
    }
    if (as){
      info=DSDPVecAddElement(vrhs2,ii,mu*as);DSDPCHKERR(info);
      info=DSDPVecAddElement(DD,ii,mu*as*as);DSDPCHKERR(info);
    }
    if (rs){
      info=DSDPVecAddR(vrhs2,mu*rs);DSDPCHKERR(info);
      info=DSDPVecAddR(DD,mu*rs*rs);DSDPCHKERR(info);
      info=DSDPSchurMatAddR(M,ii,mu*as*rs);DSDPCHKERR(info);
    }

  }
  info=DSDPSchurMatAddDiagonal(M,DD);DSDPCHKERR(info);
  
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeHessian"
static int BConeRHS(void* dcone, double mu, DSDPVec vrow,
			 DSDPVec vrhs1, DSDPVec vrhs2){
  int info,i,ii;
  BCone bcone=(BCone)dcone;
  int *ib=bcone->ib, nn=bcone->nn;
  double *us=bcone->us, *au=bcone->au,*u=bcone->u;
  double dd,cc,rr,cs,as,rs;
  double r=bcone->r;
  
  DSDPFunctionBegin;
  mu*=bcone->muscale;
  info=DSDPVecGetC(vrow,&cc);DSDPCHKERR(info);
  info=DSDPVecGetR(vrow,&rr);DSDPCHKERR(info);
  for (i=0;i<nn;i++){
    ii=ib[i];
    info=DSDPVecGetElement(vrow,ii,&dd);DSDPCHKERR(info);

    cs=cc*u[i]/us[i];
    as=dd*au[i]/us[i];
    rs=rr*r/us[i];

    if (cs){
      info=DSDPVecAddC(vrhs2,mu*cs);DSDPCHKERR(info);
    }
    if (as){
      info=DSDPVecAddElement(vrhs2,ii,mu*as);DSDPCHKERR(info);
    }
    if (rs){
      info=DSDPVecAddR(vrhs2,mu*rs);DSDPCHKERR(info);
    }

  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "BConeMultiply"
static int BConeMultiply(void* dcone,  double mu, DSDPVec vrow, DSDPVec vin, DSDPVec vout){
  int info,i,ii;
  BCone bcone=(BCone)dcone;
  int *ib=bcone->ib, nn=bcone->nn;
  double *au=bcone->au,*us=bcone->us;
  double assa,dd,vv;
   
  DSDPFunctionBegin;
  mu*=bcone->muscale;
  for (i=0;i<nn;i++){
    ii=ib[i];
    info=DSDPVecGetElement(vin,ii,&dd);DSDPCHKERR(info);
    info=DSDPVecGetElement(vrow,ii,&vv);DSDPCHKERR(info);
    if (dd==0 || vv==0) continue;
    assa=(au[i]/us[i]);
    assa=mu*vv*assa*assa;
    info=DSDPVecAddElement(vout,ii,assa);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeS"
static int BConeS(void* dcone,DSDPVec Y,DSDPDualFactorMatrix flag,DSDPTruth *psdefinite){
  int i,info;
  BCone bcone=(BCone)dcone;
  int nn=bcone->nn;
  double *us;
  DSDPFunctionBegin;

  if (flag==DUAL_FACTOR){
    us=bcone->us;
  } else {
    us=bcone->uss;
  }
  info=BConeComputeS(bcone,Y,us,nn);DSDPCHKERR(info);
  *psdefinite=DSDP_TRUE;
  for (i=0;i<nn;i++){
    if (us[i]<=0){*psdefinite=DSDP_FALSE;break;}
  }

  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeSInvert"
static int BConeSInvert(void* dcone){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "BConeSetX"
static int BConeSetX(void* dcone,double mu, DSDPVec Y,DSDPVec DY){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeX"
static int BConeX(void* dcone,double mu, DSDPVec Y,DSDPVec DY, DSDPVec AX , double *tracexs){
  int i,ii,info;
  BCone bcone=(BCone)dcone;
  int *ib=bcone->ib, nn=bcone->nn;
  double *au=bcone->au, *us=bcone->uss, *ux=bcone->ux, *uds=bcone->uds, *u=bcone->u;
  double *xuout=bcone->xuout;
  double ds,dus,dau,xx,cr,rr;
  double pobj1=0,xdots1=0;
  DSDPTruth psdefinite;

  DSDPFunctionBegin;

  info=BConeS(dcone,Y,PRIMAL_FACTOR,&psdefinite); DSDPCHKERR(info);
  info=BConeComputeS(bcone,DY,uds,nn);DSDPCHKERR(info);
  info=DSDPVecGetC(Y,&cr);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&rr);DSDPCHKERR(info);

  mu*=bcone->muscale;
  for (i=0;i<nn;i++){
    ii=ib[i];
    ds=uds[i]; dus=us[i]; dau=au[i];
    xx=(mu/dus)-(mu/dus)*(ds/dus);
    ux[i]=xx;
    info=DSDPVecAddElement(AX,ii,dau*xx);DSDPCHKERR(info);
    xdots1+=xx*us[i];
    pobj1+=xx*u[i];
    if (xuout) xuout[i]=xx;
  }

  info=DSDPVecAddC(AX,pobj1);DSDPCHKERR(info);
  *tracexs+=xdots1;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "BConeComputeMaxStepLength"
static int BConeComputeMaxStepLength(void* dcone, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
  int i,info;
  double mstep=1.0e200;
  BCone bcone=(BCone)dcone;
  int nn=bcone->nn;
  double *us,*uds=bcone->uds;
  DSDPFunctionBegin;
  
  if (bcone->nn<1) return 0;
  if (flag==DUAL_FACTOR){
    us=bcone->us;
  } else {
    us=bcone->uss;
  }

  info=BConeComputeS(bcone,DY,uds,nn);DSDPCHKERR(info);

  for (i=0;i<nn;i++){
    if (uds[i]<0){mstep=DSDPMin(mstep,-us[i]/uds[i]);}
  }

  *maxsteplength=mstep;

  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "BConePotential"
static int BConePotential(void* dcone, double *logobj, double *logdet){
  BCone bcone=(BCone)dcone;
  int i;
  double sumlog=0;
  double mu=bcone->muscale;
  int nn=bcone->nn;
  double *us=bcone->us;

  DSDPFunctionBegin;
  if (bcone->nn<1) return 0;
  for (i=0;i<nn;i++){
    sumlog+= mu*log(us[i]);
  }
  *logdet=sumlog;
  *logobj=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BConeSparsity"
static int BConeSparsity(void *dcone,int row, int *tnnz, int rnnz[], int m){
  DSDPFunctionBegin;
  *tnnz=1;rnnz[row]++;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "BConeMonitor"
static int BConeMonitor(void *dcone,int di){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LPANorm2"
static int LPANorm2( void *dcone, DSDPVec ANorm){
  BCone bcone=(BCone)dcone;
  int *ib=bcone->ib, nn=bcone->nn;
  double *u=bcone->u;
  double yy=1.0,cnorm2=0;
  int i,ii,info;
  DSDPFunctionBegin;
  for (i=0;i<nn;i++){
    ii=ib[i]; yy=1.0;
    info=DSDPVecAddElement(ANorm,ii,yy);
    yy=u[i]; cnorm2+=yy*yy;
  }
  info=DSDPVecAddC(ANorm,cnorm2);
  info=DSDPVecAddR(ANorm,bcone->r*bcone->nn);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "BConeView"
int BConeView(BCone bcone){
  int i,ii;
  int *ib=bcone->ib, nn=bcone->nn;
  double  *au=bcone->au,*u=bcone->u;
  DSDPFunctionBegin;
  BConeValid(bcone);
  ib=bcone->ib; nn=bcone->nn;
  au=bcone->au;  u=bcone->u;
  for (i=0;i<nn;i++){
    ii=ib[i];
    if (au[i]>0){
      printf("Upper Bound.  Var %d: %4.8e\n",ii,u[i]);
    } else {
      printf("Lower Bound.  Var %d: %4.8e\n",ii,u[i]);
    }
  }
  DSDPFunctionReturn(0);
}


static struct  DSDPCone_Ops kops;
static const char *bconename="VariableBounds Cone";

#undef __FUNCT__
#define __FUNCT__ "BConeOperationsInitialize"
static int BConeOperationsInitialize(struct  DSDPCone_Ops* coneops){
  int info;
  if (coneops==NULL) return 0;
  info=DSDPConeOpsInitialize(coneops); DSDPCHKERR(info);
  coneops->conehessian=BConeHessian;
  coneops->conerhs=BConeRHS;
  coneops->conesetup=BConeSetup;
  coneops->conesetup2=BConeSetup2;
  coneops->conedestroy=BConeDestroy;
  coneops->conecomputes=BConeS;
  coneops->coneinverts=BConeSInvert;
  coneops->conecomputex=BConeX;
  coneops->conesetxmaker=BConeSetX;
  coneops->conemaxsteplength=BConeComputeMaxStepLength;
  coneops->conelogpotential=BConePotential;
  coneops->conesize=BConeSize;
  coneops->conemonitor=BConeMonitor;
  coneops->conesparsity=BConeSparsity;
  coneops->conehmultiplyadd=BConeMultiply;
  coneops->coneanorm2=LPANorm2;
  coneops->id=2;
  coneops->name=bconename;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPAddBounds"
int DSDPAddBounds(DSDP dsdp,BCone bcone){
  int info;
  DSDPFunctionBegin;
  BConeValid(bcone);
  info=BConeOperationsInitialize(&kops); DSDPCHKERR(info);
  info=DSDPAddCone(dsdp,&kops,(void*)bcone); DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateBCone"
/*!
\fn int DSDPCreateBCone(DSDP dsdp, BCone *dspcone);
\brief Create a new cone that represents bounds on the y variables
\ingroup Bounds
\param dsdp the solver
\param dspcone new cone

For example, to bound the first of three y variables below by 0 and 10, 
\code
DSDP dsdp;
BCone bcone;
int xl[3],xu[3];
DSDPCreate(3,&dsdp);
DSDPCreateBCone(dsdp,&bcone);
BConeSetLowerBound(bcone,1,0.0);
BConeSetLowerBound(bcone,1,10.0);
...

DSDPComputeX(dsdp);
BConeCopyX(bcone,xl,xu,3);
printf("The sensitivity of the solution to lower bound is %4.4f\n",xl[0]);
\endcode
*/
int DSDPCreateBCone(DSDP dsdp, BCone *dspcone){
  int m,info;
  struct BCone_C *bcone;
  DSDPFunctionBegin;
  if (!dsdp){DSDPFunctionReturn(1);}
  DSDPCALLOC1(&bcone,struct BCone_C,&info);DSDPCHKERR(info);
  bcone->keyid=BKEY;
  *dspcone=bcone;
  /* info=DSDPAddBounds(dsdp,bcone);DSDPCHKERR(info); */
  info=BConeOperationsInitialize(&kops); DSDPCHKERR(info);
  info=DSDPAddCone(dsdp,&kops,(void*)bcone); DSDPCHKERR(info);
  info=DSDPGetNumberOfVariables(dsdp,&m);DSDPCHKERR(info);
  bcone->nn=0;
  bcone->m=m;
  bcone->muscale=1.0;
  bcone->r=1.0;
  bcone->nnmax=0;
  bcone->xuout=0;
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "BConeScaleBarrier"
int BConeScaleBarrier(BCone bcone,double muscale){
  DSDPFunctionBegin;
  BConeValid(bcone);
  if (muscale>0){
    bcone->muscale=muscale;
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "BConeCopyX"
/*!
\fn int BConeCopyX(BCone bcone,double xl[], double xu[], int m);
\brief Copy the variables into arrays
\param bcone Bounds
\param xl array
\param xu array
\param m length of the arrays
\ingroup Bounds
\sa DSDPComputeX()

This routine will set the values of this array to the value of the
corresponding variable.  When no bound is present, the variable will
equal zero.

*/
int BConeCopyX(BCone bcone,double xl[], double xu[], int m){
  int i,ii,*ib,nn;
  double *xx,*au;
  DSDPFunctionBegin;
  BConeValid(bcone);
  if (m!=bcone->m){ DSDPSETERR1(6,"Invalid Array Length.\n",bcone->m);}
  xx=bcone->ux; au=bcone->au; nn=bcone->nn; ib=bcone->ib;
  for (i=0;i<m;i++){
    xl[i]=0;xu[i]=0;
  }
  for (i=0;i<nn;i++){
    ii=ib[i]-1;
    if (au[i]<0){
      xl[ii]+=xx[i];
    } else {
      xu[ii]+=xx[i];
    }
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "BConeCopyXSingle"
/*!
\fn int BConeCopyXSingle(BCone bcone,double lambda[], int m);
\brief Copy the variables into arrays
\param bcone Bounds
\param lambda array
\param m length of the arrays
\ingroup Bounds
\sa DSDPComputeX()

This routine will set the values of this array to the value of the
corresponding variable.  When no bound is present, the variable will
equal zero.

*/
int BConeCopyXSingle(BCone bcone,double lambda[], int m){
  int i,ii,*ib,nn;
  double *xx,*au;
  DSDPFunctionBegin;
  BConeValid(bcone);
  if (m!=bcone->m){ DSDPSETERR1(6,"Invalid Array Length.\n",bcone->m);}
  xx=bcone->ux; au=bcone->au; nn=bcone->nn; ib=bcone->ib;
  for (i=0;i<m;i++){
    lambda[i]=0;
  }
  for (i=0;i<nn;i++){
    ii=ib[i]-1;
    lambda[ii] += ((au[i]<0) ? -1: 1) * xx[i];
  }
  DSDPFunctionReturn(0); 
}
#undef __FUNCT__
#define __FUNCT__ "BConeSetBound"
int BConeSetBound(BCone bcone,int vari, double ai, double bound){
  int spot,info;
  DSDPFunctionBegin;
  BConeValid(bcone);
  if (vari<1 || vari>bcone->m){ DSDPSETERR2(6,"Invalid Variable number 1 <= %d <= %d.\n",vari,bcone->m);}
  if (bcone->nn>=bcone->nnmax){
    DSDPLogInfo(0,19,"REALLOCATING SPACE FOR BOUNDS! %d \n",bcone->nnmax);
    info=BConeAllocateBounds(bcone,2*bcone->nn+4);DSDPCHKERR(info);
  }
  spot=bcone->nn;
  bcone->u[spot]=bound;
  bcone->au[spot]=ai;
  bcone->ib[spot]=vari;
  bcone->nn++;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "BConeSetLowerBound"
/*!
\fn int BConeSetLowerBound(BCone bcone,int vari, double lbound);
\brief Set a lower bound on a variable y. 
\ingroup Bounds
\param bcone bounds
\param vari y variable number
\param lbound lower bound
*/
int BConeSetLowerBound(BCone bcone,int vari, double lbound){
  int info;
  DSDPFunctionBegin;
  info=BConeSetBound(bcone,vari,-1.0,-lbound);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "BConeSetUpperBound"
/*!
\fn int BConeSetUpperBound(BCone bcone,int vari, double ubound);
\brief Set an upper bound on a variable y. 
\ingroup Bounds
\param bcone bounds
\param vari y variable number
\param ubound upper bound
*/
int BConeSetUpperBound(BCone bcone,int vari, double ubound){
  int info;
  DSDPFunctionBegin;
  info=BConeSetBound(bcone,vari,1.0,ubound);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "BConeSetUnboundedLower"
/*!
\fn int BConeSetUnboundedLower(BCone bcone,int vari);
\brief Set a lower bound on a variable y. 
\ingroup Bounds
\param bcone bounds
\param vari y variable number
\param lbound lower bound
*/
int BConeSetUnboundedLower(BCone bcone,int vari){
  int info;
  DSDPFunctionBegin;
  info=BConeSetBound(bcone,vari,0.0,1.0);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "BConeSetUnboundedUpper"
/*!
\fn int BConeSetUnboundedUpper(BCone bcone,int vari);
\brief Set an upper bound on a variable y. 
\ingroup Bounds
\param bcone bounds
\param vari y variable number
*/
int BConeSetUnboundedUpper(BCone bcone,int vari){
  int info;
  DSDPFunctionBegin;
  info=BConeSetBound(bcone,vari,0.0,1.0);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "BConeSetPSlackVariable"
/*!
\fn int BConeSetPSlackVariable(BCone bcone,int vari);
\brief Set a slack variable to a constraint in (P).
\ingroup Bounds
\param bcone bounds
\param vari y variable number
\sa BConeSetUpperBound()
\sa BConeSetPSurplusVariable()

When a constraint in (P) is best expressed with an inequality
constraint, a slack variable may be added.

This command is equivalent to setting an upper bound on
variable $y_i$ to zero.
*/
int BConeSetPSlackVariable(BCone bcone,int vari){
  int info;
  DSDPFunctionBegin;
  info=BConeSetUpperBound(bcone,vari,0);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "BConeSetPSurplusVariable"
/*!
\fn int BConeSetPSurplusVariable(BCone bcone,int vari);
\brief Set a surplus variable in constraint in (P).
\ingroup Bounds
\param bcone bounds
\param vari y variable number
\sa BConeSetLowerBound()
*/
int BConeSetPSurplusVariable(BCone bcone,int vari){
  int info;
  DSDPFunctionBegin;
  info=BConeSetLowerBound(bcone,vari,0);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "BConeAllocateBounds"
/*!
\fn int BConeAllocateBounds(BCone bcone,int nnz);
\brief Set a surplus variable in constraint in (P).
\ingroup Bounds
\param bcone bounds
\param nnz number of bounds that will be set.

This optional routine allows for efficient allocation
of memory for this object.
\sa BConeSetLowerBound()
\sa BConeSetUpperBound()
*/
int BConeAllocateBounds(BCone bcone, int nnz){
  int j,info,*uindex;
  double *uu,*au;

  DSDPFunctionBegin;
  BConeValid(bcone);
  if (nnz<=bcone->nnmax){DSDPFunctionReturn(0);}
  DSDPCALLOC2(&uu,double,(nnz),&info); DSDPCHKERR(info);
  DSDPCALLOC2(&au,double,(nnz),&info); DSDPCHKERR(info);
  DSDPCALLOC2(&uindex,int,(nnz),&info); DSDPCHKERR(info);
  for (j=0;j<nnz;j++){uu[j]=0; uindex[j]=0; au[j]=0;}
  if (bcone->nnmax>0){
    for (j=0;j<bcone->nn;j++){uu[j]=bcone->u[j];}
    for (j=0;j<bcone->nn;j++){uindex[j]=bcone->ib[j];}
    for (j=0;j<bcone->nn;j++){au[j]=bcone->au[j];}
    DSDPFREE(&bcone->u,&info);DSDPCHKERR(info);
    DSDPFREE(&bcone->au,&info);DSDPCHKERR(info);
    DSDPFREE(&bcone->ib,&info);DSDPCHKERR(info);
  } else {
    bcone->nn=0;
  }
  bcone->nnmax=nnz;
  bcone->u=uu;
  bcone->au=au;
  bcone->ib=uindex;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BConeSetXArray"
int BConeSetXArray(BCone bcone,double *xuout, int n){
  DSDPFunctionBegin;
  BConeValid(bcone);
  if (n==bcone->nn) bcone->xuout=xuout;
  DSDPFunctionReturn(0); 
}

