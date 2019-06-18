
#include "dsdpcone_impl.h"
#include "dsdpsys.h"
#include "dsdp5.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/*!
\file dsdplp.c
\brief Set linear inequalities in (D) and implement the DSDPCone operations
 */
typedef struct {
  int    nrow;
  int    ncol;
  int    owndata;
  
  const double *an;
  const int    *col;
  const int    *nnz;
  int    *nzrows;
  int    nnzrows;

} smatx;

/*!
struct LPCone_C ;
Internal Structure of LP variables.
 */
struct LPCone_C{
  smatx *A,*AT;
  DSDPVec C;
  DSDPVec PS,DS,X;
  double sscale;
  double r;
  double muscale;
  DSDPVec Y,WY,WY2,WX,WX2;
  double *xout;
  int n,m;
};


static int CreateSpRowMatWdata(int,int,const double[], const int[], const int[],smatx **);
static int SpRowMatMult(smatx*,const  double[], int , double[], int);
static int SpRowMatMultTrans(smatx *, const double[],int, double[],int);
static int SpRowMatGetRowVector(smatx*, int, double*,int);
static int SpRowMatGetScaledRowVector(smatx*, int, const double[], double*, int);
static int SpRowMatDestroy(smatx*);
static int SpRowMatView(smatx*);
/*
static int SpRowMatGetSize(smatx *, int *, int *);
static int SpRowMatZero(smatx*);
static int SpRowMatAddRowMultiple(smatx*, int, double, double[], int);
*/
static int SpRowIMultAdd(smatx*,int*,int,int *,int);
static int SpRowMatRowNnz(smatx*, int, int*,int);
static int SpRowMatNorm2(smatx*, int, double*);


#undef __FUNCT__  
#define __FUNCT__ "LPConeSetUp"
static int LPConeSetup(void *dcone,DSDPVec y){
  int m,info;
  LPCone lpcone=(LPCone)dcone;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  m=lpcone->m;
  info=DSDPVecCreateSeq(m+2,&lpcone->WY);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->WY,&lpcone->WY2);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->WY,&lpcone->Y);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->WX);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->WX2);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->PS);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->DS);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->X);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPConeSetUp2"
static int LPConeSetup2(void *dcone, DSDPVec Y, DSDPSchurMat M){
  LPCone lpcone=(LPCone)dcone;
  DSDPFunctionBegin;
  DSDPLogInfo(0,19,"Setup LP Cone of dimension: %d\n",lpcone->n);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LPConeDestroy"
static int LPConeDestroy(void *dcone){
  int info;
  LPCone lpcone=(LPCone)dcone;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  info=DSDPVecDestroy(&lpcone->DS);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lpcone->PS);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lpcone->C);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lpcone->X);DSDPCHKERR(info);
  info=SpRowMatDestroy(lpcone->A);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lpcone->WX2);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lpcone->WY2);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lpcone->WY);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lpcone->Y);DSDPCHKERR(info);
  info=DSDPVecDestroy(&lpcone->WX);DSDPCHKERR(info);
  DSDPFREE(&lpcone,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPConeSize"
static int LPConeSize(void *dcone, double *n){
  LPCone lpcone=(LPCone)dcone;
  DSDPFunctionBegin;
  *n=lpcone->muscale*lpcone->n;
  DSDPFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "LPComputeAX"
static int LPComputeAX( LPCone lpcone,DSDPVec X, DSDPVec Y){
  int info,m=lpcone->m,n=lpcone->n;
  double *x,*y,ppobj;
  smatx *A=lpcone->A;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  info=DSDPVecGetSize(X,&n);DSDPCHKERR(info);
  info=DSDPVecDot(lpcone->C,X,&ppobj);DSDPCHKERR(info);
  info=DSDPVecSetC(Y,ppobj);
  info=DSDPVecSum(X,&ppobj);DSDPCHKERR(info);
  info=DSDPVecSetR(Y,ppobj*lpcone->r);DSDPCHKERR(info);
  info=DSDPVecGetArray(Y,&y);DSDPCHKERR(info);
  info=DSDPVecGetArray(X,&x);DSDPCHKERR(info);
  info=SpRowMatMult(A,x,n,y+1,m); 
  info=DSDPVecRestoreArray(X,&x);DSDPCHKERR(info);
  info=DSDPVecRestoreArray(Y,&y);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPComputeATY"
static int LPComputeATY(LPCone lpcone,DSDPVec Y, DSDPVec S){
  int info,m=lpcone->m,n=lpcone->n;
  double cc,r,*s,*y;
  DSDPVec C=lpcone->C;
  smatx *A=lpcone->A;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  info=DSDPVecGetC(Y,&cc);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&r);DSDPCHKERR(info);
  info=DSDPVecGetSize(S,&n);DSDPCHKERR(info);
  info=DSDPVecGetArray(S,&s);DSDPCHKERR(info);
  info=DSDPVecGetArray(Y,&y);DSDPCHKERR(info);
  info=SpRowMatMultTrans(A,y+1,m,s,n); DSDPCHKERR(info);
  info=DSDPVecRestoreArray(S,&s);DSDPCHKERR(info);
  info=DSDPVecRestoreArray(Y,&y);DSDPCHKERR(info);
  info=DSDPVecAXPY(cc,C,S);DSDPCHKERR(info);
  info=DSDPVecShift(r*lpcone->r,S);DSDPCHKERR(info);
  info=DSDPVecScale(-1.0,S);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "LPConeHessian"
static int LPConeHessian(void* dcone, double mu, DSDPSchurMat M,
			 DSDPVec vrhs1, DSDPVec vrhs2){
  int info,i,m,n,ncols;
  double r=1.0,*wx,*wx2;
  LPCone lpcone=(LPCone)dcone;
  DSDPVec WX=lpcone->WX,WX2=lpcone->WX2,WY=lpcone->WY,WY2=lpcone->WY2,S=lpcone->DS;
  smatx *A=lpcone->A;
  
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  mu*=lpcone->muscale;
  info=DSDPVecGetSize(vrhs1,&m);DSDPCHKERR(info);
  info=DSDPVecGetSize(WX,&n);DSDPCHKERR(info);
  info=DSDPVecSet(mu,WX2);DSDPCHKERR(info);
  info=DSDPVecPointwiseDivide(WX2,S,WX2);DSDPCHKERR(info);
  info=DSDPVecPointwiseDivide(WX2,S,WX2);DSDPCHKERR(info);
  for (i=0;i<m;i++){

    info=DSDPSchurMatRowColumnScaling(M,i,WY2,&ncols); DSDPCHKERR(info); 
    if (ncols==0) continue;

    if (i==0){
      info=DSDPVecPointwiseMult(lpcone->C,WX2,WX); DSDPCHKERR(info); 
    } else if (i==m-1){
      info=DSDPVecScaleCopy(WX2,r,WX); DSDPCHKERR(info); 
    } else {
      info=DSDPVecGetArray(WX,&wx);
      info=DSDPVecGetArray(WX2,&wx2);DSDPCHKERR(info);
      info=SpRowMatGetScaledRowVector(A,i-1,wx2,wx,n);
      info=DSDPVecRestoreArray(WX,&wx);
      info=DSDPVecRestoreArray(WX2,&wx2);
    }

    info=LPComputeAX(lpcone,WX,WY);DSDPCHKERR(info);

    info=DSDPVecPointwiseMult(WY2,WY,WY);DSDPCHKERR(info);

    info=DSDPSchurMatAddRow(M,i,1.0,WY);DSDPCHKERR(info);
  }

  /* Compute AS^{-1}  */
  info=DSDPVecSet(mu,WX);DSDPCHKERR(info);
  info=DSDPVecPointwiseDivide(WX,S,WX);DSDPCHKERR(info);
  info=LPComputeAX(lpcone,WX,WY);DSDPCHKERR(info);

  info=DSDPSchurMatDiagonalScaling(M, WY2);DSDPCHKERR(info);
  info=DSDPVecPointwiseMult(WY2,WY,WY);DSDPCHKERR(info);
  info=DSDPVecAXPY(1.0,WY,vrhs2);DSDPCHKERR(info);

  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPConeHessian"
static int LPConeRHS(void* dcone, double mu, DSDPVec vrow,
			 DSDPVec vrhs1, DSDPVec vrhs2){
  int info;
  LPCone lpcone=(LPCone)dcone;
  DSDPVec WX=lpcone->WX,WY=lpcone->WY,S=lpcone->DS;
  
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  mu*=lpcone->muscale;

  /* Compute AS^{-1}  */
  info=DSDPVecSet(mu,WX);DSDPCHKERR(info);
  info=DSDPVecPointwiseDivide(WX,S,WX);DSDPCHKERR(info);
  info=LPComputeAX(lpcone,WX,WY);DSDPCHKERR(info);

  info=DSDPVecPointwiseMult(vrow,WY,WY);DSDPCHKERR(info);
  info=DSDPVecAXPY(1.0,WY,vrhs2);DSDPCHKERR(info);
  
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPConeMultiply"
static int LPConeMultiply(void* dcone,  double mu, DSDPVec vrow, DSDPVec vin, DSDPVec vout){
  int info;
  LPCone lpcone=(LPCone)dcone;
  DSDPVec WX=lpcone->WX,S=lpcone->DS,WY=lpcone->WY;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  mu*=lpcone->muscale;
  info=LPComputeATY(lpcone,vin,WX);DSDPCHKERR(info);
  info=DSDPVecPointwiseDivide(WX,S,WX);DSDPCHKERR(info);
  info=DSDPVecScale(mu,WX);DSDPCHKERR(info);
  info=DSDPVecPointwiseDivide(WX,S,WX);DSDPCHKERR(info);
  info=LPComputeAX(lpcone,WX,WY);DSDPCHKERR(info);
  info=DSDPVecPointwiseMult(WY,vrow,WY);DSDPCHKERR(info);
  info=DSDPVecAXPY(1.0,WY,vout);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPConeSetX"
static int LPConeSetX(void* dcone,double mu, DSDPVec Y,DSDPVec DY){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPConeX"
static int LPConeX(void* dcone,double mu, DSDPVec Y,DSDPVec DY, 
		   DSDPVec AX , double *tracexs){
  int info;
  double dtracexs;
  LPCone lpcone=(LPCone)dcone;
  DSDPVec S=lpcone->PS,WX=lpcone->WX,X=lpcone->X,DS=lpcone->DS,WY=lpcone->WY;
  double *xx,*xout=lpcone->xout;
  DSDPFunctionBegin;

  if (lpcone->n<1) return 0;
  mu*=lpcone->muscale;
  info=LPComputeATY(lpcone,Y,S);DSDPCHKERR(info);

  info=DSDPVecSet(1.0,WX);
  info=DSDPVecPointwiseDivide(WX,S,WX);DSDPCHKERR(info);

  info=LPComputeATY(lpcone,DY,DS);DSDPCHKERR(info);
  info=DSDPVecPointwiseMult(WX,DS,X);DSDPCHKERR(info);

  info=DSDPVecScale(-mu,WX);DSDPCHKERR(info);
  info=DSDPVecPointwiseMult(WX,X,X);DSDPCHKERR(info);
  info=DSDPVecAXPY(-1.0,WX,X);DSDPCHKERR(info);
  info=DSDPVecGetArray(X,&xx);DSDPCHKERR(info);
  for (info=0;info<lpcone->n;info++){
    if (xx[info]<0) xx[info]=0;
  }
  info=DSDPVecRestoreArray(X,&xx);DSDPCHKERR(info);
  info=LPComputeAX(lpcone,X,WY);DSDPCHKERR(info);
  info=DSDPVecAXPY(1.0,WY,AX);DSDPCHKERR(info);
  info=DSDPVecDot(S,X,&dtracexs);DSDPCHKERR(info);  
  *tracexs+=dtracexs;
  info=DSDPVecGetArray(X,&xx);DSDPCHKERR(info);
  if (xout){
    for (info=0;info<lpcone->n;info++){
      if (xout){ xout[info]=xx[info]; }
    }
  }
  info=DSDPVecRestoreArray(X,&xx);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LPConeS"
static int LPConeS(void* dcone, DSDPVec Y, DSDPDualFactorMatrix flag, 
		   DSDPTruth *psdefinite){
  int i,n,info;
  double *s;
  LPCone lpcone=(LPCone)dcone;
  DSDPVec S;

  DSDPFunctionBegin;
  
  if (lpcone->n<1) return 0;
  if (flag==DUAL_FACTOR){
    S=lpcone->DS;
  } else {
    S=lpcone->PS;
  }

  info=DSDPVecCopy(Y,lpcone->Y);DSDPCHKERR(info);
  info=LPComputeATY(lpcone,Y,S);DSDPCHKERR(info);
  info=DSDPVecGetC(Y,&lpcone->sscale);
  info=DSDPVecGetSize(S,&n);DSDPCHKERR(info);
  info=DSDPVecGetArray(S,&s);DSDPCHKERR(info);
  *psdefinite=DSDP_TRUE;
  for (i=0;i<n;i++){ if (s[i]<=0) *psdefinite=DSDP_FALSE;}
  info=DSDPVecRestoreArray(S,&s);DSDPCHKERR(info);
  
  DSDPFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "LPConeInvertS"
static int LPConeInvertS(void* dcone){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LPConeComputeMaxStepLength"
static int LPConeComputeMaxStepLength(void* dcone, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
  int i,n,info;
  double *s,*ds,mstep=1.0e200;
  LPCone lpcone=(LPCone)dcone;
  DSDPVec S,DS=lpcone->WX;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  if (flag==DUAL_FACTOR){
    S=lpcone->DS;
  } else {
    S=lpcone->PS;
  }

  info=LPComputeATY(lpcone,DY,DS);DSDPCHKERR(info);

  info=DSDPVecGetSize(DS,&n);DSDPCHKERR(info);
  info=DSDPVecGetArray(S,&s);DSDPCHKERR(info);
  info=DSDPVecGetArray(DS,&ds);DSDPCHKERR(info);
  for (i=0;i<n;i++) if (ds[i]<0){mstep=DSDPMin(mstep,-s[i]/ds[i]);}
  info=DSDPVecRestoreArray(S,&s);DSDPCHKERR(info);
  info=DSDPVecRestoreArray(DS,&ds);DSDPCHKERR(info);

  *maxsteplength=mstep;
  
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LPConePotential"
static int LPConePotential(void* dcone, double *logobj, double *logdet){
  int i,n,info;
  double *s,mu,sumlog=0;
  LPCone lpcone=(LPCone)dcone;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  mu=lpcone->muscale;
  info=DSDPVecGetArray(lpcone->DS,&s);DSDPCHKERR(info);
  info=DSDPVecGetSize(lpcone->DS,&n);DSDPCHKERR(info);
  for (i=0;i<n;i++){
    sumlog+= mu*log(s[i]);
  }
  info=DSDPVecRestoreArray(lpcone->DS,&s);DSDPCHKERR(info);
  *logdet=sumlog;
  *logobj=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPConeSparsity"
static int LPConeSparsity(void *dcone,int row, int *tnnz, int rnnz[], int m){
  int info,*wi,n;
  double *wd;
  LPCone lpcone=(LPCone)dcone;
  DSDPVec W=lpcone->WX;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  if (row==0) return 0;
  if (row==m-1){
    return 0;
  }
  info=DSDPVecGetSize(W,&n);DSDPCHKERR(info);
  info=DSDPVecGetArray(W,&wd);DSDPCHKERR(info);
  wi=(int*)wd;
  info=SpRowMatRowNnz(lpcone->A,row-1,wi,n);DSDPCHKERR(info);
  info=SpRowIMultAdd(lpcone->A,wi,n,rnnz+1,m-2);DSDPCHKERR(info);
  info=DSDPVecRestoreArray(W,&wd);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "LPConeMonitor"
static int LPConeMonitor( void *dcone, int tag){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "LPANorm2"
static int LPANorm2( void *dcone, DSDPVec ANorm){
  int i,info;
  double dd;
  LPCone lpcone=(LPCone)dcone;
  DSDPFunctionBegin;
  if (lpcone->n<1) return 0;
  info=DSDPVecNorm22(lpcone->C,&dd);DSDPCHKERR(info);
  info=DSDPVecAddC(ANorm,dd);DSDPCHKERR(info);
  for (i=0;i<lpcone->m;i++){
    info=SpRowMatNorm2(lpcone->A,i,&dd);DSDPCHKERR(info);
    info=DSDPVecAddElement(ANorm,i+1,dd);DSDPCHKERR(info);
  }
  info=DSDPVecAddR(ANorm,1.0);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


static struct  DSDPCone_Ops kops;
static const char *lpconename="LP Cone";

#undef __FUNCT__
#define __FUNCT__ "LPConeOperationsInitialize"
static int LPConeOperationsInitialize(struct  DSDPCone_Ops* coneops){
  int info;
  if (coneops==NULL) return 0;
  info=DSDPConeOpsInitialize(coneops); DSDPCHKERR(info);
  coneops->conehessian=LPConeHessian;
  coneops->conerhs=LPConeRHS;
  coneops->conesetup=LPConeSetup;
  coneops->conesetup2=LPConeSetup2;
  coneops->conedestroy=LPConeDestroy;
  coneops->conecomputes=LPConeS;
  coneops->coneinverts=LPConeInvertS;
  coneops->conesetxmaker=LPConeSetX;
  coneops->conecomputex=LPConeX;
  coneops->conemaxsteplength=LPConeComputeMaxStepLength;
  coneops->conelogpotential=LPConePotential;
  coneops->conesize=LPConeSize;
  coneops->conesparsity=LPConeSparsity;
  coneops->conehmultiplyadd=LPConeMultiply;
  coneops->conemonitor=LPConeMonitor;
  coneops->coneanorm2=LPANorm2;
  coneops->id=2;
  coneops->name=lpconename;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPAddLP"
int DSDPAddLP(DSDP dsdp,LPCone lpcone){
  int info;
  DSDPFunctionBegin;
  info=LPConeOperationsInitialize(&kops); DSDPCHKERR(info);
  info=DSDPAddCone(dsdp,&kops,(void*)lpcone); DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateLPCone"
/*!
\fn int DSDPCreateLPCone(DSDP dsdp, LPCone *dspcone){
\brief Create a new object for linear programs and scalar inequalities.
\ingroup LPRoutines
\param dsdp the solver
\param *dspcone new LP Cone object

\note LP data can be specified using one or more LPCone objects.  Although
the use multiple cones is allowed, efficiency is best when all LP data
is set in the same cone.

\code
DSDP dsdp;
LPCone lpcone;
DSDPCreate(3,&dsdp);
DSDPCreateLPCone(dsdp,&lpcone);
\endcode

\sa DSDPCreate()
 */
int DSDPCreateLPCone(DSDP dsdp, LPCone *dspcone){
  int m,info;
  struct LPCone_C *lpcone;
  DSDPFunctionBegin;
  DSDPCALLOC1(&lpcone,struct LPCone_C,&info);DSDPCHKERR(info);
  *dspcone=lpcone;
  /*
  info=DSDPAddLP(dsdp,lpcone);DSDPCHKERR(info);
  */
  info=LPConeOperationsInitialize(&kops); DSDPCHKERR(info);
  info=DSDPAddCone(dsdp,&kops,(void*)lpcone); DSDPCHKERR(info);
  info=DSDPGetNumberOfVariables(dsdp,&m);DSDPCHKERR(info);
  lpcone->m=m;
  lpcone->muscale=1.0;
  lpcone->n=0;
  lpcone->xout=0;
  lpcone->r=1.0;
  info=DSDPVecCreateSeq(0,&lpcone->C);DSDPCHKERR(info);
  info=DSDPVecCreateSeq(0,&lpcone->WY);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->WX);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->WX2);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->PS);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->DS);DSDPCHKERR(info);
  info=DSDPVecDuplicate(lpcone->C,&lpcone->X);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "LPConeGetXArray"
/*!
\fn int LPConeGetXArray(LPCone lpcone, double *x[],int *n);
\brief Get the array used to store the x variables
\ingroup LPRoutines
\param lpcone LP Cone
\param x array of variables
\param n the dimension of the cone and length of the array
\sa DSDPCreateLPCone()

For example, after calling DSDPComputeX(),
\code
int i,n;
double *x;
LPConeGetXArray(lpcone,&x,&n);
for (i=0;i<n,i++) printf("x[%d] = %16.8f\n",x[i]);
\endcode
*/
int LPConeGetXArray(LPCone lpcone,double *x[], int *n){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecGetArray(lpcone->X,x);DSDPCHKERR(info);
  info=DSDPVecGetSize(lpcone->X,n);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LPConeGetSArray"
/*
\fn int LPConeGetSArray(LPCone lpcone, double *s[], int *n);
\brief Get the array used to store the s variables
\ingroup LPRoutines
\param lpcone LP Cone
\param s array of variables
\param n the dimension of the cone and length of the array
\sa DSDPCreateLPCone()
 */
int LPConeGetSArray(LPCone lpcone,double *s[], int *n){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecGetArray(lpcone->DS,s);DSDPCHKERR(info);
  info=DSDPVecGetSize(lpcone->DS,n);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LPConeCopyS"
/*!
\fn int LPConeCopyS(LPCone lpcone, double s[],int n);
\brief Copy the variables s into the spedified array.
\ingroup LPRoutines
\param lpcone LP Cone
\param s array
\param n the conic dimension and length of the array
\sa DSDPCreateLPCone()
\sa LPConeGetDimension()
 */
int LPConeCopyS(LPCone lpcone,double s[], int n){
  int i,info;
  double *ss,sscale=lpcone->sscale;
  DSDPTruth psdefinite;
  DSDPFunctionBegin;
  info=LPConeS((void*)lpcone,lpcone->Y,DUAL_FACTOR ,&psdefinite);DSDPCHKERR(info);
  info=DSDPVecGetArray(lpcone->DS,&ss);DSDPCHKERR(info);
  for (i=0;i<n;i++) s[i]=ss[i]/fabs(sscale);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LPConeGetDimension"
/*!
\fn int LPConeGetDimension(LPCone lpcone,int *n);
\brief Get the dimension is the number of variables x,
which equals the number of slack variables s.
\param lpcone LP Cone
\param n dimension
\ingroup LPRoutines
*/
int LPConeGetDimension(LPCone lpcone,int *n){
  DSDPFunctionBegin;
  *n=(int)(lpcone->n*lpcone->muscale);
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "LPConeScaleBarrier"
int LPConeScaleBarrier(LPCone lpcone,double muscale){
  DSDPFunctionBegin;
  if (muscale>0){
    lpcone->muscale=muscale;
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LPConeSetData"
/*!
\fn LPConeSetData(LPCone lpcone, int n, const int ik[],const int cols[],const double vals[])
\brief Set data into the LP cone.
\param lpcone the LP cone
\param n the number of inequalities
\param ik the number of nonzeros in each column of the matrix
\param cols array of column numbers
\param vals array of LP data
\ingroup LPRoutines
\sa DSDPSetDualObjective()

For example, consider the following problem in the form of (D):
\f[ \begin{array}{llllll}
\mbox{Maximize} & & y_1 & + & y_2 \\
\mbox{Subject to}
 & & 4 y_1 & + & 2 y_2 & \leq 6 \\
 & & 3 y_1 & + & 7 y_2 & \leq 10 \\
 & &   &   & -  y_2 & \leq 12 \\
\end{array}
\f]
In this example, there three inequalities, so the dimension of the x vector would be 3 and n=3.  The input arrays would be as follows:
\code
int ik[]={0,3,5,7};
int row[]={0,1,2,0,1,0,1,2};
double vals[]={6.0,10.0,12.0,4.0,3.0,1.0,7.0,-1.0};
LPConeSetData(lpcone,3,ik,row,vals);
DSDPSetDualObjective(dsdp,1,1);
DSDPSetDualObjective(dsdp,2,1);
\endcode

*/
int LPConeSetData(LPCone lpcone,int n, const int ik[],const int cols[],const double vals[]){
  int info,i,spot,m=lpcone->m;
  DSDPVec C;
  DSDPFunctionBegin;
  lpcone->n=n;
  info=DSDPVecCreateSeq(n,&C);DSDPCHKERR(info);
  lpcone->C=C;
  info=DSDPVecZero(C);DSDPCHKERR(info);
  lpcone->muscale=1.0;
  if (n<100) lpcone->muscale=1.0;
  if (n<10) lpcone->muscale=1.0;
  for (i=ik[0];i<ik[1];i++){
    info=DSDPVecSetElement(C,cols[i],vals[i]);
  }
  spot=ik[0];
  info=CreateSpRowMatWdata(m,n,vals+spot,cols+spot,ik+1,&lpcone->A);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

int LPConeSetDataC(LPCone lpcone,int n, const double vals[]){
  int i;
  DSDPFunctionBegin;
  for (i=0;i<n;i++){
    lpcone->C.val[i]=vals[i];
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LPConeSetData2"
/*!
\fn LPConeSetData2(LPCone lpcone, int n, const int ik[],const int cols[],const double vals[])
\brief Set data A and into the LP cone.
\param lpcone the LP cone
\param n the number of inequalities
\param ik the number of nonzeros in each column of the matrix \f$A^T\f$
\param cols array of column numbers in A
\param vals array of nonzeros in A and c
\ingroup LPRoutines
\sa DSDPSetDualObjective()

For example, consider the following problem in the form of (D):
\f[ \begin{array}{llllll}
\mbox{Maximize} & & y_1 & + & y_2 \\
\mbox{Subject to}
 & & 4 y_1 & + & 2 y_2 & \leq 6 \\
 & & 3 y_1 & + & 7 y_2 & \leq 10 \\
 & &   &   & -  y_2 & \leq 12 \\
\end{array}
\f]

\code
int ik[]={0,2,5,7};
int row[]={0,1,0,1,2,0,1,2};
double vals[]={4.0,3.0,2.0,7.0,-1.0,6.0,10.0,12.0};
LPConeSetData2(lpcone,3,ik,row,vals);
DSDPSetDualObjective(dsdp,1,1);
DSDPSetDualObjective(dsdp,2,1);
\endcode
*/
int LPConeSetData2(LPCone lpcone,int n, const int ik[],const int cols[],const double vals[]){
  int info,i,spot,m=lpcone->m;
  DSDPVec C;
  DSDPFunctionBegin;
  lpcone->n=n;
  info=DSDPVecCreateSeq(n,&C);DSDPCHKERR(info);
  lpcone->C=C;
  info=DSDPVecZero(C);DSDPCHKERR(info);
  lpcone->muscale=1.0;
  if (n<100) lpcone->muscale=1.0;
  if (n<10) lpcone->muscale=1.0;
  for (i=ik[m];i<ik[m+1];i++){
    info=DSDPVecSetElement(C,cols[i],vals[i]);
  }
  spot=ik[0];
  info=CreateSpRowMatWdata(m,n,vals+spot,cols+spot,ik,&lpcone->A);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LPConeView2"
/*!
\fn LPConeView2(LPCone lpcone)
\brief Print the data in the LP cone to the screen.
\param lpcone the LP cone
\ingroup LPRoutines
 */
int LPConeView2(LPCone lpcone){
  int info;
  DSDPFunctionBegin;
  printf("LPCone Constraint Matrix\n");
  info=SpRowMatView(lpcone->A);DSDPCHKERR(info);
  printf("LPCone Objective C vector\n");
  info=DSDPVecView(lpcone->C);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}



#undef __FUNCT__
#define __FUNCT__ "LPConeGetConstraint"
int LPConeGetConstraint(LPCone lpcone,int column,DSDPVec W){
  int n,info;
  double *w;
  DSDPFunctionBegin;
  if (column==0){
    info=DSDPVecCopy(lpcone->C,W);DSDPCHKERR(info);
  } else {
    info=DSDPVecGetSize(W,&n);DSDPCHKERR(info);
    info=DSDPVecGetArray(W,&w);DSDPCHKERR(info);
    info=SpRowMatGetRowVector(lpcone->A,column-1,w,n);info=0;DSDPCHKERR(info);
    info=DSDPVecRestoreArray(W,&w);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LPConeGetData"
/*!
\fn int LPConeGetData(LPCone lpcone, int vari, double vv[], int n);
\brief Get one column (or row) of the LP data.
\param lpcone LP cone
\param vari column of data in (D) or row of data in (P).
\param vv array of data
\param n length of array and conic dimension
*/
int LPConeGetData(LPCone lpcone,int vari,double vv[], int n){
  int info;
  DSDPVec W;
  DSDPFunctionBegin;
  info=DSDPVecCreateWArray(&W,vv,n);DSDPCHKERR(info);
  info=LPConeGetConstraint(lpcone,vari,W);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "LPConeSetXVec"
int LPConeSetXVec(LPCone lpcone,double *xout, int n){
  DSDPFunctionBegin;
  if (n==lpcone->n) lpcone->xout=xout;
  DSDPFunctionReturn(0); 
}


static int vsdot(const int*,const double *,int,const double *, int, double *);
static int checkvsparse(smatx *);


static int CreateSpRowMatWdata(int m,int n,const double vals[], const int cols[], const int ik[],
			smatx **A){
  
  smatx *V;

  V=(smatx*) malloc(1*sizeof(smatx));
  if (V==NULL) return 1;
  V->nrow=m;
  V->ncol=n;
  V->owndata=0;
  V->an=vals; V->col=cols; V->nnz=ik;

  *A=V;
  checkvsparse(V);
  return 0;
}

static int vsdot(const int ja[], const double an[], int nn0, const double vv[], int n, double *vdot){

  int i;
  double tmp=0.0;

  for (i=0; i<nn0; i++){
    /*    if (ja[i]<n) tmp = tmp + an[i] * vv[ja[i]]; */
    tmp += an[i] * vv[ja[i]];
  }
  *vdot=tmp;
  return 0;
}

static int checkvsparse(smatx *A){
  int i,k=0,m=A->nrow,tnnz=0;
  const int *nnz=A->nnz;

  for (i=0;i<m;++i){
    if (nnz[i+1]-nnz[i]>0){
      tnnz++;
    }
  }
  if (tnnz < m/2){
    A->nzrows =(int*)malloc((tnnz)*sizeof(int));
    A->nnzrows=tnnz;
    for (i=0;i<m;++i){
      if (nnz[i+1]-nnz[i]>0){
	A->nzrows[k]=i;
	k++;
      }
    }    
  } else {
    A->nzrows = 0;
    A->nnzrows = m;
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SpRowMatMult"
static int SpRowMatMult(smatx* A, const double x[], int n, double y[], int m){

  int i,k1,k2,nrow=A->nrow;
  const int *nnz=A->nnz,*col=A->col;
  const double *an=A->an;

  if (A->ncol != n) return 1;
  if (A->nrow != m) return 2;
  if (x==0 && n>0) return 3;
  if (y==0 && m>0) return 3;

  if (m>0){
    memset((void*)y,0,m*sizeof(double));
    for (i=0; i<nrow; i++){
      k1=*(nnz+i); k2=*(nnz+i+1);
      vsdot(col+k1,an+k1,k2-k1,x,n,y+i);
    }
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SpRowMatMultTrans"
static int SpRowMatMultTrans(smatx * A, const double x[], int m, double y[], int n){

  int i,j,k1,k2,nrow=A->nrow;
  const int *col=A->col,*nnz=A->nnz;
  const double *an=A->an;
  if (A->ncol != n) return 1;
  if (A->nrow != m) return 2;
  if (x==0 && m>0) return 3;
  if (y==0 && n>0) return 3;

  memset((void*)y,0,n*sizeof(double));
  for (i=0; i<nrow; i++){
    k1=nnz[i]; k2=nnz[i+1];
    for (j=k1; j<k2; j++){
      y[col[j]] += an[j]*x[i]; 
    }
  }

  return 0;
}


static int SpRowMatRowNnz(smatx *A, int row, int* wi,int n){
  int i,k1,k2;
  const int *col=A->col;
  DSDPFunctionBegin;
  memset((void*)wi,0,n*sizeof(double));
  k1=A->nnz[row];
  k2=A->nnz[row+1];
  for (i=k1; i<k2; i++){
    wi[col[i]]=1;
  }
  DSDPFunctionReturn(0);
}

static int SpRowIMultAdd(smatx *A,int *wi,int n,int *rnnz,int m){

  int i,j,k1,k2,nrow=A->nrow;
  const int *nnz=A->nnz,*col=A->col;
  DSDPFunctionBegin;
  for (i=0; i<nrow; i++){
    k1=nnz[i];
    k2=nnz[i+1];
    for (j=k1; j<k2; j++){
      if (wi[col[j]]){
	rnnz[i]++;
      }
    }
  }
  DSDPFunctionReturn(0);
}
/*
static int SpRowMatAddRowMultiple(smatx* A, int nrow, double ytmp, double row[], int n){
  int k;
  int *col=A->col, *nnz=A->nnz;
  double *an=A->an;

  for (k=nnz[nrow]; k<nnz[nrow+1]; k++){
    row[col[k]] += ytmp * an[k];
  }
  
  return 0;
}
*/
static int SpRowMatNorm2(smatx* A, int nrow, double *norm22){
  int k;
  const int *nnz=A->nnz;
  double tt=0;
  const double *an=A->an;

  for (k=nnz[nrow]; k<nnz[nrow+1]; k++){
    tt+=an[k]*an[k];
  }
  *norm22=tt;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "SpRowMatGetRowVector"
static int SpRowMatGetRowVector(smatx* M, int row, double r[], int m){

  int i,k1,k2;
  const int *col=M->col;
  const double *an=M->an;
  /*  
  if (M->ncol != m) return 1;
  if (row<0 || row>= M->nrow) return 2;
  if (r==0) return 3;
  */
  memset((void*)r,0,m*sizeof(double));
  k1=M->nnz[row];
  k2=M->nnz[row+1];

  for (i=k1; i<k2; i++){
    r[col[i]]=an[i];
  }

  return 0;
}
#undef __FUNCT__
#define __FUNCT__ "SpRowMatGetScaledRowVector"
static int SpRowMatGetScaledRowVector(smatx* M, int row, const double ss[], double r[], int m){

  int i,k1,k2;
  const int *col=M->col;
  const double *an=M->an;
  /*  
  if (M->ncol != m) return 1;
  if (row<0 || row>= M->nrow) return 2;
  if (r==0) return 3;
  */
  memset((void*)r,0,m*sizeof(double));
  k1=M->nnz[row];
  k2=M->nnz[row+1];

  for (i=k1; i<k2; i++){
    r[col[i]]=ss[col[i]]*an[i];
  }

  return 0;
}


/*
#undef __FUNCT__
#define __FUNCT__ "SpRowMatZero"
static int SpRowMatZero(smatx* M){

  int nnz=M->nnz[M->nrow];
  memset(M->an,0,(nnz)*sizeof(double));

  return 0;
}



#undef __FUNCT__
#define __FUNCT__ "SpRowGetSize"
static int SpRowMatGetSize(smatx *M, int *m, int *n){
  *m=M->nrow;
  *n=M->ncol;

  return 0;
}
*/
#undef __FUNCT__
#define __FUNCT__ "SpRowMatDestroy"
static int SpRowMatDestroy(smatx* A){

  if (A->owndata){
    printf("Can't free array");
    return 1;
    /*
    if (A->an)  free(A->an);
    if (A->col) free(A->col);
    if (A->nnz) free(A->nnz);
    */
  }
  if (A->nzrows) free(A->nzrows);
  free(A);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SpRowMatView"
static int SpRowMatView(smatx* M){

  int i,j,k1,k2;

  for (i=0; i<M->nrow; i++){
    k1=M->nnz[i]; k2=M->nnz[i+1];

    if (k2-k1 >0){
      printf("Row %d, (Variable y%d) :  ",i,i+1);
      for (j=k1; j<k2; j++)
	printf(" %4.2e x%d + ",M->an[j],M->col[j]);
      printf("= dobj%d \n",i+1);
    }
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "LPConeView"
/*!
  \fn LPConeView(LPCone lpcone)
 \brief Print the data in the LP cone to the screen
\param lpcone the LP cone
\ingroup LPRoutines
 */
int LPConeView(LPCone lpcone){

  smatx* A=lpcone->A;
  int i,j,jj,info;
  const int *row=A->col,*nnz=A->nnz;
  int  n=A->ncol,m=A->nrow;
  const double *an=A->an;
  double cc;
  DSDPVec C=lpcone->C;
  DSDPFunctionBegin;
  printf("LPCone Constraint Matrix\n");
  printf("Number y variables 1 through %d\n",m);
  for (i=0; i<n; i++){
    printf("Inequality %d:  ",i);
    for (j=0;j<m;j++){
      for (jj=nnz[j];jj<nnz[j+1];jj++){
	if (row[jj]==i){
	  printf("%4.2e y%d + ",an[jj],j+1);
	}
      }
    }
    info=DSDPVecGetElement(C,i,&cc);DSDPCHKERR(info);
    printf(" <= %4.2e\n",cc);
  }
  DSDPFunctionReturn(0); 
}


