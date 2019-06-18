#include "dsdpcone_impl.h"
#include "dsdpsys.h"
/*!
\file dsdpobjcone.c
\brief Apply a constraint that the objective solution (D) is greater than
some lower bound.
 */
typedef struct {
  DSDPVec     b,bb,T;
  double      dmin;
  double      pss,dss;
  DSDP dsdp;
  DSDPTruth useit;
} BDCone;


static int BComputeS(BDCone *K, DSDPVec v, double *ss){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecDot(K->bb,v,ss);DSDPCHKERR(info); 
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPRHessian"
static int DSDPRHessian( void *dspcone, double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2){
  BDCone *K=(BDCone*)dspcone;
  double bb,dd,ss=K->dss;
  int info,i,m,ncols;
  DSDPVec T=K->T,b=K->bb;
  DSDPFunctionBegin;
  if (K->useit){
    info=DSDPVecGetSize(b,&m);DSDPCHKERR(info);
    for (i=0;i<m;i++){
      info=DSDPVecGetElement(b,i,&bb);DSDPCHKERR(info);
      if (bb==0) continue;

      info=DSDPSchurMatRowColumnScaling(M,i,T,&ncols); DSDPCHKERR(info); 
      if (ncols==0) continue;

      info=DSDPVecGetElement(T,i,&dd);DSDPCHKERR(info);
      info=DSDPVecAddElement(vrhs2,i,-bb*dd*mu/ss);DSDPCHKERR(info);

      info=DSDPVecPointwiseMult(T,b,T);DSDPCHKERR(info);

      info=DSDPVecScale(mu*bb/(ss*ss),T);DSDPCHKERR(info);
      info=DSDPSchurMatAddRow(M,i,1.0,T);DSDPCHKERR(info);
      /*
      info=DSDPSchurMatAddRow(M,i,mu*bb/(ss*ss),T);DSDPCHKERR(info);
      */
      if (i==-m-1){info=DSDPVecView(T);}
    }
  }
  DSDPFunctionReturn(0);
}

static int DSDPRRHS( void *dcone, double mu, DSDPVec vrow, DSDPVec vrhs1, DSDPVec vrhs2){
  BDCone *K=(BDCone*)dcone;
  double bb,dd,ss=K->dss;
  int info,i,m;
  DSDPVec b=K->bb;
  DSDPFunctionBegin;
  if (K->useit){
    info=DSDPVecGetSize(b,&m);DSDPCHKERR(info);
    for (i=0;i<m;i++){
      info=DSDPVecGetElement(b,i,&bb);DSDPCHKERR(info);
      info=DSDPVecGetElement(vrow,i,&dd);DSDPCHKERR(info);
      info=DSDPVecAddElement(vrhs2,i,-bb*dd*mu/ss);DSDPCHKERR(info);
    }
    /*
    info=DSDPVecGetR(b,&bb);DSDPCHKERR(info);
    info=DSDPVecGetR(vrhs3,&dd);DSDPCHKERR(info);
    info=DSDPVecPointwiseMult(vrow,b,T);DSDPCHKERR(info);
    info=DSDPVecScale(mu*bb/(ss*ss),T);DSDPCHKERR(info);
    info=DSDPVecAXPY(dd,T,vrhs3);DSDPCHKERR(info);
    */
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPSetupBCone"
static int DSDPSetupBCone(void* dspcone,DSDPVec y){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDRData"
static int DSDPSetDRData(BDCone *K){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecCopy(K->b,K->bb);DSDPCHKERR(info);
  info=DSDPVecSetC(K->bb,K->dmin);DSDPCHKERR(info);
  info=DSDPVecSetR(K->bb,-1.0);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPSetupBCone2"
static int DSDPSetupBCone2(void* dspcone, DSDPVec y, DSDPSchurMat M){
  BDCone *K=(BDCone*)dspcone;
  int info;
  DSDPFunctionBegin;
  info=DSDPVecDuplicate(K->b,&K->T);DSDPCHKERR(info);
  info=DSDPVecDuplicate(K->b,&K->bb);DSDPCHKERR(info);
  info=DSDPSetDRData(K);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPDestroyBCone"
static int DSDPDestroyBCone(void* dspcone){
  BDCone *K=(BDCone*)dspcone;
  int info;
  DSDPFunctionBegin;
  info=DSDPVecDestroy(&K->T);DSDPCHKERR(info);
  info=DSDPVecDestroy(&K->bb);DSDPCHKERR(info);
  DSDPFREE(&dspcone,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPRSize"
static int DSDPRSize(void*dspcone,double*n){
  DSDPFunctionBegin;
  *n=1.0;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPRSparsity"
static int DSDPRSparsity(void*dspcone,int row, int *tnnz, int rnnz[], int m){
  BDCone *K=(BDCone*)dspcone;
  int i,info;double dd;
  DSDPFunctionBegin;
  *tnnz=0;

  info=DSDPVecGetElement(K->b,row,&dd);DSDPCHKERR(info);
  if (dd){
    for (i=0;i<m;i++){
      info=DSDPVecGetElement(K->b,i,&dd);DSDPCHKERR(info);
      if (dd!=0){rnnz[i]++; (*tnnz)++;}
    }
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPComputeRS"
static int DSDPComputeRS(void *dspcone, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
  BDCone *K=(BDCone*)dspcone;
  int info;
  double ss;
  DSDPFunctionBegin;
  info=BComputeS(K,Y,&ss);DSDPCHKERR(info);
  if (ss>0){ *ispsdefinite=DSDP_TRUE; } else { *ispsdefinite=DSDP_FALSE;} 
  if (flag==DUAL_FACTOR){ K->dss=ss; } else { K->pss=ss;}
  DSDPLogInfo(0,2,"DOBJCone SS: %4.4e \n",ss);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPInvertRS"
static int DSDPInvertRS(void *dspcone){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPComputeRStepLength"
static int DSDPComputeRStepLength(void *dspcone, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
  BDCone *K=(BDCone*)dspcone;
  double ds,ss,rt=1.0e30;
  int info;

  DSDPFunctionBegin;
  info=BComputeS(K,DY,&ds);DSDPCHKERR(info);
  if (flag==DUAL_FACTOR){ ss=K->dss; } else { ss=K->pss;}
  if (ds<0) rt=-ss/ds;
  if (K->useit){
    *maxsteplength=rt;
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSetX"
static int DSDPSetX( void *dspcone, double mu, DSDPVec y, DSDPVec dy){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DSDPRX"
static int DSDPRX( void *dspcone, double mu, DSDPVec y, DSDPVec dy, DSDPVec AX,double*tracxs){
  BDCone *K=(BDCone*)dspcone;
  double x,dss,ss=K->dss;
  int info;
  DSDPFunctionBegin;
  
  info=BComputeS(K,y,&ss);DSDPCHKERR(info);
  ss=1.0/ss;
  info=BComputeS(K,dy,&dss);DSDPCHKERR(info);
  x=mu*(ss+ss*dss*ss);
  DSDPLogInfo(0,2,"DOBJCone SS: %4.4e, RESIDUAL X: %4.4e\n",1.0/ss,x);
  if (fabs(x*ss)>1.0 && mu < 1) printf("Check Dual Min Bound\n");
  info=DSDPVecAXPY(-x,K->bb,AX);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPComputeRLog"
static int DSDPComputeRLog(void *dspcone, double *logobj, double *logdet){
  BDCone *K=(BDCone*)dspcone;
  DSDPFunctionBegin;
  *logobj=0;
  *logdet=log(K->dss);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPRANorm2"
static int DSDPRANorm2(void *dspcone, DSDPVec Anorm2){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPRMultiplyAdd"
static int DSDPRMultiplyAdd(void *dspcone,double mu,DSDPVec vrow,DSDPVec vin,DSDPVec vout){
  BDCone *K=(BDCone*)dspcone;
  DSDPVec T=K->T;
  int info;
  double dd,ss=K->dss;
  DSDPFunctionBegin;
  info=DSDPVecDot(vin,K->bb,&dd);DSDPCHKERR(info);
  dd=-mu*dd/(ss*ss);
  info=DSDPVecPointwiseMult(K->bb,vrow,T);DSDPCHKERR(info);
  info=DSDPVecAXPY(dd,T,vout);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPRMonitor"
static int DSDPRMonitor( void *dspcone, int tag){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

static struct  DSDPCone_Ops kops;
static const char* matname="Dual Obj Cone";

#undef __FUNCT__
#define __FUNCT__ "BConeOperationsInitialize"
static int BConeOperationsInitialize(struct  DSDPCone_Ops* coneops){
  int info;
  if (coneops==NULL) return 0;
  info=DSDPConeOpsInitialize(coneops); DSDPCHKERR(info);
  coneops->conehessian=DSDPRHessian;
  coneops->conesetup=DSDPSetupBCone;
  coneops->conesetup2=DSDPSetupBCone2;
  coneops->conedestroy=DSDPDestroyBCone;
  coneops->conecomputes=DSDPComputeRS;
  coneops->coneinverts=DSDPInvertRS;
  coneops->conecomputex=DSDPRX;
  coneops->conesetxmaker=DSDPSetX;
  coneops->conemaxsteplength=DSDPComputeRStepLength;
  coneops->conelogpotential=DSDPComputeRLog;
  coneops->conesize=DSDPRSize;
  coneops->conesparsity=DSDPRSparsity;
  coneops->coneanorm2=DSDPRANorm2;
  coneops->conemonitor=DSDPRMonitor;
  coneops->conehmultiplyadd=DSDPRMultiplyAdd;
  coneops->conerhs=DSDPRRHS;
  coneops->id=119;
  coneops->name=matname;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPAddBCone"
int DSDPAddBCone(DSDP dsdp, DSDPVec bb, double dmin){
  BDCone *rcone;
  int info;
  DSDPFunctionBegin;
  info=BConeOperationsInitialize(&kops); DSDPCHKERR(info);
  DSDPCALLOC1(&rcone,BDCone,&info); DSDPCHKERR(info);
  rcone->b=bb;
  rcone->dmin=dmin;
  rcone->dsdp=dsdp;
  rcone->useit=DSDP_TRUE;
  info=DSDPAddCone(dsdp,&kops,(void*)rcone); DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#include "dsdp.h"
#include "dsdp5.h"

#undef __FUNCT__
#define __FUNCT__ "DSDPSetDualLowerBound"
int DSDPSetDualLowerBound(DSDP dsdp, double dobjmin){
  int info;
  DSDPFunctionBegin;
  info = DSDPAddBCone(dsdp,dsdp->b,dobjmin);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}
