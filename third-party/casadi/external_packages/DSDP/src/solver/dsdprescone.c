#include "dsdpcone_impl.h"
#include "dsdpsys.h"
#include "dsdpbasictypes.h"
/*!
\file dsdprescone.c
\brief Variable r must be nonnegative.
*/

struct RDCone{
  double primalr,dualr,x,logr;
  DSDPPenalty UsePenalty;
  DSDP dsdp;  /* Really only need the Penalty flag, which is here. */
};

typedef struct RDCone RCone;

#undef __FUNCT__  
#define __FUNCT__ "DSDPRHessian"
static int DSDPRHessian( void *dspcone, double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2){
  RCone *K=(RCone*)dspcone;
  int info,m;
  double rr,sr,rssr;
  DSDPFunctionBegin;
  if (K->dualr ){
    info=DSDPVecGetSize(vrhs2,&m);DSDPCHKERR(info);
    info=DSDPSchurMatVariableCompute(M,m-1,&rr);DSDPCHKERR(info);
    if (rr){
      sr=-mu*rr/(K->dualr);
      rssr=mu*rr/(K->dualr*K->dualr);
      info=DSDPVecAddR(vrhs2,sr);DSDPCHKERR(info);
      info=DSDPSchurMatAddDiagonalElement(M,m-1,rssr);DSDPCHKERR(info);
    }
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPRHS"
static int DSDPRHS( void *dspcone, double mu, DSDPVec vrow,DSDPVec vrhs1,DSDPVec vrhs2){
  RCone *K=(RCone*)dspcone;
  int info;
  double rr,sr;
  DSDPFunctionBegin;
  if (K->dualr ){
    sr=-mu/(K->dualr);
    info=DSDPVecGetR(vrow,&rr);DSDPCHKERR(info);
    info=DSDPVecAddR(vrhs2,rr*sr);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPSetupRCone"
static int DSDPSetupRCone(void* dspcone,DSDPVec y){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPSetupRCone2"
static int DSDPSetupRCone2(void* dspcone, DSDPVec y, DSDPSchurMat M){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPDestroyRCone"
static int DSDPDestroyRCone(void* dspcone){
  int info;
  DSDPFunctionBegin;
  DSDPFREE(&dspcone,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPRSize"
static int DSDPRSize(void*dspcone,double *n){
  RCone *K=(RCone*)dspcone;
  DSDPFunctionBegin;
  if (K->dualr){*n=1;}
  else {*n=0;}
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPRSparsity"
static int DSDPRSparsity(void*dspcone,int row, int *tnnz, int rnnz[], int m){
  DSDPFunctionBegin;
  *tnnz=0;
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPComputeRS"
static int DSDPComputeRS(void *dspcone, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
  RCone *K=(RCone*)dspcone;
  int info;
  double rbeta;
  DSDPFunctionBegin;
  info=DSDPVecGetR(Y,&rbeta); DSDPCHKERR(info);
  if (K->UsePenalty==DSDPAlways){
    if (rbeta<0){ *ispsdefinite=DSDP_TRUE; } else { *ispsdefinite=DSDP_FALSE;} 
  } else {
    if (rbeta>0) rbeta=0;
    *ispsdefinite=DSDP_TRUE;
  }
  if (flag==DUAL_FACTOR){ K->dualr=rbeta; } else { K->primalr=rbeta;}
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
  RCone *K=(RCone*)dspcone;
  double r,rbeta,msteplength=1.0e100,rt=1.0e30;
  int info;
  DSDPFunctionBegin;

  info=DSDPVecGetR(DY,&rbeta); DSDPCHKERR(info);
  if (flag==DUAL_FACTOR){ r=K->dualr; } else { r=K->primalr;}
  if (r * rbeta<0) rt=-r/rbeta;

  if (K->UsePenalty==DSDPAlways){msteplength=rt;}
  else if (flag==PRIMAL_FACTOR){ msteplength=rt;}
  else if (flag==DUAL_FACTOR){msteplength=rt/0.94;}

  *maxsteplength=msteplength;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPRX"
static int DSDPRX( void *dspcone, double mu, DSDPVec y, DSDPVec dy, DSDPVec AX,double *tracexs){
  RCone *K=(RCone*)dspcone;
  int info;
  double rr,dr,trxs,r=K->dualr;
  DSDPFunctionBegin;
  
  info=DSDPVecGetR(y,&rr); DSDPCHKERR(info);
  info=DSDPVecGetR(dy,&dr); DSDPCHKERR(info);
  if (K->dualr){
    r=-1.0/rr;
    K->x=mu*(r-r*dr*r);
    trxs=K->x/r;
    DSDPLogInfo(0,2,"RESIDUAL X (Minimum Penalty Parameter): %4.4e, Trace(XS): %4.4e\n",K->x,trxs);
    /* *tracexs=trxs */
  } else {
    K->x=0.0;
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPSetX"
static int DSDPSetX( void *dspcone, double mu, DSDPVec y, DSDPVec dy){
  RCone *K=(RCone*)dspcone;
  int info;
  double rr,dr,trxs,r;
  DSDPFunctionBegin;
  
  info=DSDPVecGetR(y,&rr); DSDPCHKERR(info);
  info=DSDPVecGetR(dy,&dr); DSDPCHKERR(info);
  if (rr){
    r=-1.0/rr;
    K->x=mu*(r-r*dr*r);
    trxs=K->x/r;
    DSDPLogInfo(0,2,"RESIDUAL X (Minimum Penalty Parameter): %4.4e, Trace(XS): %4.4e\n",K->x,trxs);
    /* *tracexs=trxs */
  } else {
    K->x=0.0;
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPComputeRLog"
static int DSDPComputeRLog(void *dspcone, double *logobj, double *logdet){
  RCone *K=(RCone*)dspcone;
  DSDPFunctionBegin;
  *logdet=K->logr;
  *logobj=0;
  if (K->dualr<0){
    *logdet=log(-K->dualr);
    K->logr=log(-K->dualr);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPRANorm2"
static int DSDPRANorm2(void *dspcone,DSDPVec Anorm2){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPRMultiplyAdd"
static int DSDPRMultiplyAdd(void *dspcone,double mu,DSDPVec vrow,DSDPVec vin,DSDPVec vout){
  RCone *K=(RCone*)dspcone;
  int info;
  double v1,v2,rssr;
  DSDPFunctionBegin;
  if (K->dualr){
    info=DSDPVecGetR(vrow,&v1);DSDPCHKERR(info);
    info=DSDPVecGetR(vin,&v2);DSDPCHKERR(info);
    rssr=v1*v2*mu/(K->dualr*K->dualr);
    info=DSDPVecAddR(vout,rssr);DSDPCHKERR(info);  
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPRMonitor"
static int DSDPRMonitor( void *dspcone, int tag){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

static struct  DSDPCone_Ops kops;
static const char* matname="R Cone";

#undef __FUNCT__
#define __FUNCT__ "RConeOperationsInitialize"
static int RConeOperationsInitialize(struct  DSDPCone_Ops* coneops){
  int info;
  if (coneops==NULL) return 0;
  info=DSDPConeOpsInitialize(coneops); DSDPCHKERR(info);
  coneops->conehessian=DSDPRHessian;
  coneops->conesetup=DSDPSetupRCone;
  coneops->conesetup2=DSDPSetupRCone2;
  coneops->conedestroy=DSDPDestroyRCone;
  coneops->conecomputes=DSDPComputeRS;
  coneops->coneinverts=DSDPInvertRS;
  coneops->conesetxmaker=DSDPSetX;
  coneops->conecomputex=DSDPRX;
  coneops->conerhs=DSDPRHS;
  coneops->conemaxsteplength=DSDPComputeRStepLength;
  coneops->conelogpotential=DSDPComputeRLog;
  coneops->conesize=DSDPRSize;
  coneops->conesparsity=DSDPRSparsity;
  coneops->coneanorm2=DSDPRANorm2;
  coneops->conemonitor=DSDPRMonitor;
  coneops->conehmultiplyadd=DSDPRMultiplyAdd;
  coneops->id=19;
  coneops->name=matname;
  return 0;
}

/*
\fn int RConeSetType(RCone *rcone, DSDPPenalty UsePenalty);
\brief Set penalty type.

\param dsdp the solver
\param UsePenalty true or false
*/
#undef __FUNCT__
#define __FUNCT__ "RConeSetType"
int RConeSetType(RCone *rcone, DSDPPenalty UsePenalty){
  DSDPFunctionBegin;
  rcone->UsePenalty=UsePenalty;
  DSDPFunctionReturn(0); 
}
/*
\fn int RConeGetRX(RCone *rcone, double *rx);
\brief Get slack of trace of matrix

\param rcone cone
\param rx dual of r.
Accurate only when r > 0.

*/
#undef __FUNCT__
#define __FUNCT__ "RConeGetRX"
int RConeGetRX(RCone *rcone, double *xtrace){
  DSDPFunctionBegin;
  *xtrace=rcone->x;
  DSDPFunctionReturn(0); 
}

/*!
\fn int DSDPAddRCone(DSDP dsdp, RCone **rrcone);
\brief A separate cone specifies that r must be nonnegative.

\param dsdp the solver
\param rrcone new r cone
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPAddRCone"
int DSDPAddRCone(DSDP dsdp, RCone **rrcone){
  RCone *rcone;
  DSDPPenalty UsePenalty=DSDPInfeasible;
  int info;
  DSDPFunctionBegin;
  info=RConeOperationsInitialize(&kops); DSDPCHKERR(info);
  DSDPCALLOC1(&rcone,RCone,&info); DSDPCHKERR(info);
  info=RConeSetType(rcone,UsePenalty); DSDPCHKERR(info);
  rcone->dsdp=dsdp;
  rcone->logr=0;
  *rrcone=rcone;
  info=DSDPAddCone(dsdp,&kops,(void*)rcone); DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

