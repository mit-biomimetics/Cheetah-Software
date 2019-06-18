#include "dsdpsdp.h"
#include "dsdpcone_impl.h"
#include "dsdpsys.h"
/*! \file sdpkcone.c 
\brief Implement the DSDPCone operations using the SDPCone subroutines.
 */
static int SDPConeOperationsInitialize(struct  DSDPCone_Ops*);

static int KSDPConeSetup(void*,DSDPVec);
static int KSDPConeSetup2(void*, DSDPVec, DSDPSchurMat);
static int KSDPConeSize(void*,double*);
static int KSDPConeSparsity(void*,int,int[],int[],int);
static int KSDPConeComputeHessian(void*,double,DSDPSchurMat,DSDPVec,DSDPVec);
static int KSDPConeComputeMaxStepLength(void*, DSDPVec, DSDPDualFactorMatrix, double *);
static int KSDPConeComputeSS(void*, DSDPVec, DSDPDualFactorMatrix, DSDPTruth *);
static int KSDPConeComputeLogSDeterminant(void *, double *, double*);
static int KSDPConeComputeXX(void*, double, DSDPVec,DSDPVec,DSDPVec,double*);
static int KSDPConeDestroy(void*);

#undef __FUNCT__  
#define __FUNCT__ "KSDPConeComputeHessian"
static int KSDPConeComputeHessian( void *K, double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2){
  int info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  info=SDPConeComputeHessian(sdpcone,mu,M,vrhs1,vrhs2);DSDPCHKERR(info);
  DSDPFunctionReturn(0);   
}

#undef __FUNCT__  
#define __FUNCT__ "KDPConeMultiply"
static int KSDPConeMultiply( void *K, double mu, DSDPVec vrow, DSDPVec vin, DSDPVec vout){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  SDPConeValid(sdpcone);
  DSDPFunctionBegin;
  for (kk=0; kk<sdpcone->nblocks; kk++){      
    info=SDPConeMultiply(sdpcone,kk,mu,vrow,vin,vout);DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0);   
}

#undef __FUNCT__  
#define __FUNCT__ "KDPConeRHS  "
static int KSDPConeRHS( void *K, double mu, DSDPVec vrow, DSDPVec vrhs1, DSDPVec vrhs2){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0; kk<sdpcone->nblocks; kk++){      
    if (sdpcone->blk[kk].n<1) continue;
    info=SDPConeComputeRHS(sdpcone,kk,mu,vrow,vrhs1,vrhs2); DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0);   
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSetup"
static int KSDPConeSetup(void* K, DSDPVec y){
  int info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=SDPConeSetup(sdpcone,y);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSetup2"
static int KSDPConeSetup2(void* K, DSDPVec y, DSDPSchurMat M){
  int info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  info=SDPConeSetup2(sdpcone,y,M); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeDestroy"
static int KSDPConeDestroy(void* K){
  int info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=SDPConeDestroy(sdpcone);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSize"
static int KSDPConeSize(void* K,double *n){
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  *n=sdpcone->nn;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSparsity"
static int KSDPConeSparsity(void *K,int row, int *tnnz, int rnnz[], int m){
  SDPCone sdpcone=(SDPCone)K;
  SDPblk *blk=sdpcone->blk;
  int info,j,kk;
  int nnzblocks=sdpcone->ATR.nnzblocks[row],*nzblocks=sdpcone->ATR.nzblocks[row];
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (j=0; j<nnzblocks; j++){
    kk=nzblocks[j];
    if (blk[kk].n<1) continue;
    info=DSDPBlockDataMarkNonzeroMatrices(&blk[kk].ADATA,rnnz);DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "KSDPConeComputeSS"
static int KSDPConeComputeSS(void *K, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  SDPblk *blk=sdpcone->blk;
  DSDPTruth psdefinite;
  DSDPDualMat SS;

  DSDPFunctionBegin;
  psdefinite = DSDP_TRUE;
  for (kk=sdpcone->nblocks-1; kk>=0 && psdefinite == DSDP_TRUE; kk--){
    if (blk[kk].n<1) continue;
    if (flag==DUAL_FACTOR) SS=blk[kk].S;
    else SS=blk[kk].SS;

    switch (sdpcone->optype){
    default:
      info=SDPConeComputeSS(sdpcone,kk,Y,blk[kk].T);DSDPCHKBLOCKERR(kk,info);
      info=DSDPDualMatSetArray(SS,blk[kk].T); DSDPCHKBLOCKERR(kk,info);
      info=DSDPDualMatCholeskyFactor(SS,&psdefinite); DSDPCHKBLOCKERR(kk,info);
      if (psdefinite == DSDP_FALSE){
	if (flag==DUAL_FACTOR){
	  DSDPLogInfo(0,2,"Dual SDP Block %2.0f not PSD\n",kk);
	} else {
	  DSDPLogInfo(0,2,"Primal SDP Block %2.0f not PSD\n",kk);
	}
      }
      break;
    }
  }
  *ispsdefinite=psdefinite;
  if (flag==DUAL_FACTOR && psdefinite==DSDP_TRUE){
    info=DSDPVecCopy(Y,sdpcone->YY);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSDPConeInvertSS"
static int KSDPConeInvertSS(void *K){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPDualMat SS;
  
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0;kk<sdpcone->nblocks;kk++){
    if (sdpcone->blk[kk].n<1) continue;
    SS=sdpcone->blk[kk].S;
    info=DSDPDualMatInvert(SS);DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "KSDPConeComputeMaxStepLength"
static int KSDPConeComputeMaxStepLength(void *K, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
  int kk,info;
  double smaxstep,maxmaxstep=1.0e20;
  SDPCone sdpcone=(SDPCone)K;
  DSDPDualMat SS;
  SDPblk *blk=sdpcone->blk;
  DSDPDSMat DS;
  DSDPVMat T;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0; kk<sdpcone->nblocks; kk++){      
    if (blk[kk].n<1) continue;
    if (flag==DUAL_FACTOR) SS=blk[kk].S;
    else SS=blk[kk].SS;
    DS=blk[kk].DS; T=blk[kk].T;

    info=SDPConeComputeSS(sdpcone,kk,DY,T);DSDPCHKBLOCKERR(kk,info);
    info=DSDPDSMatSetArray(DS,T); DSDPCHKBLOCKERR(kk,info);

    info=DSDPLanczosStepSize( &blk[kk].Lanczos,blk[kk].W,blk[kk].W2,SS,DS,&smaxstep );DSDPCHKBLOCKERR(kk,info);
    DSDPLogInfo(0,12,"Block %d, PD %d, maxstepsize: %4.4e\n",kk,flag,smaxstep);
    maxmaxstep=DSDPMin(smaxstep,maxmaxstep); 
  }
  *maxsteplength=maxmaxstep;
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "KSDPConeAddANorm2"
static int KSDPConeAddANorm2(void *K, DSDPVec ANorm2){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  SDPblk *blk=sdpcone->blk;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0; kk<sdpcone->nblocks; kk++){      
    if (blk[kk].n<1) continue;
    info=DSDPBlockANorm2( &blk[kk].ADATA,ANorm2,blk[kk].n); DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSetX"
static int KSDPConeSetX(void *K, double mu, DSDPVec Y,DSDPVec DY){
  SDPCone sdpcone=(SDPCone)K;
  int info;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=DSDPVecCopy(Y,sdpcone->YX);DSDPCHKERR(info);
  info=DSDPVecCopy(DY,sdpcone->DYX);DSDPCHKERR(info);
  sdpcone->xmakermu=mu;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeComputeXX"
static int KSDPConeComputeXX(void *K, double mu, DSDPVec Y,DSDPVec DY,DSDPVec AX,double* tracexs){

  SDPCone sdpcone=(SDPCone)K;
  int info,kk;
  double xnorm,trxs,xtrace;
  DSDPVMat X;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=KSDPConeSetX(K,mu,Y,DY);DSDPCHKERR(info);
  for (kk=0; kk<sdpcone->nblocks; kk++){
    if (sdpcone->blk[kk].n<1) continue;
    X=sdpcone->blk[kk].T;
    info=SDPConeComputeX3(sdpcone,kk,mu,Y,DY,X);DSDPCHKBLOCKERR(kk,info);
    info=SDPConeComputeXDot(sdpcone,kk,Y,X,AX,&xtrace,&xnorm,&trxs);DSDPCHKBLOCKERR(kk,info);
    *tracexs+=trxs;
    DSDPLogInfo(0,10,"SDP Block %d: norm(X): %4.2e, trace(X): %4.2e, trace(XS): %4.2e.\n",kk,xnorm,xtrace,trxs);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "KSDPConeComputeLogSDeterminant"
static int KSDPConeComputeLogSDeterminant(void *K, double *logdetobj, double *logdet){
  int kk,info;
  double dlogdet=0,dlogdet2=0,dd;
  SDPCone sdpcone=(SDPCone)K;
  SDPblk *blk=sdpcone->blk;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0; kk<sdpcone->nblocks; kk++){
    if (blk[kk].n<1) continue;
    info=DSDPDualMatLogDeterminant(blk[kk].S,&dd);DSDPCHKBLOCKERR(kk,info);
    dlogdet+=dd*blk[kk].gammamu;
    dlogdet2+=dd*blk[kk].bmu;
  }
  *logdet=dlogdet;
  *logdetobj=dlogdet2;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeMonitor"
int KSDPConeMonitor(void *K, int tag){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

static struct DSDPCone_Ops kops;
static const char *sdpconename ="SDP Cone";

#undef __FUNCT__
#define __FUNCT__ "SDPConeOperationsInitialize"
static int SDPConeOperationsInitialize(struct  DSDPCone_Ops* coneops){
  int info;
  if (coneops==NULL) return 0;
  info=DSDPConeOpsInitialize(coneops); DSDPCHKERR(info);
  coneops->conehessian=KSDPConeComputeHessian;
  coneops->conerhs=KSDPConeRHS;
  coneops->conesetup=KSDPConeSetup;
  coneops->conesetup2=KSDPConeSetup2;
  coneops->conedestroy=KSDPConeDestroy;
  coneops->conecomputes=KSDPConeComputeSS;
  coneops->coneinverts=KSDPConeInvertSS;
  coneops->conesetxmaker=KSDPConeSetX;
  coneops->conecomputex=KSDPConeComputeXX;
  coneops->conemaxsteplength=KSDPConeComputeMaxStepLength;
  coneops->conelogpotential=KSDPConeComputeLogSDeterminant;
  coneops->conesize=KSDPConeSize;
  coneops->conesparsity=KSDPConeSparsity;
  coneops->conehmultiplyadd=KSDPConeMultiply;
  coneops->coneanorm2=KSDPConeAddANorm2;
  coneops->conemonitor=KSDPConeMonitor;
  coneops->id=1;
  coneops->name=sdpconename;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPAddSDP"
/*!
\fn int DSDPAddSDP(DSDP dsdp, SDPCone sdpcone);
\brief Pass a semidefinite cone to the solver.
\param dsdp solver
\param sdpcone semidefinite cone
*/
int DSDPAddSDP(DSDP dsdp,SDPCone sdpcone){
  int info;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=SDPConeOperationsInitialize(&kops); DSDPCHKERR(info);
  info=DSDPAddCone(dsdp,&kops,(void*)sdpcone); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

