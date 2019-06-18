#include "dsdpcg.h"
#include "dsdpschurmat.h"
#include "dsdpvec.h"
#include "dsdpsys.h"
#include "dsdp.h"
/*!
  \file dsdpcg.c
  \brief Apply Conjugate Gradient method to solve the Schur matrix.
*/


typedef enum { DSDPNoMatrix=1, DSDPUnfactoredMatrix=2, DSDPFactoredMatrix=3} DSDPCGType;

typedef struct{
  DSDPCGType type;
  DSDPSchurMat M;
  DSDPVec Diag;
  DSDP dsdp;
} DSDPCGMat;

#undef __FUNCT__  
#define __FUNCT__ "DSDPCGMatMult"
int DSDPCGMatMult(DSDPCGMat M, DSDPVec X, DSDPVec Y){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecZero(Y); DSDPCHKERR(info);
  if (M.type==DSDPUnfactoredMatrix){
    info=DSDPSchurMatMultiply(M.M,X,Y); DSDPCHKERR(info);
  } else if (M.type==DSDPFactoredMatrix){
    info=DSDPSchurMatMultR(M.M,X,Y); DSDPCHKERR(info);
    info=DSDPVecAXPY(-0*M.dsdp->Mshift,X,Y); DSDPCHKERR(info);
  } else if (M.type==DSDPNoMatrix){
    info=DSDPHessianMultiplyAdd(M.dsdp,X,Y);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPCGMatPre"
int DSDPCGMatPre(DSDPCGMat M, DSDPVec X, DSDPVec Y){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecZero(Y); DSDPCHKERR(info);
  if (M.type==DSDPUnfactoredMatrix){
    info=DSDPVecPointwiseMult(X,M.Diag,Y);DSDPCHKERR(info);
    info=DSDPVecPointwiseMult(Y,M.Diag,Y);DSDPCHKERR(info);
  } else if (M.type==DSDPFactoredMatrix){
    info=DSDPSchurMatSolve(M.M,X,Y);DSDPCHKERR(info);
  } else if (M.type==DSDPNoMatrix){
    info=DSDPVecCopy(X,Y);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPCGMatPreLeft"
int DSDPCGMatPreLeft(DSDPCGMat M, DSDPVec X, DSDPVec Y){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecZero(Y); DSDPCHKERR(info);
  if (M.type==DSDPUnfactoredMatrix){
    info=DSDPVecPointwiseMult(X,M.Diag,Y);DSDPCHKERR(info);
  } else if (M.type==DSDPFactoredMatrix){
    info=DSDPSchurMatSolve(M.M,X,Y);DSDPCHKERR(info);
  } else if (M.type==DSDPNoMatrix){
    info=DSDPVecCopy(X,Y);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPCGMatPreRight"
int DSDPCGMatPreRight(DSDPCGMat M, DSDPVec X, DSDPVec Y){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecZero(Y); DSDPCHKERR(info);
  if (M.type==DSDPNoMatrix){
    info=DSDPVecPointwiseMult(X,M.Diag,Y);DSDPCHKERR(info);
  } else if (M.type==DSDPFactoredMatrix){
    info=DSDPVecCopy(X,Y);DSDPCHKERR(info);
  } else if (M.type==DSDPUnfactoredMatrix){
    info=DSDPVecCopy(X,Y);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPConjugateResidual"
int DSDPConjugateResidual(DSDPCGMat B, DSDPVec X, DSDPVec D, DSDPVec R, DSDPVec BR, DSDPVec P, DSDPVec BP, DSDPVec TT3, int maxits, int *iter){

  int i,n,info;
  double zero=0.0,minus_one=-1.0;
  double alpha,beta,bpbp,rBr,rBrOld;
  double r0,rerr=1.0e20;

  DSDPFunctionBegin;
  info=DSDPVecNorm2(X,&rBr); DSDPCHKERR(info);
  if (rBr>0){
    info=DSDPVecCopy(X,P); DSDPCHKERR(info);
    info=DSDPCGMatPreRight(B,P,X);DSDPCHKERR(info);
    info=DSDPCGMatMult(B,X,R); DSDPCHKERR(info);
  } else {
    info=DSDPVecSet(zero,R); DSDPCHKERR(info);
  }
  info=DSDPVecAYPX(minus_one,D,R); DSDPCHKERR(info);

  info=DSDPCGMatPreLeft(B,D,R);DSDPCHKERR(info);
  info=DSDPVecCopy(R,P); DSDPCHKERR(info);

  info=DSDPCGMatPreRight(B,R,BR);DSDPCHKERR(info);
  info=DSDPCGMatMult(B,BR,TT3); DSDPCHKERR(info);
  info=DSDPCGMatPreRight(B,TT3,BR);DSDPCHKERR(info);

  info=DSDPVecCopy(BR,BP); DSDPCHKERR(info);
  info=DSDPVecDot(BR,R,&rBr); DSDPCHKERR(info);
  info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
  r0=rBr;

  for (i=0;i<maxits;i++){

    if (rerr/n < 1.0e-30 || rBr/n <= 1.0e-30  || rBr< 1.0e-12 * r0 ) break;

    info=DSDPVecDot(BP,BP,&bpbp); DSDPCHKERR(info);
    alpha = rBr / bpbp;
    info= DSDPVecAXPY(alpha,P,X); DSDPCHKERR(info);
    alpha = -alpha;
    info=DSDPVecAXPY(alpha,BP,R); DSDPCHKERR(info);

    info=DSDPCGMatPreRight(B,R,BR);DSDPCHKERR(info);
    info=DSDPCGMatMult(B,BR,TT3); DSDPCHKERR(info);
    info=DSDPCGMatPreLeft(B,TT3,BR);DSDPCHKERR(info);
    
    rBrOld=rBr;
    info=DSDPVecNorm2(R,&rerr); DSDPCHKERR(info);
    info=DSDPVecDot(BR,R,&rBr); DSDPCHKERR(info);

    DSDPLogInfo(0,11,"CG: rerr: %4.4e, rBr: %4.4e \n",rerr,rBr);  

    beta = rBr/rBrOld;
    info=DSDPVecAYPX(beta,R,P); DSDPCHKERR(info);
    info=DSDPVecAYPX(beta,BR,BP); DSDPCHKERR(info);
    
  }
  info=DSDPVecCopy(X,BR);DSDPCHKERR(info);
  info=DSDPCGMatPreRight(B,BR,X);DSDPCHKERR(info);

  DSDPLogInfo(0,2,"Conjugate Residual, Initial rMr, %4.2e, Final rMr: %4.2e, Iterates: %d\n",r0,rBr,i);

  *iter=i;

  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPConjugateGradient"
int DSDPConjugateGradient(DSDPCGMat B, DSDPVec X, DSDPVec D, DSDPVec R, DSDPVec BR, DSDPVec P, DSDPVec BP, DSDPVec Z, double cgtol, int maxits, int *iter){

  int i,n,info;
  double alpha,beta=0,bpbp;
  double r0,rerr=1.0e20;
  double rz,rzold,rz0;

  DSDPFunctionBegin;
  if (maxits<=0){
    *iter=0;
    DSDPFunctionReturn(0);
  } 
  info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
  info=DSDPVecNorm2(X,&rerr); DSDPCHKERR(info);
  info=DSDPVecZero(R); DSDPCHKERR(info);
  if (rerr>0){
    info=DSDPCGMatMult(B,X,R);DSDPCHKERR(info);
  }
  alpha=-1.0;
  info=DSDPVecAYPX(alpha,D,R); DSDPCHKERR(info);
  info=DSDPVecNorm2(R,&rerr); DSDPCHKERR(info);    
  if (rerr*sqrt(1.0*n)<1e-11 +0*cgtol*cgtol){
    *iter=1;
    DSDPFunctionReturn(0);
  } 

  if (rerr>0){
    info=DSDPVecCopy(R,P); DSDPCHKERR(info);
    info=DSDPVecSetR(P,0);DSDPCHKERR(info);
    info=DSDPVecNorm2(P,&rerr); DSDPCHKERR(info);    
    info=DSDPCGMatPre(B,R,Z);DSDPCHKERR(info);
  }
    
  info=DSDPVecCopy(Z,P); DSDPCHKERR(info);
  info=DSDPVecDot(R,Z,&rz); DSDPCHKERR(info);
  r0=rerr;rz0=rz;
  
  for (i=0;i<maxits;i++){
    
    info=DSDPCGMatMult(B,P,BP);DSDPCHKERR(info);
    info=DSDPVecDot(BP,P,&bpbp); DSDPCHKERR(info);
    if (bpbp==0) break;
    alpha = rz/bpbp;
    info=DSDPVecAXPY(alpha,P,X); DSDPCHKERR(info);
    info=DSDPVecAXPY(-alpha,BP,R); DSDPCHKERR(info);
    info=DSDPVecNorm2(R,&rerr); DSDPCHKERR(info);

    DSDPLogInfo(0,15,"CG: iter: %d, rerr: %4.4e, alpha: %4.4e, beta=%4.4e, rz: %4.4e \n",i,rerr,alpha,beta,rz);  
    if (i>1){
      if (rerr*sqrt(1.0*n) < cgtol){ break;}
      if (rz*n < cgtol*cgtol){ break;}
      if (rerr/r0 < cgtol){ break;}
    }
    if (rerr<=0){ break;}
    info=DSDPCGMatPre(B,R,Z);DSDPCHKERR(info);
    rzold=rz;
    info=DSDPVecDot(R,Z,&rz); DSDPCHKERR(info);
    beta=rz/rzold;
    info= DSDPVecAYPX(beta,Z,P); DSDPCHKERR(info);
  }
  DSDPLogInfo(0,2,"Conjugate Gradient, Initial ||r||: %4.2e, Final ||r||: %4.2e, Initial ||rz||: %4.4e, ||rz||: %4.4e, Iterates: %d\n",r0,rerr,rz0,rz,i+1);

  *iter=i+1;

  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success);

\brief Apply CG to solve for the step directions.

\param dsdp the solver
\param MM matrix
\param RHS right-hand side
\param X solution
\param cgtol accuracy
\param success output whether a solution of suitable accuracy was found

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPCGSolve"
int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
  int iter=0,n,info,maxit=10;
  double dd,ymax;
  DSDPCG *sles=dsdp->sles; 
  DSDPCGMat    CGM;
  DSDPFunctionBegin;

  info=DSDPEventLogBegin(dsdp->cgtime);
  info=DSDPVecZero(X);DSDPCHKERR(info);
  info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
  *success=DSDP_TRUE;
  if (0){
    maxit=0;
  } else if (dsdp->slestype==1){

    CGM.type=DSDPNoMatrix;
    CGM.M=MM;
    CGM.dsdp=dsdp;
    cgtol*=1000;
    maxit=5;

  } else if (dsdp->slestype==2){
    CGM.type=DSDPUnfactoredMatrix;
    CGM.M=MM;
    CGM.Diag=sles->Diag;
    CGM.dsdp=dsdp;
    cgtol*=100;
    maxit=(int)sqrt(1.0*n)+10;
    if (maxit>20) maxit=20;
    info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
    /*
      info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
      info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
    */
    
  } else if (dsdp->slestype==3){
    CGM.type=DSDPFactoredMatrix;
    CGM.M=MM;
    CGM.dsdp=dsdp;
    maxit=0;
    info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
    if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
    if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
      maxit=6;
    } else if (dsdp->rgap<1e-5){
      maxit=3;
    }

    info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);

  } else if (dsdp->slestype==4){
    CGM.type=DSDPFactoredMatrix;
    CGM.M=MM;
    CGM.dsdp=dsdp;
    maxit=3;
    info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
  }
  if (n<maxit) maxit=n;
  
  info=DSDPConjugateGradient(CGM,X,RHS,
			     sles->R,sles->BR,sles->P,sles->BP,
			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);

  if (iter>=maxit) *success=DSDP_FALSE;
  info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
  if (dd<0) *success=DSDP_FALSE;
  info=DSDPEventLogEnd(dsdp->cgtime);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPCGInitialize"
int DSDPCGInitialize(DSDPCG **sles){
  DSDPCG *ssles;
  int info;
  DSDPFunctionBegin;
  DSDPCALLOC1(&ssles,DSDPCG,&info);DSDPCHKERR(info);
  ssles->setup2=0;
  *sles=ssles;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPCGSetup"
int DSDPCGSetup(DSDPCG *sles,DSDPVec X){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecGetSize(X,&sles->m);DSDPCHKERR(info);
  if (sles->setup2==0){
    info = DSDPVecDuplicate(X,&sles->R); DSDPCHKERR(info);
    info = DSDPVecDuplicate(X,&sles->P); DSDPCHKERR(info);
    info = DSDPVecDuplicate(X,&sles->BP); DSDPCHKERR(info);
    info = DSDPVecDuplicate(X,&sles->BR); DSDPCHKERR(info);
    info = DSDPVecDuplicate(X,&sles->Diag); DSDPCHKERR(info);
    info = DSDPVecDuplicate(X,&sles->TTT); DSDPCHKERR(info);
  }
  sles->setup2=1;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPCGDestroy"
int DSDPCGDestroy(DSDPCG **ssles){
  int info;
  DSDPCG *sles=*ssles;
  DSDPFunctionBegin;
  if (ssles==0 || sles==0){DSDPFunctionReturn(0);}
  if (sles->setup2==1){
    info = DSDPVecDestroy(&sles->R); DSDPCHKERR(info);
    info = DSDPVecDestroy(&sles->P); DSDPCHKERR(info);
    info = DSDPVecDestroy(&sles->BP); DSDPCHKERR(info);
    info = DSDPVecDestroy(&sles->BR); DSDPCHKERR(info);
    info = DSDPVecDestroy(&sles->Diag); DSDPCHKERR(info);
    info = DSDPVecDestroy(&sles->TTT); DSDPCHKERR(info);
  }
  DSDPFREE(ssles,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
