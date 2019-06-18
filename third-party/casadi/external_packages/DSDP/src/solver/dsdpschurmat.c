#include "dsdpschurmat_impl.h"
#include "dsdpschurmat.h"
#include "dsdpbasictypes.h"
#include "dsdpsys.h"

/*!
\file dsdpschurmat.c
\brief Solve the Schur matrix for a step direction.
*/

static int hfactorevent=0,hsolveevent=0;

#define DSDPNoOperationError(a);  { DSDPSETERR1(10,"Schur matrix type: %s, Operation not defined\n",(a).dsdpops->matname); }
#define DSDPChkMatError(a,b);  { if (b){ DSDPSETERR1(b,"Schur matrix type: %s,\n",(a).dsdpops->matname);} }

static int DSDPApplySMW(DSDPSchurMat, DSDPVec, DSDPVec);
static int DSDPSchurMatSolveM(DSDPSchurMat, DSDPVec, DSDPVec);

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatSetData"
/*!
\fn int DSDPSchurMatSetData(DSDPSchurMat *M, struct DSDPSchurMat_Ops* ops,  void*data);
\brief Set the Schur matrix with an opaque pointer and structure of function pointers.
\param M Schur matrix object
\param ops address of a structure of function pointers.
\param data opaque pointer to its internal data structure.
*/
int DSDPSchurMatSetData(DSDPSchurMat *M, struct DSDPSchurMat_Ops* ops,  void*data){
  DSDPFunctionBegin;
  (*M).dsdpops=ops;
  (*M).data=data;
  DSDPFunctionReturn(0); 
}

static const char* schurmatname="NOT NAMED YET";

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatOpsInitialize"
/*!
\fn int DSDPSchurMatOpsInitialize(struct DSDPSchurMat_Ops* dops);
\brief Initialize function pointers to 0.
\param dops address of a structure of function pointers.
*/
int DSDPSchurMatOpsInitialize(struct  DSDPSchurMat_Ops* dops){
  DSDPFunctionBegin;
  if (dops==NULL) return 0;
  dops->matzero=0;
  dops->matrownonzeros=0;
  dops->mataddrow=0;
  dops->mataddelement=0;
  dops->matadddiagonal=0;
  dops->matshiftdiagonal=0;
  dops->matassemble=0;
  dops->matscaledmultiply=0;
  dops->matmultr=0;
  dops->matfactor=0;
  dops->matsolve=0;
  dops->pmatonprocessor=0;
  dops->pmatwhichdiag=0;
  dops->pmatdistributed=0;
  dops->matdestroy=0;
  dops->matview=0; 
  dops->matsetup=0;
  dops->id=0;
  dops->matname=schurmatname;
  DSDPFunctionReturn(0); 
}

static struct  DSDPSchurMat_Ops dsdpmops;


#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatOpsInitialize"
/*!
\fn int DSDPSchurMatInitialize(DSDPSchurMat *M);
\brief Initialize pointers to null.
\param M Schur matrix object
*/
int DSDPSchurMatInitialize(DSDPSchurMat *M){
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatOpsInitialize(&dsdpmops); DSDPCHKERR(info);
  info=DSDPSchurMatSetData(M,&dsdpmops,0); DSDPCHKERR(info);
  DSDPCALLOC1(&M->schur,DSDPSchurInfo,&info);DSDPCHKERR(info);
  M->schur->m=0;  M->schur->r=0;  M->schur->dd=0;
  info=DSDPInitializeFixedVariable(&M->schur->fv);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatZeroEntries"
/*!
\fn int DSDPSchurMatZeroEntries(DSDPSchurMat M);
\brief Zero all element in the matrix.
\param M matrix
*/
int DSDPSchurMatZeroEntries(DSDPSchurMat M){
  int info;
  DSDPFunctionBegin;
  if (M.dsdpops->matzero){
    info=(M.dsdpops->matzero)(M.data); DSDPChkMatError(M,info);
  } else {
    DSDPNoOperationError(M);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatShiftDiagonal"
/*!
\fn int DSDPSchurMatShiftDiagonal(DSDPSchurMat M, double dd);

\brief Add a scalar to each diagonal element of the matrix.

\param M matrix
\param dd diagonal shift

*/
int DSDPSchurMatShiftDiagonal(DSDPSchurMat M, double dd){
  int info;
  DSDPFunctionBegin;
  if (dd==0){DSDPFunctionReturn(0);}
  M.schur->dd=dd;
  if (M.dsdpops->matshiftdiagonal){
    /* if(M.schur->r){info=DSDPVecAddR(M.schur->rhs3,dd);DSDPCHKERR(info);} */
    info=(M.dsdpops->matshiftdiagonal)(M.data,dd); DSDPChkMatError(M,info);
    DSDPLogInfo(0,2,"Add %4.4e to the Diagonal of Schur Matrix\n",dd);
  } else {
    DSDPNoOperationError(M);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatInParallel"
/*!
\fn int DSDPSchurMatInParallel(DSDPSchurMat M, DSDPTruth *flag);

\brief Determine whether M is computed in parallel.

\param M matrix
\param flag true or false

Important in parallel implementation.

*/
int DSDPSchurMatInParallel(DSDPSchurMat M, DSDPTruth *flag){
  int info,flg;
  DSDPFunctionBegin;
  if (M.dsdpops->pmatdistributed){
    info=(M.dsdpops->pmatdistributed)(M.data,&flg); DSDPChkMatError(M,info);
    if (flg) *flag=DSDP_TRUE; else *flag=DSDP_FALSE;
  } else {
    *flag=DSDP_FALSE;
  }
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatAssemble"
/*!
\fn int DSDPSchurMatAssemble(DSDPSchurMat M);

\brief Final assembly of M.

\param M matrix

Important in parallel implementation.
*/
int DSDPSchurMatAssemble(DSDPSchurMat M){
  int info;
  DSDPFunctionBegin;
  if (M.dsdpops->matassemble){
    info=(M.dsdpops->matassemble)(M.data); DSDPChkMatError(M,info);
  } else {
    DSDPNoOperationError(M);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatFactor"
/*!
\fn int DSDPSchurMatFactor(DSDPSchurMat M, DSDPTruth *successful);

\brief Factor M

\param M matrix
\param successful false if factorization failed.

*/
int DSDPSchurMatFactor(DSDPSchurMat M, DSDPTruth *successful){
  int info,flag=0;
  DSDPVec rhs3=M.schur->rhs3,dy3=M.schur->dy3;
  DSDPFunctionBegin;
  *successful=DSDP_TRUE;
  DSDPEventLogBegin(hfactorevent);
  if (M.dsdpops->matfactor){
    info=(M.dsdpops->matfactor)(M.data,&flag); DSDPChkMatError(M,info);
    if (flag){ 
      *successful=DSDP_FALSE;
      DSDPLogInfo(0,2,"Indefinite Schur Matrix -- Bad Factorization\n");
    }
  } else {
    DSDPNoOperationError(M);
  }
  DSDPEventLogEnd(hfactorevent);
  if (M.schur->r){
    info=DSDPSchurMatSolveM(M,rhs3,dy3);DSDPCHKERR(info);}
  else {info=DSDPVecZero(dy3);DSDPCHKERR(info);}
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatMultiply"
/*!
\fn int DSDPSchurMatMultiply(DSDPSchurMat M, DSDPVec x, DSDPVec y);

\brief Multiply M by a vector. y = M x

\param M matrix
\param x in vector.
\param y product.

*/
int DSDPSchurMatMultiply(DSDPSchurMat M, DSDPVec x, DSDPVec y){
  int info,n;
  double *xx,*yy,r=M.schur->r;
  double r1,r2,dd;
  DSDPVec rhs3;
  DSDPFunctionBegin;

  if (M.dsdpops->matscaledmultiply){
    info=DSDPVecGetSize(x,&n); DSDPCHKERR(info);
    info=DSDPVecGetArray(x,&xx); DSDPCHKERR(info);
    info=DSDPVecGetArray(y,&yy); DSDPCHKERR(info);
    info=(M.dsdpops->matscaledmultiply)(M.data,xx+1,yy+1,n-2); DSDPChkMatError(M,info);
    yy[0]=0;
    yy[n-1]=0;
    info=DSDPVecRestoreArray(y,&yy); DSDPCHKERR(info);
    info=DSDPVecRestoreArray(x,&xx); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(M);
  }
  if (r){
    rhs3=M.schur->rhs3;
    info=DSDPVecGetR(rhs3,&r2);DSDPCHKERR(info);
    info=DSDPVecGetR(x,&r1);DSDPCHKERR(info);
    info=DSDPVecAXPY(r1,rhs3,y);DSDPCHKERR(info);
    info=DSDPVecDot(rhs3,x,&dd);DSDPCHKERR(info);
    info=DSDPVecAddR(y,dd-r1*r2);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatMultR"
int DSDPSchurMatMultR(DSDPSchurMat M, DSDPVec x, DSDPVec y){
  int info,n;
  double *xx,*yy,r=M.schur->r;
  double r1,r2,dd;
  DSDPVec rhs3;
  DSDPFunctionBegin;

  if (M.dsdpops->matmultr){
    info=DSDPVecGetSize(x,&n); DSDPCHKERR(info);
    info=DSDPVecGetArray(x,&xx); DSDPCHKERR(info);
    info=DSDPVecGetArray(y,&yy); DSDPCHKERR(info);
    info=(M.dsdpops->matmultr)(M.data,xx+1,yy+1,n-2); DSDPChkMatError(M,info);
    yy[0]=0;
    yy[n-1]=0;
    info=DSDPVecRestoreArray(y,&yy); DSDPCHKERR(info);
    info=DSDPVecRestoreArray(x,&xx); DSDPCHKERR(info);
    if (r){
      rhs3=M.schur->rhs3;
      info=DSDPVecGetR(rhs3,&r2);DSDPCHKERR(info);
      info=DSDPVecGetR(x,&r1);DSDPCHKERR(info);
      info=DSDPVecAXPY(r1,rhs3,y);DSDPCHKERR(info);
      info=DSDPVecDot(rhs3,x,&dd);DSDPCHKERR(info);
      info=DSDPVecAddR(y,dd-r1*r2);DSDPCHKERR(info);
    }
  } else {
    info=DSDPVecZero(y);DSDPCHKERR(info);
    /*    DSDPNoOperationError(M); */
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatReducePVec"
/*!
\fn int DSDPSchurMatReducePVec(DSDPSchurMat M, DSDPVec x);

\brief Collect elements of the vector

\param M matrix
\param x vector.

Important in parallel implementation. Otherwise does nothing.

*/
int DSDPSchurMatReducePVec(DSDPSchurMat M, DSDPVec x){
  int info,n;
  double *xx;
  DSDPTruth flag;
  DSDPFunctionBegin;

  if (M.dsdpops->pmatreduction){
    info=DSDPVecGetSize(x,&n); DSDPCHKERR(info);
    info=DSDPVecGetArray(x,&xx); DSDPCHKERR(info);
    info=(M.dsdpops->pmatreduction)(M.data,xx+1,n-2); DSDPChkMatError(M,info);
    info=DSDPVecRestoreArray(x,&xx); DSDPCHKERR(info);
  } else {
    info=DSDPSchurMatInParallel(M,&flag);DSDPChkMatError(M,info);
    if (flag==DSDP_TRUE){
      DSDPNoOperationError(M);
    }
  }
  info=DSDPZeroFixedVariables(M,x);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatSetR"
/*!
\fn int DSDPSchurMatSetR(DSDPSchurMat M, double rr);
\brief Set up the data structure
\param M matrix
\param rr residual
*/
int DSDPSchurMatSetR(DSDPSchurMat M, double rr){
  DSDPFunctionBegin;
  M.schur->r=rr;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatSetup"
/*!
\fn int DSDPSchurMatSetup(DSDPSchurMat M, DSDPVec Y);
\brief Set up the data structure
\param M matrix
\param Y variable vector.
*/
int DSDPSchurMatSetup(DSDPSchurMat M, DSDPVec Y){
  int info,m;
  DSDPFunctionBegin;
  info=DSDPVecDuplicate(Y,&M.schur->rhs3);
  info=DSDPVecDuplicate(Y,&M.schur->dy3);
  info=DSDPVecGetSize(Y,&m);DSDPCHKERR(info);
  if (M.dsdpops->matsetup){
    info=(M.dsdpops->matsetup)(M.data,m-2); DSDPChkMatError(M,info);
  } else {
    DSDPNoOperationError(M);
  }
  DSDPEventLogRegister("Factor Newton Eq.",&hfactorevent);
  DSDPEventLogRegister("Solve Newton Eq.",&hsolveevent);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatView"
/*!
\fn int DSDPSchurMatView(DSDPSchurMat M);
\brief Print the matrix
\param M matrix
*/
int DSDPSchurMatView(DSDPSchurMat M){
  int info;
  DSDPFunctionBegin;
  if (M.dsdpops->matview){
    info=(M.dsdpops->matview)(M.data); DSDPChkMatError(M,info);
  } else {
    DSDPNoOperationError(M);
  }
  info=DSDPVecView(M.schur->rhs3);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatRowScaling"
/*!
\fn int DSDPSchurMatRowScaling(DSDPSchurMat M, DSDPVec D);

\brief Identify which rows on on this processor.
\param M matrix
\param D Variables marked with 0 or 1.

*/
int DSDPSchurMatRowScaling(DSDPSchurMat M, DSDPVec D){
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatDiagonalScaling(M,D);DSDPCHKERR(info);
  info=DSDPZeroFixedVariables(M,D);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatDestroy"
/*!
\fn int DSDPSchurMatDestroy(DSDPSchurMat *M);
\brief Free the memory in the data structure
\param M matrix
*/
int DSDPSchurMatDestroy(DSDPSchurMat *M){
  int info;
  DSDPFunctionBegin;
  if ((*M).dsdpops->matdestroy){
    info=((*M).dsdpops->matdestroy)((*M).data); DSDPChkMatError(*M,info);
  } else {
    /*
    DSDPNoOperationError(*M);
    */
  }
  info=DSDPVecDestroy(&M->schur->rhs3);DSDPCHKERR(info);
  info=DSDPVecDestroy(&M->schur->dy3);DSDPCHKERR(info);
  /*  info=DSDPSchurMatSetData(M,0,0); DSDPCHKERR(info); */
  info=DSDPSchurMatOpsInitialize(&dsdpmops); DSDPCHKERR(info);
  info=DSDPSchurMatSetData(M,&dsdpmops,0); DSDPCHKERR(info);
  DSDPFREE(&M->schur,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatSolveM"
static int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x){
  int info,n;
  double *xx,*bb;
  DSDPFunctionBegin;
  info=DSDPEventLogBegin(hsolveevent);
  if (M.dsdpops->matsolve){
    info=DSDPVecGetArray(b,&bb); DSDPCHKERR(info);
    info=DSDPVecGetSize(x,&n); DSDPCHKERR(info);
    info=DSDPVecZero(x);DSDPCHKERR(info);
    info=DSDPVecGetArray(x,&xx); DSDPCHKERR(info);
    info=(M.dsdpops->matsolve)(M.data,bb+1,xx+1,n-2); DSDPChkMatError(M,info);
    info=DSDPVecRestoreArray(b,&bb); DSDPCHKERR(info);
    info=DSDPVecRestoreArray(x,&xx); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(M);
  }
  info=DSDPVecSetR(x,0.0);DSDPCHKERR(info);
  info=DSDPVecSetC(x,0.0);DSDPCHKERR(info);
  info=DSDPEventLogEnd(hsolveevent);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatSolve"
/*!
\fn int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x);
\brief Solve the linear system
\param M matrix
\param b the right-hand side
\param x solution
*/
int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
  info=DSDPApplySMW(M,b,x);DSDPCHKERR(info);
  info=DSDPZeroFixedVariables(M,x);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPApplySMW"
static int DSDPApplySMW(DSDPSchurMat M, DSDPVec rhs, DSDPVec dy){
  int info;
  double r=M.schur->r,rr,dr,rhsr,rssr;
  double rhsnorm,rhsnorm3,rhs1mrhs3=0,rhs3mrhs3=0;
  DSDPVec rhs3=M.schur->rhs3,dy3=M.schur->dy3;
  DSDPFunctionBegin;
  
  info=DSDPVecNormInfinity(rhs,&rhsnorm);DSDPCHKERR(info);
  info=DSDPVecNormInfinity(rhs3,&rhsnorm3);DSDPCHKERR(info);
  if (r==0 || rhsnorm==0){
    info=DSDPVecSetR(dy,0); DSDPCHKERR(info);
    info=DSDPVecSetR(rhs,0); DSDPCHKERR(info);
  } else if (0 && rhsnorm3==0){ /* dsdp->UsePenalty==DSDPNever */
    info=DSDPVecGetR(rhs,&rr); DSDPCHKERR(info);
    info=DSDPVecSetR(dy,rr); DSDPCHKERR(info);
  } else {
    /* Use bigM penalty method and Sherman-Morrison-Woodbury */
    info=DSDPVecGetR(rhs,&rhsr); DSDPCHKERR(info); 
    info=DSDPVecGetR(rhs3,&rssr); DSDPCHKERR(info); 
    info=DSDPVecDot(rhs3,dy,&rhs1mrhs3); DSDPCHKERR(info);
    info=DSDPVecDot(rhs3,dy3,&rhs3mrhs3); DSDPCHKERR(info);
    if (rssr-rhs3mrhs3==0) rssr*=(1.00001);
    dr=-(rhs1mrhs3-rhsr    )/(rssr-rhs3mrhs3);
    info=DSDPVecAXPY(-dr,dy3,dy);DSDPCHKERR(info);
    info=DSDPVecSetR(dy,dr); DSDPCHKERR(info);
    info=DSDPVecSetR(rhs,rhsr); DSDPCHKERR(info);
    info=DSDPVecDot(rhs,dy,&rhs3mrhs3); DSDPCHKERR(info);
    if (rhs3mrhs3 <=0){ 
      DSDPLogInfo(0,3,"DSDP Step Direction Not Descent, Adjusting. \n");
      info=DSDPVecAddR(rhs3,rssr*0.1);DSDPCHKERR(info);
      info=DSDPVecAXPY(dr,dy3,dy);DSDPCHKERR(info);
      info=DSDPVecSetR(dy,0); DSDPCHKERR(info);
      info=DSDPApplySMW(M,rhs,dy);DSDPCHKERR(info);
    }
  }   
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPZeroFixedVariables"
int DSDPZeroFixedVariables( DSDPSchurMat M, DSDPVec dy){
  int i,info; 
  FixedVariables *fv=&M.schur->fv;
  DSDPFunctionBegin;
  for (i=0;i<fv->nvars;i++){
    info=DSDPVecSetElement(dy,fv->var[i],0.0);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPApplyFixedVariables"
int DSDPApplyFixedVariables( DSDPSchurMat M, DSDPVec y){
  int i,jj,info;
  double vv,scl;
  FixedVariables *fv=&M.schur->fv;
  info=DSDPVecGetC(y,&scl);DSDPCHKERR(info);
  DSDPFunctionBegin;
  for (i=0;i<fv->nvars;i++){
    vv=fv->fval[i]*fabs(scl);
    jj=fv->var[i];
    info=DSDPVecSetElement(y,jj,vv);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPFixedVariableNorm"
int DSDPFixedVariablesNorm( DSDPSchurMat M, DSDPVec y){
  int i,jj,info;
  double vv;
  FixedVariables *fv=&M.schur->fv;
  DSDPFunctionBegin;
  for (i=0;i<fv->nvars;i++){
    jj=fv->var[i]; vv=fv->fval[i];
    info=DSDPVecAddC(y,1.0);DSDPCHKERR(info);
    info=DSDPVecAddElement(y,jj,vv*vv);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeFixedYX"
int DSDPComputeFixedYX( DSDPSchurMat M, DSDPVec berr){
  int i,jj,info;
  double vv;
  FixedVariables *fv=&M.schur->fv;
  DSDPFunctionBegin;
  for (i=0;i<fv->nvars;i++){
    jj=fv->var[i];
    info=DSDPVecGetElement(berr,jj,&vv);DSDPCHKERR(info);
    info=DSDPVecSetElement(berr,jj,0);DSDPCHKERR(info);
    info=DSDPVecAddC(berr,-vv*fv->fval[i]);DSDPCHKERR(info);
    info=DSDPVecAddR(berr,fabs(vv));DSDPCHKERR(info);
    fv->fdual[i]=-vv;
    if (fv->xout) fv->xout[i]=-vv;
    DSDPLogInfo(0,2,"FIXED VAR DUAL: %d %4.4f, ADD %4.4f to objective.\n",jj,vv,-vv*fv->fval[i]);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPIsFixed"
int DSDPIsFixed( DSDPSchurMat M, int vari, DSDPTruth *flag){
  int i;
  FixedVariables *fv=&M.schur->fv;
  DSDPFunctionBegin;
  *flag=DSDP_FALSE;
  for (i=0;i<fv->nvars;i++){
    if (fv->var[i]==vari){
      *flag=DSDP_TRUE;
      break;
    }
  }
  DSDPFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "DSDPInitializeFixedVariables"
int DSDPInitializeFixedVariable( FixedVariables *fv){
  DSDPFunctionBegin;
  fv->nmaxvars=0;
  fv->nvars=0;
  fv->fval=0;
  fv->var=0;
  fv->fdual=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPAddFixedVariables"
int DSDPAddFixedVariable( DSDPSchurMat M, int vari, double val){
  int i,t,*iinew,info,nvars;
  double *ddnew,*vvnew;
  FixedVariables *fv=&M.schur->fv;
  DSDPFunctionBegin;
  nvars=fv->nvars;
  if (nvars>=fv->nmaxvars){
    t=2*nvars + 2;
    DSDPCALLOC2(&iinew,int,t,&info);
    DSDPCALLOC2(&ddnew,double,t,&info);
    DSDPCALLOC2(&vvnew,double,t,&info);
    for (i=0;i<nvars;i++){
      iinew[i]=fv->var[i];
      ddnew[i]=fv->fval[i];
      vvnew[i]=fv->fdual[i];
    }
    DSDPFREE(&fv->var,&info);DSDPCHKERR(info);
    DSDPFREE(&fv->fval,&info);DSDPCHKERR(info);
    DSDPFREE(&fv->fdual,&info);DSDPCHKERR(info);
    fv->var=iinew;
    fv->fval=ddnew;
    fv->fdual=vvnew;
    fv->nmaxvars=t;
  }
  fv->var[fv->nvars]=vari;
  fv->fval[fv->nvars]=val;
  fv->nvars++;
  DSDPFunctionReturn(0);
}
 
#include "dsdp.h"

#undef __FUNCT__
#define __FUNCT__ "DSDPSparsityInSchurMat"
/*!
\fn int DSDPSparsityInSchurMat(DSDP dsdp, int row, int rnnz[], int mm);
\brief Identify nonzero elements in a row of the Schur complement.
\param dsdp solver
\param row corresponding to variable y.
\param rnnz array to be marked nonzero if nonzero.
\param mm dimension of M matrix
*/
int DSDPSparsityInSchurMat(DSDP dsdp, int row, int rnnz[], int mm){
  int info,*iptr,m=mm+2;
  double *dd;
  DSDPVec R=dsdp->M.schur->rhs3;
  DSDPFunctionBegin;
  info=DSDPVecZero(R);DSDPCHKERR(info);
  info=DSDPVecGetArray(R,&dd);DSDPCHKERR(info);
  iptr=(int*)dd;
  info=DSDPSchurSparsity(dsdp,row+1,iptr,m);DSDPCHKERR(info);
  memcpy((void*)rnnz,(void*)(iptr+1),(mm)*sizeof(int));
  info=DSDPVecRestoreArray(R,&dd);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#include "dsdp5.h"
#undef __FUNCT__
#define __FUNCT__ "DSDPSetFixedVariable"
/*!
\fn int DSDPSetFixedVariable(DSDP dsdp, int vari, double val);
\brief Fix variable y to exact value.
\param dsdp solver
\param vari variables y
\param val fixed value
\ingroup DSDPSolver
\sa DSDPSetFixedVariables()
 */
int DSDPSetFixedVariable(DSDP dsdp, int vari, double val){
   int info;
   DSDPFunctionBegin;
   DSDPLogInfo(0,2,"Set Fixed Variable: %d, %12.8f\n",vari,val);
   info= DSDPAddFixedVariable(dsdp->M,vari,val);DSDPCHKERR(info);
   DSDPFunctionReturn(0);
 }

#undef __FUNCT__
#define __FUNCT__ "DSDPSetFixedVariables"
/*!
\fn int DSDPSetFixedVariables(DSDP dsdp, double vars[], double vals[], double xout[], int nvars);
\brief Fix variable y to exact values.
\param dsdp solver
\param vars variables y ( integer valued from 1 through m)
\param vals fixed values
\param xout array for the dual variables
\param nvars length of the previous arrays.
\ingroup DSDPSolver
 */
int DSDPSetFixedVariables(DSDP dsdp, double vars[], double vals[], double xout[], int nvars){
   int i,info;
   DSDPFunctionBegin;
   for (i=0;i<nvars;i++){
     info=DSDPSetFixedVariable(dsdp,(int)vars[i],vals[i]);
     dsdp->M.schur->fv.xout=xout;
   }
   DSDPFunctionReturn(0);
 }

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetFixedYX"
int DSDPGetFixedYX( DSDP dsdp, int vari, double *dd){
  int i;
  FixedVariables *fv=&dsdp->M.schur->fv;
  DSDPFunctionBegin;
  for (i=0;i<fv->nvars;i++){
    if (vari==fv->var[i]){
      *dd=fv->fdual[i]; break;
    }
  }
  DSDPFunctionReturn(0);
}
