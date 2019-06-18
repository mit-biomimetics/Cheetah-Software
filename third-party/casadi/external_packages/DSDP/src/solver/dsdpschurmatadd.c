#include "dsdpschurmat_impl.h"
#include "dsdpschurmat.h"
#include "dsdpbasictypes.h"
#include "dsdpsys.h"

/*!
\file dsdpschurmatadd.c
\brief Cones to assemble the Schur matrix with these routines
*/

#define DSDPNoOperationError(a);  { DSDPSETERR1(10,"Schur matrix type: %s, Operation not defined\n",(a).dsdpops->matname); }
#define DSDPChkMatError(a,b);  { if (b){ DSDPSETERR1(b,"Schur matrix type: %s,\n",(a).dsdpops->matname);} }


#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatRowColumnScaling"
/*!
\fn int DSDPSchurMatRowColumnScaling(DSDPSchurMat M, int row, DSDPVec V, int *nzcols);
\brief Get the scaling and nonzero pattern of each column in this row of the matrix.
\param M matrix
\param row correponding to variable y
\param V multply each element of the row by this quantity.
\param nzcols how many nonzeros. Check for a 0!
Conic object call this routine when evaluating the Hessian of the barrier
term.  The vector V identifies sparsity, whether its using upper
or lower half of Schur, and also used to distribute rows over processros
The elements will be a 0 or a 1.  This routine is to be used with
DSDPSchurMatAddRow().

\sa DSDPSchurMatAddRow()

*/
int DSDPSchurMatRowColumnScaling(DSDPSchurMat M,int row, DSDPVec V, int *nzcols){
  int info,m;
  double *cols,r=M.schur->r;
  DSDPTruth flag;
  DSDPFunctionBegin;
  info=DSDPVecSet(0.0,V);DSDPCHKERR(info);
  info=DSDPVecGetSize(V,&m);DSDPCHKERR(info);
  if (row==0){info=DSDPVecZero(V);DSDPCHKERR(info);*nzcols=0;}
  else if (row==m-1){
    info=DSDPVecZero(V);DSDPCHKERR(info);*nzcols=0;
    if (r){info=DSDPVecSetR(V,1.0);DSDPCHKERR(info);*nzcols=1;}
  } else if (M.dsdpops->matrownonzeros){
    info=DSDPVecGetSize(V,&m);DSDPCHKERR(info);
    info=DSDPVecGetArray(V,&cols);DSDPCHKERR(info);
    info=(M.dsdpops->matrownonzeros)(M.data,row-1,cols+1,nzcols,m-2); DSDPChkMatError(M,info);
    info=DSDPVecRestoreArray(V,&cols);DSDPCHKERR(info);
    info=DSDPZeroFixedVariables(M,V);DSDPCHKERR(info);
    info=DSDPVecSetC(V,0.0);DSDPCHKERR(info);
    if (r){info=DSDPVecSetR(V,1.0);DSDPCHKERR(info);}
    info=DSDPIsFixed(M,row,&flag);DSDPCHKERR(info); 
    if (flag==DSDP_TRUE&&*nzcols>0){info=DSDPVecZero(V);*nzcols=0;DSDPFunctionReturn(0);}
  } else {
    DSDPNoOperationError(M);
  }

  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatAddRow"
/*!
\fn int DSDPSchurMatAddRow(DSDPSchurMat M, int row, double alpha, DSDPVec R);
\brief Add elements to a row of the Schur matrix.
\param M matrix
\param row correponding to variable y
\param alpha multiply elements in R by this scalar.
\param R a row of elements.

Conic object call this routine when evaluating the Hessian of the barrier
term.  This routine is to be used with DSDPSchurMatRowColumnScaling().

\sa DSDPSchurMatRowColumnScaling()
*/
int DSDPSchurMatAddRow(DSDPSchurMat M, int row, double alpha, DSDPVec R){
  int info,j,m;
  double *v,rr,dd=1e-1*M.schur->dd;
  DSDPVec rhs3=M.schur->rhs3;
  DSDPTruth flag;
  DSDPFunctionBegin;
  info=DSDPVecGetSize(R,&m); DSDPCHKERR(info);
  if (row==0){
  } else if (row==m-1){
    info=DSDPVecGetR(R,&rr);DSDPCHKERR(info);
    info=DSDPVecAddR(rhs3,alpha*rr);DSDPCHKERR(info);
  } else if (M.dsdpops->mataddrow){
    info=DSDPVecGetArray(R,&v); DSDPCHKERR(info);
    /*    v[row]=DSDPMax(0,v[row]); v[row]+=1.0e-15; */
    for (j=0;j<m;j++){ if (fabs(v[j]) < 1e-25 && row!=j){ v[j]=0.0;} }
    v[row]*=(1.0+dd);
    info=DSDPZeroFixedVariables(M,R);DSDPCHKERR(info);
    info=DSDPIsFixed(M,row,&flag);DSDPCHKERR(info); 
    if (flag==DSDP_TRUE){info=DSDPVecSetBasis(R,row);DSDPCHKERR(info);}
    info=(M.dsdpops->mataddrow)(M.data,row-1,alpha,v+1,m-2); DSDPChkMatError(M,info);
    info=DSDPVecRestoreArray(R,&v); DSDPCHKERR(info);  
    info=DSDPVecGetR(R,&rr); DSDPCHKERR(info);
    info=DSDPVecAddElement(rhs3,row,alpha*rr); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(M);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatVariableCompute"
/*!
\fn int DSDPSchurMatVariableCompute(DSDPSchurMat M, int row, double *rcv);

\brief Determine with the cone should compute this diagonal element of M and RHS.

\param M matrix
\param row correponding the y variable
\param rcv

Used to evaluate M. Important in parallel implementation.
\sa DSDPSchurMatAddDiagonalElement()

*/
int DSDPSchurMatVariableCompute(DSDPSchurMat M, int row, double *rcv){
  int info,m,dd=1;
  double r=M.schur->r;
  DSDPTruth flag;
  DSDPFunctionBegin;
  info=DSDPVecGetSize(M.schur->rhs3,&m);
  if (row==0){ *rcv=0.0;
  } else if (row==m-1){ 
    if (r){*rcv=1.0;}
    else {*rcv=0.0;}
  } else if (M.dsdpops->pmatonprocessor){
    info=(M.dsdpops->pmatonprocessor)(M.data,row-1,&dd); DSDPChkMatError(M,info);
    if (dd){*rcv=1.0;} else {*rcv=0;}
  } else {
    info=DSDPSchurMatInParallel(M,&flag);DSDPChkMatError(M,info);
    if (flag==DSDP_FALSE){ *rcv=1.0;
    } else {
      DSDPNoOperationError(M);
    }
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatAddDiagonalElement"
/*!
\fn int DSDPSchurMatAddDiagonalElement(DSDPSchurMat M, int row, double dd);

\brief Determine with the cone should compute this diagonal element of M and RHS.

\param M matrix
\param row correponding the y variable
\param dd zero or one..

\sa DSDPSchurMatVariableCompute()

*/
int DSDPSchurMatAddDiagonalElement(DSDPSchurMat M, int row, double dd){
  int info,m;
  DSDPVec rhs3=M.schur->rhs3;
  DSDPFunctionBegin;
  info=DSDPVecGetSize(rhs3,&m);DSDPCHKERR(info);
  if (dd==0){
  } else if (row==0){  
  } else if (row==m-1){info=DSDPVecAddR(rhs3,dd);
  } else if (M.dsdpops->mataddelement){
    info=(M.dsdpops->mataddelement)(M.data,row-1,dd); DSDPChkMatError(M,info);
  } else {
    DSDPNoOperationError(M);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatVariableComputeR"
/*!
\fn int DSDPSchurMatVariableComputeR(DSDPSchurMat M, double *rcv);
\brief Add an element to the Schur matrix correponding the variable r.
\param M matrix
\param *rcv zero or one
*/
int DSDPSchurMatVariableComputeR(DSDPSchurMat M, double *rcv){
  DSDPFunctionBegin;
  *rcv=0.0;
  if (M.schur->r) *rcv=1.0;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatAddR"
/*!
\fn int DSDPSchurMatAddR(DSDPSchurMat M, int row, double dd);
\brief Add an element to the Schur matrix correponding the variable r.
\param M matrix
\param row corresponding to variable y.
\param dd element
*/
int DSDPSchurMatAddR(DSDPSchurMat M, int row, double dd){
  int info;
  DSDPFunctionBegin;
  if (dd==0){DSDPFunctionReturn(0);}
  info=DSDPVecAddElement(M.schur->rhs3,row,dd);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatVariableComputeC"
int DSDPSchurMatVariableComputeC(DSDPSchurMat M, double *rcv){
  DSDPFunctionBegin;
  *rcv=0.0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatAddC"
int DSDPSchurMatAddC(DSDPSchurMat M, int row, double dd){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatDiagonalScaling"
/*!
\fn int DSDPSchurMatDiagonalScaling(DSDPSchurMat M, DSDPVec D);

\brief Get the scaling and nonzero pattern of each diagonal element of the matrix.
\param M matrix
\param D multply each element of the diagonal by this quantity.

Conic object call this routine when evaluating the Hessian of the barrier
term. The elements will be a 0 or a 1. Important for parallel version.

\sa DSDPSchurMatAddDiagonal()
*/
int DSDPSchurMatDiagonalScaling(DSDPSchurMat M, DSDPVec D){
  int m,info;
  double *vars;
  DSDPTruth flag;
  DSDPFunctionBegin;
  info=DSDPVecSet(1.0,D);DSDPCHKERR(info);
  info=DSDPVecGetSize(D,&m);DSDPCHKERR(info);
  if (M.dsdpops->pmatlocalvariables){
    info=DSDPVecGetArray(D,&vars);DSDPCHKERR(info);
    info=(M.dsdpops->pmatlocalvariables)(M.data,vars+1,m-2); DSDPChkMatError(M,info);
    info=DSDPVecRestoreArray(D,&vars);DSDPCHKERR(info);
  } else {
    info=DSDPSchurMatInParallel(M,&flag);DSDPChkMatError(M,info);
    if (flag==DSDP_TRUE){
      DSDPNoOperationError(M);
    }
  }
  info=DSDPVecSetC(D,0.0);DSDPCHKERR(info);
  if (M.schur->r==0){info=DSDPVecSetR(D,0.0);DSDPCHKERR(info);}
  info=DSDPZeroFixedVariables(M,D);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSchurMatAddDiagonal"
/*!
\fn int DSDPSchurMatAddDiagonal(DSDPSchurMat M, DSDPVec D);

\brief Add elements to a row of the Schur matrix.
\param M matrix
\param D diagonal elements.

Conic object call this routine when evaluating the Hessian of the barrier
term.  

\sa DSDPSchurMatDiagonalScaling()
*/
int DSDPSchurMatAddDiagonal(DSDPSchurMat M, DSDPVec D){
  int m,info;
  double *dd;
  DSDPFunctionBegin;
  if (M.dsdpops->matadddiagonal){
    info=DSDPZeroFixedVariables(M,D);DSDPCHKERR(info);
    info=DSDPVecGetSize(D,&m); DSDPCHKERR(info);
    info=DSDPVecGetArray(D,&dd); DSDPCHKERR(info);
    info=(M.dsdpops->matadddiagonal)(M.data,dd+1,m-2); DSDPChkMatError(M,info);
    info=DSDPVecAddR(M.schur->rhs3,dd[m-1]);DSDPCHKERR(info);
    info=DSDPVecRestoreArray(D,&dd); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(M);
  }
  DSDPFunctionReturn(0);
}


