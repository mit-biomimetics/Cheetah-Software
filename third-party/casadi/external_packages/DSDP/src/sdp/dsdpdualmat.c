#include "dsdpdualmat_impl.h"
#include "dsdpdualmat.h"
#include "dsdpsys.h"

/*!
\file dsdpdualmat.c
\brief Call an implementation of the S matrix operations.
*/

#define DSDPNoOperationError(a);  { DSDPSETERR1(1,"Dual natrix type: %s, Operation not defined\n",(a).dsdpops->matname);}
#define DSDPChkDMatError(a,b);  { if (b){ DSDPSETERR1(b,"Dual natrix type: %s,\n",(a).dsdpops->matname);} }

static int sdpdualsolve=0,sdpdualinvert=0;

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatEventZero"
int DSDPDualMatEventZero(void){
  DSDPFunctionBegin;
  sdpdualinvert=0;sdpdualsolve=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatEventInitialize"
int DSDPDualMatEventInitialize(void){
  DSDPFunctionBegin;
  if (sdpdualsolve==0){DSDPEventLogRegister("SDP SSolve",&sdpdualsolve);}
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatGetType"
int DSDPDualMatGetType(DSDPDualMat S, int *id){
  DSDPFunctionBegin;
  *id=S.dsdpops->id;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatSetData"
/*!
\fn int DSDPDualMatSetData(DSDPDualMat *S, struct DSDPDualMat_Ops* ops,  void*data);

\brief Set the opaque pointer and function pointers to the matrix.
\param S dual matrix
\param ops pointer to a structure of function pointers
\param data pointer to a matrix structure
*/
int DSDPDualMatSetData(DSDPDualMat *S, struct DSDPDualMat_Ops* ops,  void*data){
  int info;
  DSDPFunctionBegin;
  (*S).dsdpops=ops;
  (*S).matdata=data;
  info=DSDPDualMatTest(*S);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatDestroy"
/*!
\fn int DSDPDualMatDestroy(DSDPDualMat *S);
\brief Free the matrix structure
\param S dual matrix
*/
int DSDPDualMatDestroy(DSDPDualMat *S){
  int info;
  DSDPFunctionBegin;
  if ( S && (*S).dsdpops && (*S).dsdpops->matdestroy){
    info=((*S).dsdpops->matdestroy)((*S).matdata); DSDPChkDMatError(*S,info);
  } else {
    /*
    DSDPNoOperationError(*S);
    */
  }
  info=DSDPDualMatSetData(S,0,0); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatGetSize"
/*!
\fn int DSDPDualMatGetSize(DSDPDualMat S, int *n);
\brief Free the matrix structure
\param S dual matrix
\param n dimension
*/
int DSDPDualMatGetSize(DSDPDualMat S,int*n){
  int info;
  DSDPFunctionBegin;
  if (S.dsdpops->matgetsize){
    info=(S.dsdpops->matgetsize)(S.matdata,n); DSDPChkDMatError(S,info);
  } else {
    DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatGetArray"
int DSDPDualMatGetArray(DSDPDualMat S, double **v, int *n){
  int info;
  DSDPFunctionBegin;
  if (S.dsdpops->matgetarray){
    info=(S.dsdpops->matgetarray)(S.matdata,v,n); DSDPChkDMatError(S,info);
  } else {
    *v=0;
    *n=0;
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatLogDeterminant"
/*!
\fn int DSDPDualMatLogDeterminant(DSDPDualMat S, double *logdet);
\brief Free the matrix structure
\param S dual matrix
\param logdet logarithm of the determinant
Assumes Cholesky factorization was successful.
*/
int DSDPDualMatLogDeterminant(DSDPDualMat S,double *logdet){
  int info;
  DSDPFunctionBegin;
  if (S.dsdpops->matlogdet){
    info=(S.dsdpops->matlogdet)(S.matdata,logdet); DSDPChkDMatError(S,info);
  } else {
     DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatView"
/*!
\fn int DSDPDualMatView(DSDPDualMat S);
\brief Print the matrix.
\param S dual matrix
*/
int DSDPDualMatView(DSDPDualMat S){
  int info;
  DSDPFunctionBegin;
  if (S.dsdpops->matview){
    info=(S.dsdpops->matview)(S.matdata); DSDPChkDMatError(S,info);
  } else {
     DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatSetArray"
/*!
\fn int DSDPDualMatSetArray(DSDPDualMat S, DSDPVMat T);
\brief Print the matrix.
\param S dual matrix
\param T Dense array matrix.
*/
int DSDPDualMatSetArray(DSDPDualMat S, DSDPVMat T){
  double *ss;
  int info,n,nn;
  DSDPFunctionBegin;
  if (S.dsdpops->matseturmat){
    info=DSDPVMatGetSize(T,&n); DSDPCHKERR(info);
    info=DSDPVMatGetArray(T,&ss,&nn); DSDPCHKERR(info);
    info=(S.dsdpops->matseturmat)(S.matdata,ss,nn,n); DSDPChkDMatError(S,info);
    info=DSDPVMatRestoreArray(T,&ss,&nn); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatInvert"
/*!
\fn int DSDPDualMatInvert(DSDPDualMat S);
\brief Invert the matrix.
\param S dual matrix
Assumes Cholesky factorization was successful.  This routine
may not actually invert the matrix. It give the matrix the
opportunity to invert it.
*/
int DSDPDualMatInvert(DSDPDualMat S){
  int info;
  DSDPFunctionBegin;
  /*  DSDPEventLogBegin(sdpdualinvert); */
  if (S.dsdpops->matinvert){
    info=(S.dsdpops->matinvert)(S.matdata); DSDPChkDMatError(S,info);
  } else {
    DSDPNoOperationError(S);
  }
  /*  DSDPEventLogEnd(sdpdualinvert); */
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatInverseAdd"
/*!
\fn int DSDPDualMatInverseAdd(DSDPDualMat S, double alpha, DSDPVMat T);
\brief Add a multiple of the inverse to T.
\param S dual matrix
\param alpha scalar
\param T destination.
Assumes matrix already inverted.
*/
int DSDPDualMatInverseAdd(DSDPDualMat S, double alpha, DSDPVMat T){
  int info,n,nn;
  double *ss;
  DSDPFunctionBegin;
  if (S.dsdpops->matinverseadd){
    info=DSDPVMatGetSize(T,&n); DSDPCHKERR(info);
    info=DSDPVMatGetArray(T,&ss,&nn); DSDPCHKERR(info);
    info=(S.dsdpops->matinverseadd)(S.matdata,alpha,ss,nn,n); DSDPChkDMatError(S,info);
    info=DSDPVMatRestoreArray(T,&ss,&nn); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatInverseMultiply"
/*!
\fn int DSDPDualMatInverseMultiply(DSDPDualMat S, DSDPIndex IS, SDPConeVec B, SDPConeVec X);
\brief Multiply the inverse by a vector or solve the system of equations.
\param S dual matrix
\param IS Sparsity pattern of B
\param B Right-hand side of linear system
\param X product, or solution to linear system.
Assumes matrix already inverted.
*/
int DSDPDualMatInverseMultiply(DSDPDualMat S, DSDPIndex IS, SDPConeVec B, SDPConeVec X){
  int info,n;
  double *bb,*xx;
  DSDPFunctionBegin;
  DSDPEventLogBegin(sdpdualsolve);
  if (S.dsdpops->matinversemultiply){
    info=SDPConeVecGetSize(X,&n); DSDPCHKERR(info);
    info=SDPConeVecGetArray(B,&bb); DSDPCHKERR(info);
    info=SDPConeVecGetArray(X,&xx); DSDPCHKERR(info);
    info=(S.dsdpops->matinversemultiply)(S.matdata,IS.indx+1,IS.indx[0],bb,xx,n); DSDPChkDMatError(S,info);
    info=SDPConeVecRestoreArray(X,&xx); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(B,&bb); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(S);
  }
  DSDPEventLogEnd(sdpdualsolve);
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatCholeskySolveForward"
/*!
\fn int DSDPDualMatCholeskySolveForward(DSDPDualMat S, SDPConeVec B, SDPConeVec X);

\brief Forward triangular solve.
\param S dual matrix
\param B Right-hand side of triangular system
\param X Solution to triangular system.
Assumes Cholesky factorization successful.
*/
int DSDPDualMatCholeskySolveForward(DSDPDualMat S, SDPConeVec B, SDPConeVec X){
  int info,n;
  double *bb,*xx;
  DSDPFunctionBegin;
  if (S.dsdpops->matsolveforward){
    info=SDPConeVecGetSize(X,&n); DSDPCHKERR(info);
    info=SDPConeVecGetArray(B,&bb); DSDPCHKERR(info);
    info=SDPConeVecGetArray(X,&xx); DSDPCHKERR(info);
    info=(S.dsdpops->matsolveforward)(S.matdata,bb,xx,n); DSDPChkDMatError(S,info);
    info=SDPConeVecRestoreArray(X,&xx); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(B,&bb); DSDPCHKERR(info);
  } else {
     DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatDualMatCholeskySolveBackward"
/*!
\fn int DSDPDualMatCholeskySolveBackward(DSDPDualMat S, SDPConeVec B, SDPConeVec X);

\brief Backward triangular solve.
\param S dual matrix
\param B Right-hand side of triangular system
\param X Solution to triangular system.
Assumes Cholesky factorization successful.
*/
int DSDPDualMatCholeskySolveBackward(DSDPDualMat S, SDPConeVec B, SDPConeVec X){
  int info,n;
  double *bb,*xx;
  DSDPFunctionBegin;
  if (S.dsdpops->matsolvebackward){
    info=SDPConeVecGetSize(X,&n); DSDPCHKERR(info);
    info=SDPConeVecGetArray(B,&bb); DSDPCHKERR(info);
    info=SDPConeVecGetArray(X,&xx); DSDPCHKERR(info);
    info=(S.dsdpops->matsolvebackward)(S.matdata,bb,xx,n); DSDPChkDMatError(S,info);
    info=SDPConeVecRestoreArray(X,&xx); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(B,&bb); DSDPCHKERR(info);
  } else {
     DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatCholeskyFactor"
/*!
\fn int DSDPDualMatCholeskyFactor(DSDPDualMat S, DSDPTruth *psdefinite);
\brief Factor the matrix
\param S dual matrix
\param psdefinite true if S is positive definite and factorization successful.
*/
int DSDPDualMatCholeskyFactor(DSDPDualMat S,DSDPTruth *psdefinite){
  int info;
  int flag;
  DSDPFunctionBegin;
  if (S.dsdpops->matcholesky){
    info=(S.dsdpops->matcholesky)(S.matdata,&flag); DSDPChkDMatError(S,info);
  } else {
     DSDPNoOperationError(S);
  }
  if (flag) *psdefinite=DSDP_FALSE;
  else  *psdefinite=DSDP_TRUE;
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatCholeskyForwardMultiply"
/*!
\fn int DSDPDualMatCholeskyForwardMultiply(DSDPDualMat S, SDPConeVec B, SDPConeVec X);

\brief Multiply by triangular matrix.
\param S dual matrix
\param B in vector
\param X product
Assumes Cholesky factorization successful.
*/
int DSDPDualMatCholeskyForwardMultiply(DSDPDualMat S, SDPConeVec B, SDPConeVec X){
  int info,n;
  double *bb,*xx;
  DSDPFunctionBegin;
  if (S.dsdpops->matforwardmultiply){
    info=SDPConeVecGetSize(B,&n); DSDPCHKERR(info);
    info=SDPConeVecGetArray(B,&bb); DSDPCHKERR(info);
    info=SDPConeVecGetArray(X,&xx); DSDPCHKERR(info);
    info=(S.dsdpops->matforwardmultiply)(S.matdata,bb,xx,n); DSDPChkDMatError(S,info);
    info=SDPConeVecRestoreArray(X,&xx); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(B,&bb); DSDPCHKERR(info);
  } else {
     DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0); 
}
#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatCholeskyBackwardMultiply"
/*!
\fn int DSDPDualMatCholeskyBackwardMultiply(DSDPDualMat S, SDPConeVec B, SDPConeVec X);

\brief Multiply by triangular matrix.
\param S dual matrix
\param B in vector
\param X product
Assumes Cholesky factorization successful.
*/
int DSDPDualMatCholeskyBackwardMultiply(DSDPDualMat S, SDPConeVec B, SDPConeVec X){
  int info,n;
  double *bb,*xx;
  DSDPFunctionBegin;
  if (S.dsdpops->matbackwardmultiply){
    info=SDPConeVecGetSize(B,&n); DSDPCHKERR(info);
    info=SDPConeVecGetArray(B,&bb); DSDPCHKERR(info);
    info=SDPConeVecGetArray(X,&xx); DSDPCHKERR(info);
    info=(S.dsdpops->matbackwardmultiply)(S.matdata,bb,xx,n); DSDPChkDMatError(S,info);
    info=SDPConeVecRestoreArray(X,&xx); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(B,&bb); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(S);
  }
  DSDPFunctionReturn(0); 
}
#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatIsFull"
/*!
\fn int DSDPDualMatIsFull(DSDPDualMat S, DSDPTruth *full);
\brief Factor the matrix
\param S dual matrix
\param full true if S is a dense structure.
*/
int DSDPDualMatIsFull(DSDPDualMat S, DSDPTruth *full){
  int info,flag=0;
  DSDPFunctionBegin;
  *full=DSDP_FALSE;
  if (S.dsdpops->matfull){
    info=(S.dsdpops->matfull)(S.matdata,&flag); DSDPChkDMatError(S,info);
  } else {
     DSDPNoOperationError(S);
  }
  if (flag) *full=DSDP_TRUE;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatCheck"
int DSDPDualMatCheck(DSDPDualMat SS, SDPConeVec W1, SDPConeVec W2, DSDPIndex IS, DSDPVMat XX){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0); 
}

static const char* dualmatname="NOT SET YET";
/*!
\fn int DSDPDualMatOpsInitialize(struct  DSDPDualMat_Ops* sops){
\brief Set pointers to null.
\param sops function pointers
*/
int DSDPDualMatOpsInitialize(struct  DSDPDualMat_Ops* sops){
  if (sops==NULL) return 0;
  sops->matseturmat=0;
  sops->matgetarray=0;
  sops->matcholesky=0;
  sops->matsolveforward=0;
  sops->matsolvebackward=0;
  sops->matinvert=0;
  sops->matinverseadd=0;
  sops->matinversemultiply=0;
  sops->matforwardmultiply=0;
  sops->matbackwardmultiply=0;
  sops->matfull=0;
  sops->matdestroy=0;
  sops->matgetsize=0;
  sops->matview=0;
  sops->matlogdet=0;
  sops->matname=dualmatname;
  return 0;
 }


static struct  DSDPDualMat_Ops dsdpdualmatopsdefault;

#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatTest"
int DSDPDualMatTest(DSDPDualMat S){
  int info;
  DSDPFunctionBegin;
  if (S.dsdpops==0 || S.dsdpops==&dsdpdualmatopsdefault){
  } else if (S.dsdpops->mattest){
    info=(S.dsdpops->mattest)(S.matdata); DSDPChkDMatError(S,info);
  } else {
    /*
     DSDPNoOperationError(S);
    */
  }
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDualMatInitialize"
/*!
\fn int DSDPDualMatInitialize(DSDPDualMat *S);
\brief Set pointers to null.
\param S dual matrix
*/
int DSDPDualMatInitialize(DSDPDualMat *S){
  int info;
  DSDPFunctionBegin;
  info=DSDPDualMatOpsInitialize(&dsdpdualmatopsdefault);DSDPCHKERR(info);
  info=DSDPDualMatSetData(S,&dsdpdualmatopsdefault,0); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

