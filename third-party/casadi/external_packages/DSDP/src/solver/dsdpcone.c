#include "dsdpcone_impl.h"
#include "dsdpcone.h"
#include "dsdpsys.h"

/*!
\file dsdpcone.c
\brief Methods of a cone.
*/

#define DSDPNoOperationError(a);  { DSDPSETERR1(10,"Cone type: %s, Operation not defined\n",(a).dsdpops->name); }
#define DSDPChkConeError(a,b);  { if (b){DSDPSETERR1(b,"Cone type: %s,\n",(a).dsdpops->name); } }

/*!
\fn int DSDPConeSetUp(DSDPCone K, DSDPVec y);

\brief Factor the data and allocate data structures.
\param K the cone
\param y initial solution vector
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeSetUp"
int DSDPConeSetUp(DSDPCone K,DSDPVec y){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->conesetup){
    info=K.dsdpops->conesetup(K.conedata,y);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeSetUp2(DSDPCone K, DSDPVec yy0, DSDPSchurMat M);

\brief Factor the data and allocate data structures.
\param K the cone
\param yy0 initial solution vector
\param M Schur matrix
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeSetUp2"
int DSDPConeSetUp2(DSDPCone K, DSDPVec yy0, DSDPSchurMat M){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->conesetup2){
    info=K.dsdpops->conesetup2(K.conedata,yy0,M);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPConeDestroy(DSDPCone *K);

\brief Free the internal memory of the cone.
\param K the cone

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeDestroy"
int DSDPConeDestroy(DSDPCone *K){
  int info;
  DSDPFunctionBegin;
  if ((*K).dsdpops->conedestroy){
    info=(*K).dsdpops->conedestroy((*K).conedata);DSDPChkConeError(*K,info);
    info=DSDPConeInitialize(K); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(*K);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPConeComputeHessian(DSDPCone K, double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2);

\brief Compute Hessian and gradient of barrier function.
\param K the cone
\param mu barrier parameter
\param M Schur matrix
\param vrhs1 objective gradient
\param vrhs2 barrier gradient

This routine assumes that the dual matrix has already been factored and inverted.
\sa SDPConeComputeHessian()
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeComputeHessian"
int DSDPConeComputeHessian( DSDPCone K , double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->conehessian){
    info=K.dsdpops->conehessian(K.conedata,mu,M,vrhs1,vrhs2);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPConeMultiplyAdd(DSDPCone K, double mu, DSDPVec vrow, DSDPVec v, DSDPVec vv);

\brief Multiply Hessian by a vector and add the result.
\param K the cone
\param mu barrier parameter
\param vrow scaling for each element in the product.
\param v input vector gradient
\param vv output vector

This routine assumes that the dual matrix has already been factored and inverted.
If M is the hessian, then vv += vrow .* Mv
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeMultiplyAdd"
int DSDPConeMultiplyAdd( DSDPCone K , double mu, DSDPVec vrow, DSDPVec v, DSDPVec vv){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->conehmultiplyadd){
    info=K.dsdpops->conehmultiplyadd(K.conedata,mu,vrow,v,vv);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPConeComputeRHS(DSDPCone K, double mu, DSDPVec vrow, DSDPVec rhs1, DSDPVec rhs2);

\brief Compute gradient of barrier function.
\param K the cone
\param mu barrier parameter
\param vrow scaling for each element in the gradient.
\param rhs1 objective gradient
\param rhs2 barrier gradient

This routine assumes that the dual matrix has already been factored and inverted.
Define rhs2 += mu * vrow .* A(S^{-1})
\sa SDPConeComputeRHS()
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeComputeRHS"
int DSDPConeComputeRHS( DSDPCone K , double mu, DSDPVec vrow,DSDPVec rhs1,DSDPVec rhs2){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->conerhs){
    info=K.dsdpops->conerhs(K.conedata,mu,vrow,rhs1,rhs2);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeANorm2(DSDPCone K, DSDPVec anorm2);

\brief Add square of 2-norm of data correponding to each variable y.
\param K the cone
\param anorm2 norm of constraint data for each varibles
\sa DSDPBlockANorm2
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeANorm2"
int DSDPConeANorm2( DSDPCone K , DSDPVec anorm2){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->coneanorm2){
    info=K.dsdpops->coneanorm2(K.conedata,anorm2);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeSetXMaker(DSDPCone K, double mu, DSDPVec y, DSDPVec dy);

\brief Pass information needed to construct X.
\param K the cone
\param mu barrier parameter
\param y solution
\param dy step direction

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeSetXMaker"
int DSDPConeSetXMaker( DSDPCone K, double mu, DSDPVec y, DSDPVec dy){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->conesetxmaker){
    info=K.dsdpops->conesetxmaker(K.conedata,mu,y,dy);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeComputeX(DSDPCone K, double mu, DSDPVec y, DSDPVec dy, DSDPVec AX, double *tracexs);

\brief Given y,dy, and mu, construct X and add its inner product with the data and S
\param K the cone
\param mu barrier parameter
\param y solution
\param dy step direction
\param AX add the inner product of the data with X
\param tracexs inner product of X and S.
\sa SDPConeComputeXX()
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPConeComputeX"
int DSDPConeComputeX( DSDPCone K, double mu, DSDPVec y, DSDPVec dy, DSDPVec AX, double *tracexs){
  int info;
  double trxs;
  DSDPFunctionBegin;
  if (K.dsdpops->conecomputex){
    trxs=0;
    info=K.dsdpops->conecomputex(K.conedata,mu,y,dy,AX,&trxs);DSDPChkConeError(K,info);
    *tracexs+=trxs;
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeComputeS(DSDPCone K, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite);

\brief Given y, compute S and determine whether its in the cone.
\param K the cone
\param Y solution
\param flag identifies which of two S matrix structures should be used.
\param ispsdefinite true if S is positive definite or an element of the cone.

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeComputeS"
int DSDPConeComputeS(DSDPCone K, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->conecomputes){
    info=K.dsdpops->conecomputes(K.conedata,Y,flag,ispsdefinite);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPConeInvertS(DSDPCone K);

\brief Invert the dual matrix S.
\param K the cone

Assumes that the matrix has already been factored.

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeInvertS"
int DSDPConeInvertS(DSDPCone K){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->coneinverts){
    info=K.dsdpops->coneinverts(K.conedata);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeComputeMaxStepLength(DSDPCone K, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength);

\brief Determine distance to the edge of the cone.
\param K the cone
\param DY step direction
\param flag identifies which of two S matrix structures should be used.
\param maxsteplength distance to the edge of the cone.

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeComputeMaxStepLength"
int DSDPConeComputeMaxStepLength(DSDPCone K, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
  int info;
  double conesteplength=1.0e20;
  DSDPFunctionBegin;
  conesteplength=1.0e30;
  if (K.dsdpops->conemaxsteplength){
    info=K.dsdpops->conemaxsteplength(K.conedata,DY,flag,&conesteplength);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  *maxsteplength=conesteplength;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeGetDimension(DSDPCone K, double *n);

\brief Provide the dimension of the cone.
\param K the cone
\param n conic dimension (an integer value)

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeGetDimension"
int DSDPConeGetDimension(DSDPCone K, double *n){
  int info;
  double nn=0;
  DSDPFunctionBegin;
  if (K.dsdpops->conesize){
    info=K.dsdpops->conesize(K.conedata,&nn);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  *n=nn;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeSparsityInSchurMat(DSDPCone K, int row, int rnnz[], int m);

\brief Identify sparsity pattern in a row of the Hessian term
\param K the cone
\param row between 1 and m
\param rnnz mark elements nonzero for nonzeros in Hessian of barrier.
\param m number of y variables, length of array, and size of M matrix
\sa DSDPSparsityInSchurMat() 
\sa DSDPSchurSparsity()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSparsityInSchurMat"
int DSDPConeSparsityInSchurMat(DSDPCone K, int row, int rnnz[], int m){
  int info,tt;
  DSDPFunctionBegin;
  if (K.dsdpops->conesparsity){
    info=K.dsdpops->conesparsity(K.conedata,row,&tt,rnnz,m);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeView(DSDPCone K);

\brief View contents of the cone.
\param K the cone

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeView"
int DSDPConeView(DSDPCone K){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->coneview){
    info=K.dsdpops->coneview(K.conedata);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPConeMonitor(DSDPCone K, int tag);

\brief Do anything at in the cone at each iteration.
\param K the cone
\param tag allows for multiple types of monitors.

This routine has be used to visualize data, print some statistics, ...
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeMonitor"
int DSDPConeMonitor(DSDPCone K, int tag){
  int info;
  DSDPFunctionBegin;
  if (K.dsdpops->conemonitor){
    info=K.dsdpops->conemonitor(K.conedata,tag);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPConeComputeLogSDeterminant(DSDPCone K, double *logdetobj, double *logdet);

\brief Evaluate logrithmic barrier function.
\param K the cone
\param logdetobj used term.
\param logdet logarithmic barrier of cone
Assumes S is in cone.
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeComputeLogSDeterminant"
int DSDPConeComputeLogSDeterminant(DSDPCone K, double *logdetobj, double *logdet){
  int info;
  double conepotential1=0,conepotential2=0;
  DSDPFunctionBegin;
  if (K.dsdpops->conelogpotential){
      info=K.dsdpops->conelogpotential(K.conedata,&conepotential1,&conepotential2);DSDPChkConeError(K,info);
  } else {
    DSDPNoOperationError(K);
  }
  *logdetobj=conepotential1;
  *logdet=conepotential2;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetConeName(DSDPCone K, char *cname, int maxlength);

\brief Get name of the cone.
\param K the cone
\param cname string to copy the string
\param maxlength maximum length of the string.
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPGetConeName"
int DSDPGetConeName(DSDPCone K, char *cname, int maxlength){
  DSDPFunctionBegin;
  strncpy(cname,K.dsdpops->name,maxlength);
  DSDPFunctionReturn(0);
}



/*!
\fn int DSDPConeOpsInitialize(struct DSDPCone_Ops* dops);

\brief Initialize the function pointers to 0.
\param dops address of a structure of function pointers.
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeOpsInitialize"
int DSDPConeOpsInitialize(struct  DSDPCone_Ops* dops){
  DSDPFunctionBegin;
  if (dops==NULL) return 0;
  
  dops->conesetup=0;
  dops->conesetup2=0;
  dops->conedestroy=0;
  dops->coneanorm2=0;
  dops->conehessian=0;
  dops->conehmultiplyadd=0;
  dops->conerhs=0;
  dops->conesetxmaker=0;
  dops->conecomputex=0;
  dops->conecomputes=0;
  dops->coneinverts=0;
  dops->conemaxsteplength=0;
  dops->conesparsity=0;
  dops->conelogpotential=0;
  dops->conemonitor=0;
  dops->coneview=0;
  dops->id=0;
  DSDPFunctionReturn(0); 
}

/*!
\fn int DSDPConeSetData(DSDPCone *K, struct DSDPCone_Ops* ops,  void* data);

\brief Initialize the pointers to 0.
\param K the cone
\param ops address of a structure of function pointers.
\param data address of a structure representing a cone
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeSetData"
int DSDPConeSetData(DSDPCone *K, struct DSDPCone_Ops* ops,  void* data){
  DSDPFunctionBegin;
  (*K).dsdpops=ops;
  (*K).conedata=data;
  DSDPFunctionReturn(0); 
}


static struct  DSDPCone_Ops dsdpcops;
/*!
\fn int DSDPConeInitialize(DSDPCone *K);

\brief Initialize the pointers to 0.
\param K the cone

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPConeOpsInitialize"
int DSDPConeInitialize(DSDPCone *K){
  int info;
  DSDPFunctionBegin;
  info=DSDPConeOpsInitialize(&dsdpcops); DSDPCHKERR(info);
  info=DSDPConeSetData(K,&dsdpcops,0); DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

