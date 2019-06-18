#include "dsdpxmat_impl.h"
#include "dsdpxmat.h"
#include "dsdpsys.h"
/*! 
\file dsdpxmat.c
\brief Call an implementation of the basic dense matrix array operations.
 */

#define DSDPNoOperationError(a);  { DSDPSETERR1(1,"X Matrix type: %s, Operation not defined. Perhaps no X matrix has been set.\n",(a).dsdpops->matname); }
#define DSDPChkMatError(a,b);  { if (b){ DSDPSETERR1(b,"X Matrix type: %s,\n",(a).dsdpops->matname);} }

static int sdpxmatevent=0;

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatEventZero"
int DSDPVMatEventZero(void){
  DSDPFunctionBegin;
  sdpxmatevent=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatEventInitialize"
int DSDPVMatEventInitialize(void){
  DSDPFunctionBegin;
  if (sdpxmatevent==0){DSDPEventLogRegister("SDP X+vv'",&sdpxmatevent);}
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatSetData"
/*!
\fn int DSDPVMatSetData(DSDPVMat *X, struct DSDPVMat_Ops* ops,  void*data);
\brief Set opaque pointer an function pointers
\param X dense symmetric matrix
\param ops function pointers
\param data pointer to a matrix structure.
*/
int DSDPVMatSetData(DSDPVMat *X, struct DSDPVMat_Ops* ops,  void*data){
  int info;
  DSDPFunctionBegin;
  (*X).dsdpops=ops;
  (*X).matdata=data;
  info=DSDPVMatTest(*X);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatGetType"
int DSDPVMatGetType(DSDPVMat X, int *id){
  DSDPFunctionBegin;
  *id=X.dsdpops->id;
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatGetSize"
/*!
\fn int DSDPVMatGetSize(DSDPVMat X, int *n);
\brief Get number of rows and columns.
\param X dense symmetric matrix
\param n number of rows and columns
*/
int DSDPVMatGetSize(DSDPVMat X,int*n){
  int info;
  DSDPFunctionBegin;
  if (X.dsdpops->matgetsize){
    info=(X.dsdpops->matgetsize)(X.matdata,n); DSDPChkMatError(X,info);
  } else {
    /*
    DSDPNoOperationError(X);
    */
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatDestroy"
/*!
\fn int DSDPVMatDestroy(DSDPVMat *X);
\brief Deallocate matrix.
\param X dense symmetric matrix
*/
int DSDPVMatDestroy(DSDPVMat *X){
  int info;
  DSDPFunctionBegin;
  if (!(*X).dsdpops){ return 0;}
  if ((*X).dsdpops->matdestroy){
    info=((*X).dsdpops->matdestroy)((*X).matdata); DSDPChkMatError(*X,info);
  } else {
    /*  DSDPNoOperationError(*X); */
  }
  info=DSDPVMatInitialize(X); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatView"
/*!
\fn int DSDPVMatView(DSDPVMat X);
\brief Print matrix.
\param X dense symmetric matrix
*/
int DSDPVMatView(DSDPVMat X){
  int info;
  if (X.dsdpops->matview){
    info=(X.dsdpops->matview)(X.matdata); DSDPChkMatError(X,info);
  } else {
    printf("No viewer available for matrix type: %d",X.dsdpops->id);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatZeroEntries"
/*!
\fn int DSDPVMatZeroEntries(DSDPVMat X);
\brief Zero matrix.
\param X dense symmetric matrix
*/
int DSDPVMatZeroEntries(DSDPVMat X){
  int info;
  DSDPFunctionBegin;
  if (X.dsdpops->matzeroentries){
    info=(X.dsdpops->matzeroentries)(X.matdata); DSDPChkMatError(X,info);
  } else {
    DSDPNoOperationError(X);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatScaleDiagonal"
/*!
\fn int DSDPVMatScaleDiagonal(DSDPVMat X, double dscale);
\brief Scaling diagonal is useful for inner products and norms
\param X dense symmetric matrix
\param dscale
Semidefinite blocks scale the diagonal by half before taking the
dot product with the data matrices.
*/
int DSDPVMatScaleDiagonal(DSDPVMat X, double dscale){
  int info;
  DSDPFunctionBegin;
  if (X.dsdpops->matscalediagonal){
    info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); DSDPChkMatError(X,info);
  } else {
    DSDPNoOperationError(X);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatShiftDiagonal"
/*!
\fn int DSDPVMatShiftDiagonal(DSDPVMat X, double dadd);
\brief Add something to diagonal elements.
\param X dense symmetric matrix
\param dadd
*/
int DSDPVMatShiftDiagonal(DSDPVMat X, double dadd){
  int info;
  DSDPFunctionBegin;
  if (X.dsdpops->matshiftdiagonal){
    info=(X.dsdpops->matshiftdiagonal)(X.matdata,dadd); DSDPChkMatError(X,info);
  } else {
    DSDPNoOperationError(X);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatNormF2"
/*!
\fn int DSDPVMatNormF2(DSDPVMat X, double* normf2);
\brief Compute square of Frobenius norm of matrix.
\param X dense symmetric matrix
\param normf2 square of Frobenious norm.
*/
int DSDPVMatNormF2(DSDPVMat X, double*normf2){
  int info,n;
  double *dd;
  DSDPFunctionBegin;
  if (X.dsdpops->matfnorm2){
    info=DSDPVMatGetArray(X,&dd,&n); DSDPCHKERR(info);
    info=(X.dsdpops->matfnorm2)(X.matdata,n,normf2); DSDPChkMatError(X,info);
    info=DSDPVMatRestoreArray(X,&dd,&n); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(X);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatGetArray"
/*!
\fn int DSDPVMatGetArray(DSDPVMat X, double **v, int *nn);
\brief Get the array that stores the matrix.
\param X dense symmetric matrix
\param v array representing the matrix.
\param nn length of the array
\sa SDPConeSetStorageFormat()
*/
int DSDPVMatGetArray(DSDPVMat X, double **v, int *nn){
  int info;
  DSDPFunctionBegin;
  if (X.dsdpops->matgeturarray){
    info=(X.dsdpops->matgeturarray)(X.matdata,v,nn); DSDPChkMatError(X,info);
  } else {
    *v=0;
    *nn=0;
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatRestoreArray"
/*!
\fn int DSDPVMatRestoreArray(DSDPVMat X, double **v, int *nn);
\brief Restore the array that stores the matrix.
\param X dense symmetric matrix
\param v array representing the matrix.
\param nn length of the array
\sa SDPConeSetStorageFormat()
*/
int DSDPVMatRestoreArray(DSDPVMat X, double **v, int *nn){
  int info;
  DSDPFunctionBegin;
  if (X.dsdpops->matrestoreurarray){
    info=(X.dsdpops->matrestoreurarray)(X.matdata,v,nn); DSDPChkMatError(X,info);
  } else {
    *v=0;
    *nn=0;
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatMinEigenvalue"
int DSDPVMatMinEigenvalue(DSDPVMat X, SDPConeVec W1, SDPConeVec W2, double *mineig){
  int n,info;
  double *w,*iwork;
  DSDPFunctionBegin;
  if (X.dsdpops->matmineig){
    info=SDPConeVecGetSize(W1,&n); DSDPCHKERR(info);
    info=SDPConeVecGetArray(W1,&w); DSDPCHKERR(info);
    info=SDPConeVecGetArray(W2,&iwork); DSDPCHKERR(info);
    info=(X.dsdpops->matmineig)(X.matdata,w,iwork,n,mineig); DSDPChkMatError(X,info);
    info=SDPConeVecRestoreArray(W1,&w); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(W2,&iwork); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(X);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVMatAddOuterProduct"
/*!
\fn int DSDPVMatAddOuterProduct(DSDPVMat X, double alpha, SDPConeVec V);
\brief Add outer product of a vector to the matrix.
\param X dense symmetric matrix
\param alpha scalar multiple of outer product
\param V vector.
*/
int DSDPVMatAddOuterProduct(DSDPVMat X, double alpha, SDPConeVec V){
  int info,n;
  double *v;
  DSDPFunctionBegin;
  DSDPEventLogBegin(sdpxmatevent);
  info=SDPConeVecGetSize(V,&n); DSDPCHKERR(info);
  if (X.dsdpops->mataddouterproduct){
    info=SDPConeVecGetArray(V,&v); DSDPCHKERR(info);
    info=(X.dsdpops->mataddouterproduct)(X.matdata,alpha,v,n); DSDPChkMatError(X,info);
    info=SDPConeVecRestoreArray(V,&v); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(X);
  }
  DSDPEventLogEnd(sdpxmatevent);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatMult"
/*!
\fn int DSDPVMatMult(DSDPVMat X, SDPConeVec Z, SDPConeVec Y);
\brief Multiply X by a vector.
\param X dense symmetric matrix
\param Z input vector
\param Y equals X * Z
*/
int DSDPVMatMult(DSDPVMat X, SDPConeVec Z, SDPConeVec Y){
  int info,n;
  double *z,*y;
  DSDPFunctionBegin;
  info=SDPConeVecGetSize(Y,&n); DSDPCHKERR(info);
  if (X.dsdpops->matmult){
    info=SDPConeVecGetArray(Z,&z); DSDPCHKERR(info);
    info=SDPConeVecGetArray(Y,&y); DSDPCHKERR(info);
    info=(X.dsdpops->matmult)(X.matdata,z,y,n); DSDPChkMatError(X,info);
    info=SDPConeVecRestoreArray(Z,&z); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(Y,&y); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(X);
  }
  DSDPFunctionReturn(0);
}
 
#undef __FUNCT__
#define __FUNCT__ "DSDPVMatCheck"
/*!
\fn int DSDPVMatCheck(DSDPVMat X, SDPConeVec W1, SDPConeVec W2);
\brief Test correctness of operations.
\param X dense symmetric matrix
\param W1 work vector
\param W2 work vector
*/
int DSDPVMatCheck(DSDPVMat X, SDPConeVec W1, SDPConeVec W2){
  int info,i,n,nn;
  double *xx,eig,eps=1e-13,one=1.0;
  double fnorm0,fnorm1,fnorm2,fnorm3,fnorm4;
  DSDPFunctionBegin;
  info=DSDPVMatGetSize(X,&n);DSDPCHKERR(info);
  info=SDPConeVecSet(one,W1);DSDPCHKERR(info);
  info=DSDPVMatAddOuterProduct(X,1.23456,W1);DSDPCHKERR(info);
  info=DSDPVMatZeroEntries(X);DSDPCHKERR(info);
  info=DSDPVMatNormF2(X,&fnorm0);DSDPCHKERR(info);
  if (fabs(fnorm0)>eps){ printf("Check DSDPVMatZero of DSDPVMatNorm\n");}
  
  info=SDPConeVecSet(one,W1);DSDPCHKERR(info);
  info=DSDPVMatAddOuterProduct(X,1.0,W1);DSDPCHKERR(info);
  info=DSDPVMatNormF2(X,&fnorm1);DSDPCHKERR(info);
  if (fabs(fnorm1-n*n)>eps) printf("Check DSDPVMatZero()\n");

  info=DSDPVMatGetArray(X,&xx,&nn);DSDPCHKERR(info);
  for (i=0;i<nn;i++){xx[i]=1.0;}
  info=DSDPVMatRestoreArray(X,&xx,&nn);DSDPCHKERR(info);
  info=DSDPVMatNormF2(X,&fnorm2);DSDPCHKERR(info);
  if (fabs(fnorm2-n*n)>eps) printf("Check DSDPXGetArray()\n");

  info=DSDPVMatAddOuterProduct(X,-1.0,W1);DSDPCHKERR(info);
  info=DSDPVMatNormF2(X,&fnorm3);DSDPCHKERR(info);

  info=DSDPVMatZeroEntries(X);DSDPCHKERR(info);
  info=DSDPVMatAddOuterProduct(X,1.0,W1);DSDPCHKERR(info);
  info=DSDPVMatScaleDiagonal(X,2.0);DSDPCHKERR(info);

  info=DSDPVMatZeroEntries(X);DSDPCHKERR(info);
  info=DSDPVMatAddOuterProduct(X,1.0,W1);DSDPCHKERR(info);
  info=DSDPVMatShiftDiagonal(X,1.0);DSDPCHKERR(info);
  info=DSDPVMatNormF2(X,&fnorm4);DSDPCHKERR(info);

  info=DSDPVMatMult(X,W1,W2);DSDPCHKERR(info);
  info=DSDPVMatMinEigenvalue(X,W1,W2,&eig);DSDPCHKERR(info);
  if (fabs(fnorm0)>eps) printf("Check DSDPVMatZero()\n");
  DSDPFunctionReturn(0);
}

static struct  DSDPVMat_Ops dsdpmatops2;
static const char *urmatname="NOT SET YET";
#undef __FUNCT__
#define __FUNCT__ "DSDPVMatOpsInitialize"
/*!
\fn int DSDPVMatOpsInitialize(struct  DSDPVMat_Ops*aops);
\brief Set function pointers to null
\param aops pointer to a structure of function pointers.
*/
int DSDPVMatOpsInitialize(struct  DSDPVMat_Ops*aops){
  aops->matgetsize=0;
  aops->matzeroentries=0;
  aops->matfnorm2=0;
  aops->mataddouterproduct=0;
  aops->matmult=0;
  aops->matgeturarray=0;
  aops->matrestoreurarray=0;
  aops->matview=0;
  aops->matdestroy=0;
  aops->matmineig=0;
  aops->matshiftdiagonal=0;
  aops->matscalediagonal=0;
  aops->id=0;
  aops->matname=urmatname;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatTest"
/*!
\fn int DSDPVMatTest(DSDPVMat X);
\brief Test validity of matrix
\param X dense symmetric matrix
*/
int DSDPVMatTest(DSDPVMat X){
  int info;
  DSDPFunctionBegin;
  if (X.dsdpops==0 || X.dsdpops==&dsdpmatops2){
  } else if (X.dsdpops->mattest){
    info=(X.dsdpops->mattest)(X.matdata); DSDPChkMatError(X,info);
  } else {
    /*
    DSDPNoOperationError(X);
    */
  }
  DSDPFunctionReturn(0);
}
 

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatInitialize"
/*!
\fn int DSDPVMatInitialize(DSDPVMat *B);
\brief Set pointers to null
\param B dense symmetric matrix
*/
int DSDPVMatInitialize(DSDPVMat *B){
  int info;
  DSDPFunctionBegin;
  info=DSDPVMatOpsInitialize(&dsdpmatops2); DSDPCHKERR(info);
  info=DSDPVMatSetData(B, &dsdpmatops2, 0); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVMatExist"
/*!
\fn int DSDPVMatExist(DSDPVMat X, int *flag);
\brief Answer whether the array has been allocated or not.
\param X dense symmetric matrix
\param flag true if the array has been allocated
*/
int DSDPVMatExist(DSDPVMat X,int *flag){
  DSDPFunctionBegin;
  if (X.dsdpops && X.dsdpops!=&dsdpmatops2) *flag=1;
  else *flag=0;
  DSDPFunctionReturn(0);
}

