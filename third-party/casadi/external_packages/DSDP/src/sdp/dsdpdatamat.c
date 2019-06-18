#include "dsdpdatamat_impl.h"
#include "dsdpdatamat.h"
#include "dsdpsys.h"
/*!
\file dsdpdatamat.c
\brief Call an implementation of the data matrix operations.
*/
#define DSDPNoOperationError(a);  { DSDPSETERR1(1,"Data natrix type: %s, Operation not defined\n",(a).dsdpops->matname); }
#define DSDPChkDataError(a,b);  { if (b){ DSDPSETERR1(b,"Data natrix type: %s,\n",(a).dsdpops->matname);} }


static struct  DSDPDataMat_Ops dsdpdatamatdefault;

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatSetData"
/*!
\fn int DSDPDataMatSetData(DSDPDataMat *A, struct DSDPDataMat_Ops* ops,  void*data);

\brief Set the opaque pointer and function pointers to the matrix.
\param A symmetric data matrix
\param ops pointer to a structure of function pointers
\param data pointer to a matrix structure

*/
int DSDPDataMatSetData(DSDPDataMat *A, struct DSDPDataMat_Ops* ops,  void*data){
  int info;
  DSDPFunctionBegin;
  (*A).dsdpops=ops;
  (*A).matdata=data;
  if (ops==NULL){
    (*A).dsdpops=&dsdpdatamatdefault;   
  }
  info = DSDPDataMatOpsInitialize(&dsdpdatamatdefault); DSDPCHKERR(info);
  info=DSDPDataMatTest(*A);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

static char datamatnoname[20]="MATRIX NOT SET";
#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatOpsInitialize"
/*!
\fn int DSDPDataMatOpsInitialize(struct  DSDPDataMat_Ops* dops);
\brief Initialize the table of function pointers for SDP Data matrices.

\param dops pointer to a structure of function pointers.
*/
int DSDPDataMatOpsInitialize(struct  DSDPDataMat_Ops* dops){
  DSDPFunctionBegin;
  if (dops==NULL) return 0;
  dops->matfactor1=0;
  dops->matfactor2=0;
  dops->matgetrank=0;
  dops->matgeteig=0;
  dops->matvecvec=0;
  dops->matdot=0;
  dops->mataddrowmultiple=0;
  dops->mataddallmultiple=0;
  dops->matdestroy=0;
  dops->matview=0;
  dops->matrownz=0;
  dops->matnnz=0;
  dops->matfnorm2=0;
  dops->id=0;
  dops->matname=datamatnoname;
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatInitialize" 
/*!
\fn int DSDPDataMatInitialize(DSDPDataMat *A);

\brief Set pointers to NULL;
\param A symmetric data matrix

*/
int DSDPDataMatInitialize(DSDPDataMat *A){
  int info;
  DSDPFunctionBegin;
  info = DSDPDataMatOpsInitialize(&dsdpdatamatdefault); DSDPCHKERR(info);
  info = DSDPDataMatSetData(A, &dsdpdatamatdefault,0); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatTest"
/*!
\fn int DSDPDataMatTest(DSDPDataMat A);

\brief Test validity of matrix
\param A symmetric data matrix

*/
int DSDPDataMatTest(DSDPDataMat A){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops==0 || A.dsdpops==&dsdpdatamatdefault){
  } else if (A.dsdpops->mattest){
    info=(A.dsdpops->mattest)(A.matdata); DSDPChkDataError(A,info);
  } else {
    /*
    DSDPNoOperationError(A);
    */
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatGetType"
int DSDPDataMatGetType(DSDPDataMat A, int *id){
  DSDPFunctionBegin;
  *id=A.dsdpops->id;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatGetRank"
/*!
\fn int DSDPDataMatGetRank(DSDPDataMat A, int *rank, int n);

\brief Get the number of nonzero eigenvalues/eigenvectors for the matrix
\param A symmetric data matrix
\param rank number of nonzero eigenvalues and vectors.
\param n rows and columns in matrix.

*/
int DSDPDataMatGetRank(DSDPDataMat A, int *rank, int n){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops->matgetrank){
    info=(A.dsdpops->matgetrank)(A.matdata,rank,n); DSDPChkDataError(A,info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatCountNonzeros"
/*!
\fn int DSDPDataMatCountNonzeros(DSDPDataMat A, int *nnz, int n);

\brief Compute the square of the Frobenius norm 
\param A symmetric data matrix.
\param n dimension of the matrix.
\param *nnz nonzeros in matrix.

Used to identify which of a few strategies to compute Hessian.
*/
int DSDPDataMatCountNonzeros(DSDPDataMat A, int *nnz, int n){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops->matnnz){
    info=(A.dsdpops->matnnz)(A.matdata,nnz,n); DSDPChkDataError(A,info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatFNorm2"
/*!
\fn int DSDPDataMatFNorm2(DSDPDataMat A, int n, double *fnorm2);

\brief Compute the square of the Frobenius norm 
\param A symmetric data matrix.
\param n dimension of the matrix.
\param fnorm2 square of norm.

Used to scale the problem
*/
int DSDPDataMatFNorm2(DSDPDataMat A, int n, double *fnorm2){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops->matfnorm2){
    *fnorm2=0.0;
    info=(A.dsdpops->matfnorm2)(A.matdata,n,fnorm2); DSDPChkDataError(A,info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatGetEig"
/*!
\fn int DSDPDataMatGetEig(DSDPDataMat A, int rr, SDPConeVec V, DSDPIndex S, double *eigenvalue);

\brief Get an eigenvalue/vector pair.
\param A symmetric data matrix
\param rr identifies which pair of eigs,  (0 <= rr < rank)
\param V Eigenvector
\param S Indentifies sparsity pattern in V.
\param eigenvalue the scalar associated with the vector.
These pairs do not have to be eigenvalues and eigenvectors.  What matters
is that the matrix is a sum of the outer products of these vectors.
That is, A = sum (rr * V *V') 
\sa DSDPDataMatGetRank()
*/
int DSDPDataMatGetEig(DSDPDataMat A, int rr, SDPConeVec V, DSDPIndex S, double *eigenvalue){
  int info,n;
  double *vv;
  DSDPFunctionBegin;
  if (A.dsdpops->matgeteig){
    info=SDPConeVecGetArray(V,&vv); DSDPCHKERR(info);
    info=SDPConeVecGetSize(V,&n); DSDPCHKERR(info);
    info=(A.dsdpops->matgeteig)(A.matdata,rr, eigenvalue, vv,n,S.indx+1,S.indx); DSDPChkDataError(A,info);
    info=SDPConeVecRestoreArray(V,&vv); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatFactor"
/*!
\fn int DSDPDataMatFactor(DSDPDataMat A, SDPConeVec W, double*dworknn, int nn0, double *dwork3n, int nd, int* iwork, int ni){


\brief Do eigenvalue/vector or other factorization.
\param A symmetric data matrix
\param W work vector
\param dworknn work array
\param nn0 length of dworknn array.
\param dwork3n work array
\param nd length of work array
\param iwork work array
\param ni length of iwork array

This routine is called once during the DSDPSetUp routine.

\sa DSDPDataMatGetRank()
\sa DSDPDataMatGetEig()
*/
int DSDPDataMatFactor(DSDPDataMat A, SDPConeVec W, double*dworknn, int nn0, 
		      double *dwork3n, int nd, int* iwork, int ni){
  int info,n;
  double *dvecwork;
  DSDPFunctionBegin;
  if (A.dsdpops->matfactor1){
    info=(A.dsdpops->matfactor1)(A.matdata); DSDPChkDataError(A,info);
  } else if (A.dsdpops->matfactor2){
    info=SDPConeVecGetSize(W,&n);
    info=SDPConeVecGetArray(W,&dvecwork);
    info=(A.dsdpops->matfactor2)(A.matdata,dworknn,nn0,dvecwork,n,dwork3n,nd,iwork,ni); DSDPChkDataError(A,info);
    info=SDPConeVecRestoreArray(W,&dvecwork);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0);   
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatDot"
/*!
\fn int DSDPDataMatDot(DSDPDataMat A, double x[], int nn, int n, double *v);

\brief Compute inner product of data with a dense matrix
\param A symmetric data matrix
\param x dense array matrix
\param nn length of array
\param n dimension of matrix.
\param v the inner product

\sa SDPConeSetStorageFormat()
\sa DSDPVMatGetArray()
*/
int DSDPDataMatDot(DSDPDataMat A, double x[], int nn, int n, double *v){
  int info;

  DSDPFunctionBegin;
  if (A.dsdpops->matdot){
    info=(A.dsdpops->matdot)(A.matdata,x,nn,n,v); DSDPChkDataError(A,info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatVecVec"
/*!
\fn int DSDPDataMatVecVec(DSDPDataMat A, SDPConeVec W, double *v);

\brief Compute w' A w
\param A symmetric data matrix
\param W vector
\param v the inner product


*/
int DSDPDataMatVecVec(DSDPDataMat A, SDPConeVec W, double *v){
  int info,n;
  double *x;

  DSDPFunctionBegin;
  if (A.dsdpops->matvecvec){
    info=SDPConeVecGetSize(W,&n); DSDPCHKERR(info);
    info=SDPConeVecGetArray(W,&x); DSDPCHKERR(info);
    info=(A.dsdpops->matvecvec)(A.matdata,x,n,v); DSDPChkDataError(A,info);
    info=SDPConeVecRestoreArray(W,&x); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatMultiply"
/*
\fn int DSDPDataMatMultiply(DSDPDataMat A, SDPConeVec V1, SDPConeVec V2);

\brief Compute V2 = A*V1;
\param A symmetric data matrix
\param V1 in vector
\param V2 the product
Not needed.
*/
int DSDPDataMatMultiply(DSDPDataMat A, SDPConeVec V1, SDPConeVec V2){
  int info,n;
  double *vv1,*vv2;

  DSDPFunctionBegin;
  if (A.dsdpops->matmultiply){
    info=SDPConeVecGetSize(V1,&n); DSDPCHKERR(info);
    info=SDPConeVecGetArray(V1,&vv1); DSDPCHKERR(info);
    info=SDPConeVecGetArray(V2,&vv2); DSDPCHKERR(info);
    info=(A.dsdpops->matmultiply)(A.matdata,vv1,vv2,n); DSDPChkDataError(A,info);
    info=SDPConeVecRestoreArray(V1,&vv1); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(V2,&vv2); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatGetRowNonzeros"
/*!
\fn int DSDPDataMatGetRowNonzeros(DSDPDataMat A, int nrow, int nmax, int *nz, int *nnz);

\brief Get sparsity pattern of a row of the matrix.
\param A symmetric data matrix.
\param nrow number >=0 and < n.
\param nmax dimension of the matrix.
\param nz array used to mark nonzeros.
\param nnz number of nonzeros in row.
Used to create sparse data structure for S and Delta S.
*/
int DSDPDataMatGetRowNonzeros(DSDPDataMat A, int nrow, int nmax, int *nz, int *nnz){
  int i,info;
  DSDPFunctionBegin;
  if (A.dsdpops->matrownz){
    info=(A.dsdpops->matrownz)(A.matdata,nrow,nz,nnz,nmax); DSDPChkDataError(A,info);
  } else {
    *nnz=nmax;
    for (i=0;i<nmax;i++){
      nz[i]++;
    }
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatAddRowMultipleToVector"
int DSDPDataMatAddRowMultipleToVector(DSDPDataMat A, int nrow, double ytmp, SDPConeVec R){
  int info,n;
  double *vv;
  DSDPFunctionBegin;
  if (A.dsdpops->mataddrowmultiple){
    info=SDPConeVecGetArray(R,&vv);DSDPCHKERR(info);
    info=SDPConeVecGetSize(R,&n);DSDPCHKERR(info);
    info=(A.dsdpops->mataddrowmultiple)(A.matdata,nrow,ytmp,vv,n); DSDPChkDataError(A,info);
    info=SDPConeVecRestoreArray(R,&vv);DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatAddMultiple"
/*!
\fn int DSDPDataMatAddMultiple(DSDPDataMat A, double ytmp, double *v, int nn, int n);

\brief Add a multiple the data matrix to the array.
\param A symmetric data matrix
\param ytmp scalar multiple
\param v dense array matrix
\param nn dimension of array
\param n dimension of matrix

\sa SDPConeSetStorageFormat()
\sa DSDPVMatGetArray()
*/
int DSDPDataMatAddMultiple(DSDPDataMat A, double ytmp, double *v, int nn, int n){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops->mataddallmultiple){
    info=(A.dsdpops->mataddallmultiple)(A.matdata,ytmp,v,nn,n); DSDPChkDataError(A,info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatView"
/*!
\fn int DSDPDataMatView(DSDPDataMat A);

\brief Print matrix
\param A symmetric data matrix

*/
int DSDPDataMatView(DSDPDataMat A){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops->matview){
    info=(A.dsdpops->matview)(A.matdata); DSDPChkDataError(A,info);
  } else {
    printf("No matrix view available for matrix type %s.\n",A.dsdpops->matname);
  }
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatDestroy"
/*!
\fn int DSDPDataMatDestroy(DSDPDataMat *A);

\brief Free the data structures
\param A symmetric data matrix

*/
int DSDPDataMatDestroy(DSDPDataMat* A){
  int info;
  DSDPFunctionBegin;
  if ( (*A).dsdpops->matdestroy){
    info=((*A).dsdpops->matdestroy)((*A).matdata); DSDPChkDataError(*A,info);
  } else {
    /*   DSDPNoOperationError(*A); */
  }
  info=DSDPDataMatInitialize(A); DSDPCHKERR(info);
  /*  info=DSDPZeroMatCreate(0,A); DSDPCHKERR(info); */

  DSDPFunctionReturn(0); 
}

