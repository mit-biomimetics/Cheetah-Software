#include "dsdpdsmat_impl.h"
#include "dsdpdsmat.h"
#include "dsdpsys.h"
/*! \file dsdpdsmat.c
\brief Call an implmentation of the Delta S matrix operation 
*/

#define DSDPNoOperationError(a);  { DSDPSETERR1(1,"Delta S Matrix type: %s, Operation not defined\n",(a).dsdpops->matname); }
#define DSDPChkMatError(a,b);  { if (b){ DSDPSETERR1(b,"Delta S Matrix type: %s,\n",(a).dsdpops->matname); } }

#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatGetType"
int DSDPDSMatGetType(DSDPDSMat A, int *id){
  DSDPFunctionBegin;
  *id=A.dsdpops->id;
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatSetData"
/*!
\fn int DSDPDSMatSetData(DSDPDSMat *M, struct DSDPDSMat_Ops* ops,  void*data);

\brief Set the opaque pointer and function pointers to the matrix.
\param M symmetric DS matrix
\param ops pointer to a structure of function pointers
\param data pointer to a matrix structure

*/
int DSDPDSMatSetData(DSDPDSMat *M, struct DSDPDSMat_Ops* ops,  void*data){
  int info;
  DSDPFunctionBegin;
  (*M).dsdpops=ops;
  (*M).matdata=data;
  info=DSDPDSMatTest(*M); DSDPChkMatError(*M,info);
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatGetSize"
/*!
\fn int DSDPDSMatGetSize(DSDPDSMat A,int *n);

\brief Set the opaque pointer and function pointers to the matrix.
\param A symmetric DS matrix
\param n dimension
*/
int DSDPDSMatGetSize(DSDPDSMat A,int*n){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops->matgetsize){
    info=(A.dsdpops->matgetsize)(A.matdata,n); DSDPChkMatError(A,info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatDestroy"
/*!
\fn int DSDPDSMatDestroy(DSDPDSMat *A);

\brief Free the data structure.
\param A symmetric DS matrix
*/
int DSDPDSMatDestroy(DSDPDSMat *A){
  int info;
  DSDPFunctionBegin;
  if (!(*A).dsdpops){ return 0;}
  if ((*A).dsdpops->matdestroy){
    info=((*A).dsdpops->matdestroy)((*A).matdata); DSDPChkMatError(*A,info);
  } else {
    /*    DSDPNoOperationError(1); */
  }
  info=DSDPDSMatInitialize(A); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatView"
/*!
\fn int DSDPDSMatView(DSDPDSMat A);

\brief Print the matrix
\param A symmetric DS matrix
*/
int DSDPDSMatView(DSDPDSMat A){
  int info;
  if (A.dsdpops->matview){
    info=(A.dsdpops->matview)(A.matdata); DSDPChkMatError(A,info);
  } else {
    printf("No viewer available for matrix type: %s",A.dsdpops->matname);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatZeroEntries"
/*!
\fn int DSDPDSMatZeroEntries(DSDPDSMat A);
\brief Zero the entries in the matrix.
\param A symmetric DS matrix
*/
int DSDPDSMatZeroEntries(DSDPDSMat A){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops->matzeroentries){
    info=(A.dsdpops->matzeroentries)(A.matdata); DSDPChkMatError(A,info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatSetArray"
/*!
\fn int DSDPDSMatSetArray(DSDPDSMat A, DSDPVMat T);
\brief Set values into the matrix.
\param A symmetric DS matrix
\param T Source of entries in dense format
\sa DSDPSetFormatType()
*/
int DSDPDSMatSetArray(DSDPDSMat A, DSDPVMat T){
  int info,n,nn;
  double *ds;
  DSDPFunctionBegin;
  if (A.dsdpops->matseturmat){
    info=DSDPVMatGetSize(T,&n);DSDPCHKERR(info);
    info=DSDPVMatGetArray(T, &ds, &nn); DSDPCHKERR(info);
    info=(A.dsdpops->matseturmat)(A.matdata,ds,nn,n); DSDPChkMatError(A,info);
    info=DSDPVMatRestoreArray(T, &ds, &nn); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(A);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatMult"
/*!
\fn int DSDPDSMatMult(DSDPDSMat A, SDPConeVec X, SDPConeVec Y);
\brief Set values into the matrix.
\param A symmetric DS matrix
\param X in vector
\param Y product of A and X
*/
int DSDPDSMatMult(DSDPDSMat A, SDPConeVec X, SDPConeVec Y){
  int n,info;
  double *x,*y;
  
  DSDPFunctionBegin;
  if (A.dsdpops->matmult){
    info=SDPConeVecGetArray(X,&x); DSDPCHKERR(info);
    info=SDPConeVecGetArray(Y,&y); DSDPCHKERR(info);
    info=SDPConeVecGetSize(Y,&n); DSDPCHKERR(info);
    info=(A.dsdpops->matmult)(A.matdata,x,y,n); DSDPChkMatError(A,info);
    info=SDPConeVecRestoreArray(X,&x); DSDPCHKERR(info);
    info=SDPConeVecRestoreArray(Y,&y); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(A);
   }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDSVecVec"
/*!
\fn int DSDPDSMatVecVec(DSDPDSMat A, SDPConeVec X, double *vAv);
\brief Compute the product x' A x.
\param A symmetric DS matrix
\param X vector
\param vAv the product
*/
int DSDPDSMatVecVec(DSDPDSMat A, SDPConeVec X, double *vAv){
  int n,info;
  double *x;
  
  DSDPFunctionBegin;
  if (A.dsdpops->matvecvec){
    info=SDPConeVecGetArray(X,&x); DSDPCHKERR(info);
    info=SDPConeVecGetSize(X,&n); DSDPCHKERR(info);
    info=(A.dsdpops->matvecvec)(A.matdata,x,n,vAv); DSDPChkMatError(A,info);
    info=SDPConeVecRestoreArray(X,&x); DSDPCHKERR(info);
  } else {
    DSDPNoOperationError(A);
   }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatCheck"
int DSDPDSMatCheck(DSDPDSMat DS,SDPConeVec W1,SDPConeVec W2,DSDPVMat T){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


static struct  DSDPDSMat_Ops dsdpmatops2;
static const char* dsmatname="NOT SET YET";
#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatOpsInitialize"
/*!
\fn int DSDPDSMatOpsInitialize(struct  DSDPDSMat_Ops*aops);
\brief Set pointers to null.
\param aops pointer to table.
*/
int DSDPDSMatOpsInitialize(struct  DSDPDSMat_Ops*aops){
  aops->matseturmat=0;
  aops->matview=0;
  aops->matdestroy=0;
  aops->matgetsize=0;
  aops->matzeroentries=0;
  aops->matmult=0;
  aops->mattest=0;
  aops->matvecvec=0;
  aops->id=0;
  aops->matname=dsmatname;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatTest"
int DSDPDSMatTest(DSDPDSMat A){
  int info;
  DSDPFunctionBegin;
  if (A.dsdpops==0 || A.dsdpops==&dsdpmatops2){
  } else if (A.dsdpops->mattest){
    DSDPLogInfo(0,120,"Start to set DS Matrix\n");
    info=(A.dsdpops->mattest)(A.matdata); DSDPChkMatError(A,info);
    DSDPLogInfo(0,120,"Done set DS Matrix\n");
  } else {
    /*
    DSDPNoOperationError(A);
    */
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDSMatInitialize"
/*!
\fn int DSDPDSMatInitialize(DSDPDSMat *B);
\brief Set pointers to null.
\param B pointer to matrix.
*/
int DSDPDSMatInitialize(DSDPDSMat *B){
  int info;
  DSDPFunctionBegin;
  info=DSDPDSMatOpsInitialize(&dsdpmatops2); DSDPCHKERR(info);
  info=DSDPDSMatSetData(B, &dsdpmatops2, 0); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

