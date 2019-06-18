#include "dsdpdatamat.h"
#include "dsdpsys.h"
#include "dsdp5.h"

/*!
\file dsdpadddatamat.c
\brief Set sparse or dense matrices into the cone.
*/

extern int DSDPGetZeroDataMatOps(struct DSDPDataMat_Ops**);
extern int DSDPGetConstantMat(int,double,char,struct DSDPDataMat_Ops**,void**);

extern int DSDPGetVechMat(int,int,double,const int[], const double[],int, struct  DSDPDataMat_Ops**,void**);
extern int DSDPGetVecUMat(int,int,double,const int[], const double[],int, struct  DSDPDataMat_Ops**,void**);

extern int DSDPGetIdentityDataMatP(int,double,struct DSDPDataMat_Ops**,void**);
extern int DSDPGetIdentityDataMatF(int,double,struct DSDPDataMat_Ops**,void**);

extern int DSDPGetDMat(int,double,double[],struct  DSDPDataMat_Ops**,void**);

extern int DSDPGetR1PMat(int,double,int,const int[],const double[],int,struct  DSDPDataMat_Ops**,void**);
extern int DSDPGetR1UMat(int,double,int,const int[],const double[],int,struct  DSDPDataMat_Ops**,void**);

extern int SDPConeAddDataMatrix(SDPCone,int, int, int, char, struct DSDPDataMat_Ops*, void*);
extern int SDPConeSetRMatrix(SDPCone,int,int,char,struct DSDPDataMat_Ops*,void*);


#undef __FUNCT__
#define __FUNCT__ "SDPConeAddASparseVecMat"
/*!
\fn int SDPConeAddASparseVecMat(SDPCone sdpcone,int blockj, int vari, int n, double alpha, int ishift, const int ind[], const double val[], int nnz);
\brief Add data matrix  \f$A_{i,j}\f$ in a sparse format
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param ishift index of first element (usually 0)
\param alpha multiple of the data (usually 1)
\param val array of elements in the matrix
\param ind array of indices representing the location of the elements
\param nnz length of the previous two arrays
\note DSDP will use the two arrays throughout in many routines, but it will not modify or delete them until finished with DSDP.
\sa DSDPCreateSDPCone()
\sa SDPConeGetStorageFormat()
\sa SDPConeViewDataMatrix()
\sa SDPConeSetBlockSize()
*/
int SDPConeAddASparseVecMat(SDPCone sdpcone,int blockj, int vari, int n,
			    double alpha, int ishift, 
			    const int ind[], const double val[], int nnz){
  
  int info;
  char UPLQ;
  void* dmat=0;
  struct DSDPDataMat_Ops* dmatops=0;
  
  DSDPFunctionBegin;
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  DSDPLogInfo(0,20,"Set sparse matrix:  Block: %d, Variable %d, size: %d, Nonzeros: %d .\n",blockj,vari,n,nnz);
  switch (UPLQ){
  case 'P':
    info=DSDPGetVechMat(n,ishift,alpha,ind,val,nnz,&dmatops,&dmat); DSDPCHKERR(info);
    break;
  case 'U':
    info=DSDPGetVecUMat(n,ishift,alpha,ind,val,nnz,&dmatops,&dmat); DSDPCHKERR(info);
    break;
  }
  info=SDPConeAddDataMatrix(sdpcone,blockj,vari,n,UPLQ,dmatops,dmat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeAddSparseVecMat" 
/* Needed for backward compatibility */
int SDPConeAddSparseVecMat(SDPCone sdpcone,int blockj, int vari, int n,
			   int ishift,const int ind[], const double val[], int nnz){
  
  int info;
  
  DSDPFunctionBegin;
  info= SDPConeAddASparseVecMat(sdpcone,blockj,vari,n,
				1.0,ishift,ind,val,nnz);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeSetASparseVecMat"
/*!
\fn int SDPConeSetASparseVecMat(SDPCone sdpcone,int blockj, int vari, int n, double alpha, int ishift, const int ind[], const double val[], int nnz);
\brief Set data matrix  \f$A_{i,j}\f$ in a sparse format
\ingroup SDPBasic
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param ishift index of \f$ a_{1,1} \f$ (usually 0)
\param alpha multiple of the data (usually 1)
\param val array of elements in the matrix
\param ind array of indices representing the location of the elements
\param nnz length of the previous two arrays

For example, the matrix 
\image html img208.gif
can be inserted into the cone in packed symmetric format in several ways.

Using the ordering of the packed symmetric format, 
we can index each element of the matrix
with an integer between 0 and n(n+1)/2-1 (inclusive).
If the first element in the \c val array is \f$a_{1,1}\f$, 
the first element in the \c ind array should be 0 .  
If the second element in the \c val array is
\f$a_{2,1}\f$, then
the second element in \c ind array should be 1 
(the second element in dense representation).
If the third element in the \c val array is
\f$a_{3,2}\f$, then
the third element in \c ind array should be 4.
Explicitly,
\code
double val1[]={3,2,6};
int ind1[]={0,1,4};
SDPConeSetASparseVecMat(sdpcone,j,i,3,1.0,0,ind1,val1,3);
\endcode
If we index the elements from 1 through n(n+1)/2, we can use
\code
double val2[]={3,2,6};
int ind2[]={1,2,5};
SDPConeSetASparseVecMat(sdpcone,j,i,3,1.0,1,ind2,val2,3);
\endcode
The ordering of the elements in the array can be changed and zero elements
can be included.
\code
double val3[]={6,3,0,2};
int ind3[]={4,0,2,1};
SDPConeSetASparseVecMat(sdpcone,j,i,3,1.0,0,ind3,val3,4);
\endcode
The elements can also be scaled.
\code
double val1[]={6,4,12};
int ind1[]={0,1,4};
SDPConeSetASparseVecMat(sdpcone,j,i,3,0.5,0,ind1,val1,3);
\endcode

\note DSDP will use the two arrays throughout in many routines, but it will not modify or delete them until finished with DSDP.

\sa DSDPCreateSDPCone()
\sa DSDPSetDualObjective()
\sa SDPConeViewDataMatrix()

*/
int SDPConeSetASparseVecMat(SDPCone sdpcone,int blockj, int vari, int n,
			    double alpha,int ishift,
			    const int ind[], const double val[], int nnz){
  
  int info;
  DSDPFunctionBegin;
  info=SDPConeRemoveDataMatrix(sdpcone,blockj,vari); DSDPCHKERR(info);
  info=SDPConeAddASparseVecMat(sdpcone,blockj,vari,n,alpha,ishift,ind,val,nnz); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeSetSparseVecMat"
/* Needed for backward compatibility */
int SDPConeSetSparseVecMat(SDPCone sdpcone,int blockj, int vari, int n,
			   int ishift,const int ind[], const double val[], int nnz){
  
  int info;
  DSDPFunctionBegin;
  info=SDPConeSetASparseVecMat(sdpcone,blockj,vari,n,1.0,ishift,ind,val,nnz); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeAddADenseVecMat"
/*!
\fn int SDPConeAddADenseVecMat(SDPCone sdpcone,int blockj, int vari, int n, double alpha, double val[], int nnz);
\brief Add a matrix \f$A_{i,j}\f$ in a dense format.
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param alpha multiple of data (usually 1.0)
\param val array of elements in the matrix
\param nnz length of the previous two arrays
\note DSDP will use the \c val array in many routines, but it will not modify or delete it until finished with DSDP.

The matrix 
\image html img208.gif
can be inserted into the cone in symmetric packed format as follows
\code
double val[]={3,2,0,0,6,0};
SDPConeAddADenseVecMat(sdpcone,j,i,3,1.0,val,6);
\endcode

\sa SDPConeViewDataMatrix()
\sa SDPConeSetASparseVecMat()
\sa SDPConeSetADenseVecMat()
*/
int SDPConeAddADenseVecMat(SDPCone sdpcone,int blockj, int vari,int n,
			  double alpha,double val[], int nnz){
  int info;
  char UPLQ;
  void* dmat=0;
  struct DSDPDataMat_Ops* dmatops=0;
  
  DSDPFunctionBegin;
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  DSDPLogInfo(0,20,"Set dense matrix:  Block: %d, Variable %d, size: %d, Nonzeros: %d .\n",blockj,vari,n,nnz);
  switch (UPLQ){
  case 'P':
    info=DSDPGetDMat(n,alpha,val,&dmatops,&dmat); DSDPCHKERR(info);
    break;
  case 'U':
    DSDPSETERR(1,"Dense U Mat type does not exist.\n");
    break;
  }
  info=SDPConeAddDataMatrix(sdpcone,blockj,vari,n,UPLQ,dmatops,dmat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeAddDenseVecMat"
/* Needed for backward compatibility */
int SDPConeAddDenseVecMat(SDPCone sdpcone,int blockj, int vari,int n,
			  double val[], int nnz){
  int info;  
  DSDPFunctionBegin;
  info=SDPConeAddADenseVecMat(sdpcone,blockj,vari,n,1.0,val,nnz); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeSetADenseVecMat"
/*!
\fn int SDPConeSetADenseVecMat(SDPCone sdpcone,int blockj, int vari, int n, double alpha, double val[], int nnz);
\brief Set a matrix  \f$A_{i,j}\f$ in a dense format.
\ingroup SDPBasic
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param alpha multiple of the data (usually 1.0)
\param val array of elements in the matrix
\param nnz length of the array

For example, the matrix 
\image html img208.gif
can be inserted into the cone in packed symmetric format as follows
\code
double val[]={3,2,0,0,6,0};
SDPConeSetDenseAVecMat(sdpcone,j,i,3,1.0,val,6);
\endcode

\note DSDP will use the \c val array in many routines, but it will not modify or delete it.

\sa SDPConeViewDataMatrix()
\sa SDPConeSetASparseVecMat()

*/
int SDPConeSetADenseVecMat(SDPCone sdpcone,int blockj, int vari,int n,
			  double alpha,double val[], int nnz){
  int info;
  DSDPFunctionBegin;
  info=SDPConeRemoveDataMatrix(sdpcone,blockj,vari); DSDPCHKERR(info);
  info=SDPConeAddADenseVecMat(sdpcone,blockj,vari,n,alpha,val,nnz); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeSetDenseVecMat"
/* Needed for backward compatibility */
int SDPConeSetDenseVecMat(SDPCone sdpcone,int blockj, int vari,int n,
			  double val[], int nnz){
  int info;
  DSDPFunctionBegin;
  info=SDPConeSetADenseVecMat(sdpcone,blockj,vari,n,1.0,val,nnz); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeAddIdentity"
/*!
\fn int SDPConeAddIdentity(SDPCone sdpcone,int blockj, int vari, int n, double val);
\brief Add a matrix \f$A_{i,j}\f$ that is a multiple of the identity matrix.
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param val multiple of identity matrix
\sa SDPConeAddASparseVecMat()
*/
int SDPConeAddIdentity(SDPCone sdpcone, int blockj,int vari, int n,
		       double val){
  int info;
  char UPLQ;
  struct DSDPDataMat_Ops* identitymatops=0;
  void* imat=0;

  DSDPFunctionBegin;
  DSDPLogInfo(0,20,"Set identity matrix:  Block: %d, Variable %d, size: %d, Multiple: %4.4e .\n",blockj,vari,n,val);
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  switch (UPLQ){
  case 'P':
    info=DSDPGetIdentityDataMatP(n,val,&identitymatops,&imat);DSDPCHKERR(info);
    break;
  case 'U':
    info=DSDPGetIdentityDataMatF(n,val,&identitymatops,&imat);DSDPCHKERR(info);
    break;
  }
  info=SDPConeAddDataMatrix(sdpcone,blockj,vari,n,UPLQ,identitymatops,imat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeSetIdentity"
/*!
\fn int SDPConeSetIdentity(SDPCone sdpcone,int blockj, int vari, int n, double val);
\brief Set a matrix \f$A_{i,j}\f$ to be a multiple of the identity matrix.
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param val multiple of identity matrix
\sa SDPConeAddASparseVecMat()
*/
int SDPConeSetIdentity(SDPCone sdpcone, int blockj, int vari, int n,
		       double val){
  int info;
  DSDPFunctionBegin;
  info=SDPConeRemoveDataMatrix(sdpcone,blockj,vari); DSDPCHKERR(info);
  info=SDPConeAddIdentity(sdpcone,blockj,vari,n,val); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeAddConstantMat"
/*!
\fn int SDPConeAddConstantMat(SDPCone sdpcone,int blockj, int vari, int n, double value);
\brief Add a matrix \f$A_{i,j}\f$ whose elements are all the same.
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param value the value of each element in the matrix
\sa SDPConeAddADenseVecMat()
*/
int SDPConeAddConstantMat(SDPCone sdpcone,int blockj, int vari, int n,
			  double value){
  int info;
  char UPLQ;
  struct DSDPDataMat_Ops* constantmatops=0;
  void* smat=0;

  DSDPFunctionBegin;
  DSDPLogInfo(0,20,"Add allsame matrix:  Block: %d, Variable %d, size: %d, Elements: %4.4e .\n",blockj,vari,n,value);
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  switch (UPLQ){
  case 'P':
    info=DSDPGetConstantMat(n,value,UPLQ,&constantmatops,&smat);DSDPCHKERR(info);
    break;
  case 'U':
    info=DSDPGetConstantMat(n,value,UPLQ,&constantmatops,&smat);DSDPCHKERR(info);
    break;
  }
  info=SDPConeAddDataMatrix(sdpcone,blockj,vari,n,UPLQ,constantmatops,smat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeSetConstantMat"
/*!
\fn int SDPConeSetConstantMat(SDPCone sdpcone,int blockj, int vari, int n, double value);
\brief Set a matrix \f$A_{i,j}\f$ whose elements are all the same.
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param value the value of each element in the matrix
\sa SDPConeSetADenseVecMat()
*/
int SDPConeSetConstantMat(SDPCone sdpcone,int blockj, int vari, int n,
			  double value){
  int info;
  DSDPFunctionBegin;
  info=SDPConeRemoveDataMatrix(sdpcone,blockj,vari); DSDPCHKERR(info);
  info=SDPConeAddConstantMat(sdpcone,blockj,vari,n,value); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeSetZeroMat"
/*!
\fn int SDPConeSetZeroMat(SDPCone sdpcone,int blockj, int vari, int n);
\brief Set a matrix \f$A_{i,j}\f$ whose elements are all equal zero.
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\sa SDPConeSetASparseVecMat()
\sa SDPConeSetConstantMat()
*/
int SDPConeSetZeroMat(SDPCone sdpcone,int blockj, int vari, int n){
  int info;
  char UPLQ;
  struct DSDPDataMat_Ops* zeromatops=0;
  DSDPFunctionBegin;
  DSDPLogInfo(0,20,"Add zero matrix:  Block: %d, Variable %d, size: %d .\n",blockj,vari,n);
  info=DSDPGetZeroDataMatOps(&zeromatops); DSDPCHKERR(info);
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  info=SDPConeRemoveDataMatrix(sdpcone,blockj,vari); DSDPCHKERR(info);
  info=SDPConeAddDataMatrix(sdpcone,blockj,vari,n,UPLQ,zeromatops,0); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeSetRIdentity"
/*!
\fn int SDPConeSetRIdentity(SDPCone sdpcone,int blockj, int n, double rr);

\brief Add identify matrix to dual matrix.
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param n dimension of the matrix
\param rr multiple of identity matrix.
*/
int SDPConeSetRIdentity(SDPCone sdpcone,int blockj, int n, double rr){ 
  int info;
  char UPLQ;
  struct DSDPDataMat_Ops* identitymatops=0;
  void* imat=0;
  DSDPFunctionBegin;
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  switch (UPLQ){
  case 'P':
    info=DSDPGetIdentityDataMatP(n,rr,&identitymatops,&imat);DSDPCHKERR(info); break;
  case 'U':
    info=DSDPGetIdentityDataMatF(n,rr,&identitymatops,&imat);DSDPCHKERR(info); break;
  default:
    break;
  }
  info=SDPConeSetRMatrix(sdpcone,blockj,n,UPLQ,identitymatops,imat); DSDPCHKERR(info);

  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeAddARankOneMat"
/*!
\fn int SDPConeAddARankOneMat(SDPCone sdpcone,int blockj, int vari, int n, double alpha, int ishift, const int ind[], const double val[], int nnz);
\brief Add data matrix  \f$A_{i,j}= alpha * v * v^T \f$ where v is a sparse vector
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param alpha multiple of the outer product
\param ishift index of first element in a dense vector (usually 0)
\param val array of elements in the vector
\param ind array of indices representing the location of the nonzeros
\param nnz length of the previous two arrays
\note DSDP will use the two arrays throughout in many routines, but it will not modify or delete them.
\sa DSDPCreateSDPCone()
\sa SDPConeViewDataMatrix()
\sa SDPConeAddASparseVecMat()
*/
int SDPConeAddARankOneMat(SDPCone sdpcone,int blockj, int vari, int n,
			  double alpha, int ishift,const int ind[], const double val[], int nnz){
  
  int info;
  char UPLQ;
  void* dmat=0;
  struct DSDPDataMat_Ops* dmatops=0;
  
  DSDPFunctionBegin;
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  DSDPLogInfo(0,20,"Set sparse matrix:  Block: %d, Variable %d, size: %d, Nonzeros: %d .\n",blockj,vari,n,nnz);
  switch (UPLQ){
  case 'P':
    info=DSDPGetR1PMat(n,alpha,ishift,ind,val,nnz,&dmatops,&dmat); DSDPCHKERR(info);
    break;
  case 'U':
    info=DSDPGetR1UMat(n,alpha,ishift,ind,val,nnz,&dmatops,&dmat); DSDPCHKERR(info);
    break;
  }
  info=SDPConeAddDataMatrix(sdpcone,blockj,vari,n,UPLQ,dmatops,dmat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "SDPConeSetARankOneMat"
/*!
\fn int SDPConeSetARankOneMat(SDPCone sdpcone,int blockj, int vari, int n, double alpha, int ishift, const int ind[], const double val[], int nnz);
\brief Set data matrix  \f$A_{i,j}= alpha * v * v^T \f$ where v is a sparse vector
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param alpha multiple of the outer product
\param ishift index of first element in a dense vector (usually 0)
\param val array of elements in the vector
\param ind array of indices representing the location of the nonzeros
\param nnz length of the previous two arrays
\note DSDP will use the two arrays throughout in many routines, but it will not modify or delete them.
\sa SDPConeViewDataMatrix()
\sa SDPConeAddARankOneMat()
*/
int SDPConeSetARankOneMat(SDPCone sdpcone,int blockj, int vari, int n,
			  double alpha, int ishift,const int ind[], const double val[], int nnz){
  

  int info;
  DSDPFunctionBegin;
  info=SDPConeRemoveDataMatrix(sdpcone,blockj,vari); DSDPCHKERR(info);
  info=SDPConeAddARankOneMat(sdpcone,blockj,vari,n,alpha,ishift,ind,val,nnz); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSetDataMatZero"
/*!
\fn int DSDPSetDataMatZero(DSDPDataMat *A);
\brief Make a data matrix a zero matrix.
\param A data matrix.
*/
int DSDPSetDataMatZero(DSDPDataMat *A){
  int info;
  struct DSDPDataMat_Ops* zeromatops=0;
  DSDPFunctionBegin;
  info=DSDPGetZeroDataMatOps(&zeromatops); DSDPCHKERR(info);
  info=DSDPDataMatSetData(A,zeromatops,0);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

