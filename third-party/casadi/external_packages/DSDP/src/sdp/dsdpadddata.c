#include "dsdpsdp.h"
#include "dsdpsys.h"
/*!
\file dsdpadddata.c
\brief Set block sizes, sparsity, format, and data matrices
 */
#undef __FUNCT__  
#define __FUNCT__ "SDPConeCheckI"
/*!
\fn int SDPConeCheckI(SDPCone sdpcone, int vari);
\brief Check validity of parameter.
\param sdpcone SDP cone
\param vari variable i from 0 through m
*/
int SDPConeCheckI(SDPCone sdpcone,int vari){
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  if (vari<0 || vari>sdpcone->m) {
    DSDPSETERR2(1,"Bad Data Matrix: variable: %d (Max: %d)\n",vari,sdpcone->m+1);}
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeCheckJ"
/*!
\fn int SDPConeCheckJ(SDPCone sdpcone, int blockj);
\brief Check validity of parameter.
\param sdpcone SDP cone
\param blockj from 0 to nblocks
*/
int SDPConeCheckJ(SDPCone sdpcone,int blockj){
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  if (blockj<0 || blockj>= sdpcone->nblocks) {
    DSDPSETERR2(2,"Bad Data Matrix: Block: %d (Max: %d)\n",blockj,sdpcone->nblocks-1);}
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeCheckN"
/*!
\fn int SDPConeCheckN(SDPCone sdpcone, int blockj, int n);
\brief Check validity of parameter.
\param sdpcone SDP cone
\param blockj block number
\param n dimension of block.
*/
int SDPConeCheckN(SDPCone sdpcone,int blockj, int n){
  int info;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  if (sdpcone->blk[blockj].n==0 && n>0){info=SDPConeSetBlockSize(sdpcone,blockj,n);DSDPCHKERR(info);}
  if (sdpcone->blk[blockj].n != n){
    DSDPSETERR3(3,"Check Dimension of Data Matrix: Block: %d, %d -- expecting %d\n",
		blockj,n,sdpcone->blk[blockj].n);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeCheckM"
/*!
\fn int SDPConeCheckM(SDPCone sdpcone, int m);
\brief Check validity of parameter.
\param sdpcone SDP cone
\param m number of y variables
*/
int SDPConeCheckM(SDPCone sdpcone,int m){
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  if (m!=sdpcone->m){
    DSDPSETERR1(4,"Check dimension of array. This problem has %d variables\n",sdpcone->m);}
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeValidStorageFormat"
/*!
\fn int SDPConeValidStorageFormat(SDPCone sdpcone, char format);
\brief Check validity of parameter.
\param sdpcone SDP cone
\param format such as packed symmetric or upper full symmetric
*/
int SDPConeValidStorageFormat(SDPCone sdpcone, char format){
  DSDPFunctionBegin;
  if (format!='P' && format != 'U'){
    DSDPSETERR1(4,"Check format of Block: %c is not supported! Use P or U. \n",format);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeCheckStorageFormat"
/*!
\fn int SDPConeCheckStorageFormat(SDPCone sdpcone,int blockj, char format);
\brief Check validity of parameters.
\param sdpcone SDP cone
\param blockj block number
\param format such as packed symmetric or upper full symmetric
*/
int SDPConeCheckStorageFormat(SDPCone sdpcone,int blockj, char format){
  int info;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=SDPConeValidStorageFormat(sdpcone,format);DSDPCHKERR(info);
  if (sdpcone->blk[blockj].format=='N'){
    sdpcone->blk[blockj].format = format;
  }
  if (sdpcone->blk[blockj].format != format){
    DSDPSETERR3(4,"Check format of Data Matrix: Block: %d, %c -- expecting %c\n",
		blockj,format,sdpcone->blk[blockj].format);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeRemoveDataMatrix"
/*!
\fn int SDPConeRemoveDataMatrix(SDPCone sdpcone,int blockj, int vari);
\brief Remove the data matrix \f$A_{i,j}\f$ from the cone
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\sa SDPConeView()
*/
int SDPConeRemoveDataMatrix(SDPCone sdpcone,int blockj, int vari){
  int info;
  DSDPFunctionBegin;
  info=SDPConeCheckI(sdpcone,vari);DSDPCHKERR(info);
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=DSDPBlockRemoveDataMatrix(&sdpcone->blk[blockj].ADATA,vari);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
  
#undef __FUNCT__  
#define __FUNCT__ "SDPConeAddDataMatrix"
/*!
\fn int SDPConeAddDataMatrix(SDPCone sdpcone,int blockj, int vari, int n, char format,struct DSDPDataMat_Ops* dsdpdataops, void* data);
\brief Add a data matrix \f$A_{i,j}\f$
\ingroup SDPData
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param vari variable i from 0 through m
\param n dimension of the matrix
\param format storage format 'P' (default) or 'U'
\param data address of a structure ( cast to \c void* ) with matrix data.
\param dsdpdataops address of a structure of function pointers that operate on the matrix data
\sa SDPConeSetASparseVecMat()
\sa SDPConeSetStorageFormat()
\sa SDPConeSetBlockSize()
\sa SDPConeCheckData()
*/
int SDPConeAddDataMatrix(SDPCone sdpcone,int blockj, int vari, int n, char format, struct DSDPDataMat_Ops* dsdpdataops, void* data){ 
  int info;
  DSDPFunctionBegin;
  info=SDPConeCheckI(sdpcone,vari);DSDPCHKERR(info);
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=SDPConeCheckN(sdpcone,blockj,n);DSDPCHKERR(info);
  info=SDPConeCheckStorageFormat(sdpcone,blockj,format);DSDPCHKERR(info);
  info=DSDPBlockAddDataMatrix(&sdpcone->blk[blockj].ADATA,vari,dsdpdataops,data);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeSetRMatrix"
/*!
\fn int SDPConeSetRMatrix(SDPCone sdpcone,int blockj, int n, char format,struct DSDPDataMat_Ops* dsdpdataops, void* data);
\brief Add identity to dual matrix.
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param n dimension of the matrix
\param format storage format 'P' (default) or 'U'
\param data address of a structure ( cast to \c void* ) with matrix data.
\param dsdpdataops address of a structure of function pointers that operate on the matrix data
\sa SDPConeSetASparseVecMat()
\sa SDPConeSetStorageFormat()
\sa SDPConeSetBlockSize()
\sa SDPConeCheckData()
*/
int SDPConeSetRMatrix(SDPCone sdpcone,int blockj, int n, char format, struct DSDPDataMat_Ops* dsdpdataops, void* data){ 
  int info;
  int vari=sdpcone->m+1;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=SDPConeCheckN(sdpcone,blockj,n);DSDPCHKERR(info);
  info=SDPConeCheckStorageFormat(sdpcone,blockj,format);DSDPCHKERR(info);
  info=DSDPBlockRemoveDataMatrix(&sdpcone->blk[blockj].ADATA,vari);DSDPCHKERR(info);
  info=DSDPBlockSetDataMatrix(&sdpcone->blk[blockj].ADATA,vari,dsdpdataops,data);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeViewDataMatrix"
/*!
\fn int SDPConeViewDataMatrix(SDPCone sdpcone,int blockj, int vari);
\ingroup SDPBasic
\brief Print a data matrix to the screen
\param sdpcone semidefinite cone object
\param blockj block number
\param vari variable number from 0 through m
\sa SDPConeView()
*/
int SDPConeViewDataMatrix(SDPCone sdpcone,int blockj, int vari){ 
  int info,ii,vari2,nnzmats;
  DSDPDataMat AA;
  DSDPFunctionBegin;
  info=SDPConeCheckI(sdpcone,vari);DSDPCHKERR(info);
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=DSDPBlockCountNonzeroMatrices(&sdpcone->blk[blockj].ADATA,&nnzmats);DSDPCHKERR(info);
  for (ii=0;ii<nnzmats; ii++){  /* Matrix Entries */
    info=DSDPBlockGetMatrix(&sdpcone->blk[blockj].ADATA,ii,&vari2,0,&AA);DSDPCHKVARERR(vari,info);
    if (vari2==vari){ info = DSDPDataMatView(AA);DSDPCHKERR(info);}
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeClearVMatrix"
/*!
\fn int SDPConeClearVMatrix(SDPCone sdpcone,int blockj, int n);
\brief Free V matrix.
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param n dimension of the matrix
*/
int SDPConeClearVMatrix(SDPCone sdpcone,int blockj, int n){
  int info;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=DSDPVMatDestroy(&sdpcone->blk[blockj].T);DSDPCHKERR(info);
  info=DSDPVMatInitialize(&sdpcone->blk[blockj].T);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "SDPConeSetXMat"
/*!
\fn int SDPConeSetXMat(SDPCone sdpcone,int blockj, int n);
\brief Create X matrix.
\param sdpcone SDP cone
\param blockj block number j from 0 to nblocks
\param n dimension of the matrix
*/
int SDPConeSetXMat(SDPCone sdpcone,int blockj, int n){ 
  int info;
  char UPLQ;
  DSDPVMat T;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=SDPConeClearVMatrix(sdpcone,blockj,n);DSDPCHKERR(info);
  DSDPLogInfo(0,10,"Create block X Mat:  Block: %d, size: %d.\n",blockj,n);
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  info=DSDPMakeVMat(UPLQ,n,&T);DSDPCHKERR(info);
  sdpcone->blk[blockj].T=T;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeSetXArray"
/*! 
\fn int SDPConeSetXArray(SDPCone sdpcone,int blockj, int n, double xx[], int nn);
\ingroup SDPRoutines
\brief Provide an array for the SDPCone object can use to store dense matrices.
\param sdpcone semidefinite cone object
\param blockj block number
\param n dimension of the block
\param xx array for dense matrices
\param nn length of array

This routine elimates the need to copy the solution X into a separate array.

\sa SDPConeGetXArray()
 */
int SDPConeSetXArray(SDPCone sdpcone,int blockj, int n, double xx[], int nn){ 
  int info;
  char UPLQ;
  DSDPVMat T;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=SDPConeCheckN(sdpcone,blockj,n);DSDPCHKERR(info);
  info=SDPConeClearVMatrix(sdpcone,blockj,n);DSDPCHKERR(info);
  DSDPLogInfo(0,10,"Set block X array:  Block: %d, size: %d.\n",blockj,n);
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ); DSDPCHKERR(info);
  info=DSDPMakeVMatWithArray(UPLQ,xx,nn,n,&T);DSDPCHKERR(info);
  sdpcone->blk[blockj].T=T;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeGetXArray"
/*! 
\fn int SDPConeGetXArray(SDPCone sdpcone,int blockj, double* xx[], int *nn);
\ingroup SDPBasic
\brief After applying the solver, set a pointer to the array in the object 
with the solution X.
\param sdpcone semidefinite cone object
\param blockj block number
\param xx address of an array for dense matrices
\param *nn the length of the array
\sa DSDPSolve()
\sa DSDPComputeX()
\sa SDPConeViewX()

\code
DSDP dsdp;
SDPCone sdpcone;
double *xx;
int nn;

DSDPSolve(dsdp);
DSDPComputeX(dsdp);
SDPConeGetXArray(sdpcone,0,&xx,&nn);
SDPConeViewX(sdpcone,0,xx,nn);
SDPConeRestoreXArray(sdpcone,0,&xx,&nn);
\endcode

DSDP uses a single dense array to add data matrices, 
compute the matrix X, and take the inner
product of X with the data matrices.  
Therefore, the ordering of elements in this array 
must also be used in the data matrices.

 */
int SDPConeGetXArray(SDPCone sdpcone,int blockj, double* xx[], int *nn){ 
  int info,flag;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=DSDPVMatExist(sdpcone->blk[blockj].T,&flag);DSDPCHKERR(info);
  if (flag==0){
    DSDPSETERR(6,"No X Array available, Call DSDPSetup() or SDPConeSetXArray.\n");}
  info=DSDPVMatGetArray(sdpcone->blk[blockj].T,xx,nn);DSDPCHKERR(info); 
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeRestoreXArray"
/*! 
\fn int SDPConeRestoreXArray(SDPCone sdpcone,int blockj, double* xx[], int *nn);
\ingroup SDPRoutines
\brief Restore the dense array and set these pointers to null.
\param sdpcone semidefinite cone object
\param blockj block number
\param xx address of an array for dense matrices
\param *nn the length of the array
\sa SDPConeGetXArray()
 */
int SDPConeRestoreXArray(SDPCone sdpcone,int blockj, double* xx[], int *nn){ 
  int info,flag;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=DSDPVMatExist(sdpcone->blk[blockj].T,&flag);DSDPCHKERR(info);
  if (flag==0){
    DSDPSETERR(6,"No X Array available, Call DSDPSetup() or SDPConeSetXArray.\n");}
  info=DSDPVMatRestoreArray(sdpcone->blk[blockj].T,xx,nn);DSDPCHKERR(info); 
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeMatrixView"
/*! 
\fn int SDPConeMatrixView(SDPCone sdpcone,int blockj);
\ingroup SDPRoutines
\brief Print the dense array to the screen.
\param sdpcone semidefinite cone object
\param blockj block number
\sa SDPConeGetXArray()
 */
int SDPConeMatrixView(SDPCone sdpcone, int blockj){
  int info;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=DSDPVMatView(sdpcone->blk[blockj].T);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeUseFullSymmetricFormat"
/*! 
\fn int SDPConeUseFullSymmetricFormat(SDPCone sdpcone,int blockj);
\ingroup SDPData
\brief Use full symmetric format for the dense array.
\param sdpcone semidefinite cone object
\param blockj block number

\sa SDPConeGetXArray()
\sa SDPConeSetASparseVecMat()
\sa SDPConeSetADenseVecMat()
\sa SDPConeSetStorageFormat()
\sa SDPConeUsePackedFormat()
\sa SDPConeCheckData()

In this format, an \f$n \times n\f$ symmetric matrix is represented by
an array of length \f$ n^2 \f$.  Each element of the array
corresponds to an element of the matrix and the mapping is given below.
\f[
\begin{array}{lllllllllll}
[ a_{1,1} & \ldots & a_{1,n} & a_{2,1} & \ldots & a_{2,n} & \ldots & a_{n,1} & \ldots & a_{n,n} ] \\  
\end{array}
\f]
but elements \f$a_{i,j}, i<j \f$ are not used.  

DSDP uses a single dense array to add data matrices, compute the matrix X, and take the inner
product of X with the data matrices.  Therefore, the ordering of elements in this array 
must also be used in the data matrices.

\note Since DSDP uses BLAS and LAPACK for many of its operations,
we refer to the half of matrix in use as the 'upper' triangular matrix.

*/
int SDPConeUseFullSymmetricFormat(SDPCone sdpcone, int blockj){
  int info;
  DSDPFunctionBegin;
  info=SDPConeSetStorageFormat(sdpcone,blockj,'U');DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeUsePackedFormat"
/*! 
\fn int SDPConeUsePackedFormat(SDPCone sdpcone,int blockj);
\ingroup SDPData
\brief Use packed symmetric format for the dense array.
\param sdpcone semidefinite cone object
\param blockj block number

\sa SDPConeGetXArray()
\sa SDPConeSetASparseVecMat()
\sa SDPConeSetADenseVecMat()

In this format, an \f$n \times n\f$ symmetric matrix is represented by
an array of length \f$\frac{n(n+1)}{2}\f$.  Each element of the array
corresponds to an element of the matrix and the mapping is given below.
\f[
\begin{array}{llllllll}
[ a_{1,1} & a_{2,1} & a_{2,2} & a_{3,1} & a_{3,2} & a_{3,3} & \ldots & a_{n,n} ] \\  
\end{array}
\f]

DSDP uses a single dense array to add data matrices, 
compute the matrix X, and take the inner
product of X with the data matrices.  
Therefore, the ordering of elements in this array 
must also be used in the data matrices.

\note Since DSDP uses BLAS and LAPACK for many of its operations,
we refer to the half of matrix in use as the 'upper' triangular matrix.
*/
int SDPConeUsePackedFormat(SDPCone sdpcone, int blockj){
  int info;
  DSDPFunctionBegin;
  info=SDPConeSetStorageFormat(sdpcone,blockj,'P');DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeSetStorageFormat"
/*! 
\fn int SDPConeSetStorageFormat(SDPCone sdpcone,int blockj, char format);
\ingroup SDPData
\brief Set the dense storage format of a block in the semidefinite cone.
\param sdpcone semidefinite cone object
\param blockj block number
\param format format the block
\sa SDPConeGetStorageFormat()
\sa SDPConeUsePackedFormat()
\sa SDPConeUseFullSymmetricFormat()

This routine determines the ordering of elements in the data matrices
and the X matrix.
The default format is 'P' (packed symmetric format), but the full symmetric
format 'U' is also supported.  The format used to factor the S
matrix is independent of the format used for the data and X matrix.

*/
int SDPConeSetStorageFormat(SDPCone sdpcone, int blockj, char format){
  int info;
  DSDPFunctionBegin;
  info=SDPConeValidStorageFormat(sdpcone,format);DSDPCHKERR(info);
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  sdpcone->blk[blockj].format=format;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeGetStorageFormat"
/*! 
\fn int SDPConeGetStorageFormat(SDPCone sdpcone,int blockj, char *format);
\ingroup SDPData
\brief Get the storage format for the block.
\param sdpcone semidefinite cone object
\param blockj block number
\param format format the block
\sa SDPConeSetStorageFormat()
\sa SDPConeGetXArray()
\sa SDPConeSetASparseVecMat()
\sa SDPConeSetADenseVecMat()

The default format is 'P' (packed symmetric format).

*/
int SDPConeGetStorageFormat(SDPCone sdpcone, int blockj, char *format){
  int info;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  *format=sdpcone->blk[blockj].format;
  if (*format=='N') *format='P';
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeScaleBarrier"
int SDPConeScaleBarrier(SDPCone sdpcone,int blockj, double ggamma){ 
  int info;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  sdpcone->blk[blockj].gammamu=ggamma;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeSetBlockSize"
/*! 
\fn int SDPConeSetBlockSize(SDPCone sdpcone,int blockj, int n);
\ingroup SDPRoutines
\brief Set the dimension of one block in the semidefinite cone
\param sdpcone semidefinite cone object
\param blockj block number
\param n dimension of the block
\note This routine is optional, but its use can help the object detect input errors.
 */
int SDPConeSetBlockSize(SDPCone sdpcone, int blockj, int n){
  int info,n0;
  DSDPFunctionBegin;
  DSDPLogInfo(0,10,"Set block size:  Block: %d, size: %d.\n",blockj,n);
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  n0=sdpcone->blk[blockj].n;
  if (n0==n){DSDPFunctionReturn(0);}
  if (n0!=0 &&n0!=n){
    DSDPSETERR2(5,"Block %d Size previously set to %d \n",blockj,n0);  } 
  sdpcone->blk[blockj].n=n;
  sdpcone->nn+=n-n0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeGetBlockSize"
/*! 
\fn int SDPConeGetBlockSize(SDPCone sdpcone,int blockj, int *n);
\ingroup SDPRoutines
\brief Get the dimension of one block in the semidefinite cone
\param sdpcone semidefinite cone object
\param blockj block number
\param *n set to the dimension of the block
\sa SDPConeSetBlockSize()
 */
int SDPConeGetBlockSize(SDPCone sdpcone, int blockj, int *n){
  int info;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  *n=sdpcone->blk[blockj].n;
  DSDPFunctionReturn(0);
}

/*! 
\fn int SDPConeGetNumberOfBlocks(SDPCone sdpcone, int *nblocks);
\ingroup SDPRoutines
\brief Get the number of blocks in the semidefinite cone
\param sdpcone semidefinite cone object
\param *nblocks set to the dimension of the block
\sa DSDPCreateSDPCone()
 */
#undef __FUNCT__  
#define __FUNCT__ "SDPConeGetNumberOfBlocks"
int SDPConeGetNumberOfBlocks(SDPCone sdpcone, int *nblocks){
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  *nblocks=sdpcone->nblocks;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeSetSparsity"
/*! 
\fn int SDPConeSetSparsity(SDPCone sdpcone,int blockj, int nnz);
\ingroup SDPRoutines
\brief Set the number of nonzero matrices in a block of the semidefinite cone.
\param sdpcone semidefinite cone object
\param blockj block number
\param nnz number of nonzero matrices in the block
\note This routine is optional, but its use can improve the memory allocation within the SDPCone object.
 */
int SDPConeSetSparsity(SDPCone sdpcone, int blockj, int nnz){
  int info;
  DSDPFunctionBegin;
  DSDPLogInfo(0,10,"Set block nonzeros:  Block: %d, Nonzero Matrices: %d.\n",blockj,nnz);
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  if (nnz>sdpcone->m) nnz=sdpcone->m;
  info=DSDPBlockDataAllocate(&sdpcone->blk[blockj].ADATA,nnz+2); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "SDPConeView"
/*! \fn int SDPConeView(SDPCone sdpcone);
\brief Print the SDPCone to the screen;
\param sdpcone the cone
\sa SDPConeViewDataMatrix()
\sa SDPConeView2()
\sa SDPConeView2()
\ingroup SDPBasics
*/
int SDPConeView(SDPCone sdpcone){  
  int blockj,info;
  DSDPFunctionBegin;
  for (blockj=0; blockj<sdpcone->nblocks; blockj++){
    printf("Block: %d, Dimension: %d\n",blockj,sdpcone->blk[blockj].n);
    info=DSDPBlockView(&sdpcone->blk[blockj].ADATA);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeView2"
/*! \fn int SDPConeView2(SDPCone sdpcone);
\brief Print the SDP cone to the screen in a second way
\param sdpcone the cone
\sa SDPConeViewDataMatrix()
\sa SDPConeView()
\sa SDPConeView3()
\ingroup SDPData
*/
int SDPConeView2(SDPCone sdpcone){  
  int blockj,info;
  DSDPFunctionBegin;
  for (blockj=0; blockj<sdpcone->nblocks; blockj++){
    printf("Block: %d, Dimension: %d\n",blockj,sdpcone->blk[blockj].n);
    info=DSDPBlockView2(&sdpcone->blk[blockj].ADATA);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeView3"
/*! \fn int SDPConeView3(SDPCone sdpcone);
\brief Print the SDP cone to the screen in a third way
\param sdpcone the cone
\sa SDPConeViewDataMatrix()
\sa SDPConeView()
\sa SDPConeView2()
\ingroup SDPData
*/
int SDPConeView3(SDPCone sdpcone){  
  int blockj,id,n,info,nnzmats;
  DSDPFunctionBegin;
  for (blockj=0; blockj<sdpcone->nblocks; blockj++){
    n=sdpcone->blk[blockj].n;
    printf("Block: %d \n",blockj);
    printf(" Dimension: %d\n",n);
    info=DSDPDSMatGetType(sdpcone->blk[blockj].DS,&id);
    if (id==1){
      printf(" DS Matrix Type: Dense, Using LAPACK\n");
    } else {
      printf(" DS Matrix Type: %d\n",id);
    }
    info=DSDPDualMatGetType(sdpcone->blk[blockj].S,&id);
    if (id==1){
      printf(" Dual Matrix Type: Dense, Using LAPACK\n");
    } else {
      printf(" Dual Matrix Type: %d\n",id);
    }
    info=DSDPBlockCountNonzeroMatrices(&sdpcone->blk[blockj].ADATA,&nnzmats);DSDPCHKERR(info);
    printf(" Number of Data Matrices: %d of %d\n",nnzmats-1,sdpcone->m+1);
    printf(" Number of Data Nonzeros: %d\n",sdpcone->blk[blockj].nnz);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "SDPConeCheckData"
/*! \fn int SDPConeCheckData(SDPCone sdpcone);
\brief Check the matrix operations on a data matrix;
\param sdpcone the cone
\ingroup SDPData
\sa SDPConeView
*/
int SDPConeCheckData(SDPCone sdpcone){  
  int i,ii,blockj,nnzmats,info;
  double scl=0;
  DSDPDataMat AA;
  DSDPIndex IS;
  DSDPVMat T;
  DSDPDSMat DS;
  DSDPDualMat S1,S2;
  SDPConeVec W,W2;
  DSDPFunctionBegin;
  for (blockj=0; blockj<sdpcone->nblocks; blockj++){
    T=sdpcone->blk[blockj].T;DS=sdpcone->blk[blockj].DS;
    W=sdpcone->blk[blockj].W;W2=sdpcone->blk[blockj].W2;
    S1=sdpcone->blk[blockj].S;S2=sdpcone->blk[blockj].SS;
    IS=sdpcone->blk[blockj].IS;
    printf("Block: %d\n",blockj);
    info=DSDPVMatCheck(T,W,W2);DSDPCHKERR(info);
    info=DSDPDSMatCheck(DS,W,W2,T);DSDPCHKERR(info);
    info=DSDPDualMatCheck(S1,W,W2,IS,T);DSDPCHKERR(info);
    info=DSDPDualMatCheck(S2,W,W2,IS,T);DSDPCHKERR(info);

    info=DSDPBlockCountNonzeroMatrices(&sdpcone->blk[blockj].ADATA,&nnzmats);DSDPCHKERR(info);
    for (ii=0;ii<nnzmats;ii++){
      info=DSDPBlockGetMatrix(&sdpcone->blk[blockj].ADATA,ii,&i,&scl,&AA);DSDPCHKERR(info);
      if (i==0) continue;
      printf(" Variable: %d, \n",i);
      info=DSDPDataMatCheck(AA,W,IS,T);DSDPCHKERR(info);
    }
  }
  DSDPFunctionReturn(0);
}
