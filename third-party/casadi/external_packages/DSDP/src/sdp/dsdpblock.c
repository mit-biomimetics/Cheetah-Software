#include "dsdpsdp.h"
#include "dsdpsys.h"

static int sdpvecvecevent=0,sdpdotevent=0;
/*!
\file dsdpblock.c
\brief Operations on a single SDP block.
*/

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockASum"
/*!
\fn int DSDPBlockASum(DSDPBlockData *ADATA, double aa, DSDPVec Yk, DSDPVMat XX){
\brief Sum the data matrices.
\param ADATA block of data.
\param aa scalar
\param Yk scalar.
\param XX equals aa * sum ( YK[i]* A[i] )
*/
int DSDPBlockASum(DSDPBlockData *ADATA, double aa, DSDPVec Yk, DSDPVMat XX){

  double *xx,ytmp,scl=ADATA->scl;
  int    ii,vari,n,nn,info;

  DSDPFunctionBegin;
  info=DSDPVMatGetSize(XX, &n); DSDPCHKERR(info);
  info=DSDPVMatGetArray(XX, &xx, &nn); DSDPCHKERR(info);
  for (ii=0;ii<ADATA->nnzmats;ii++){
    vari=ADATA->nzmat[ii];
    info=DSDPVecGetElement(Yk,vari,&ytmp);DSDPCHKVARERR(vari,info);
    if (ytmp==0) continue;
    info = DSDPDataMatAddMultiple(ADATA->A[ii], -aa*scl*ytmp, xx,nn,n); DSDPCHKVARERR(vari,info);
  }
  info=DSDPVMatRestoreArray(XX, &xx, &nn); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockADot"
/*!
\fn int DSDPBlockADot(DSDPBlockData *ADATA, double aa, DSDPVec Alpha, DSDPVMat X, DSDPVec AX){
\brief Compute inner product of XX with data matrices.
\param ADATA block of data.
\param aa scalar
\param Alpha scalar.
\param X Dense symmetric matrix
\param AX Set AX[i] to aa * Alpha[i] * Dot( A[i] * X)
*/
int DSDPBlockADot(DSDPBlockData *ADATA, double aa, DSDPVec Alpha, DSDPVMat X, DSDPVec AX){

  int    ii,vari,n,nn,info;
  double *x,sum=0,aalpha=0,scl=ADATA->scl;

  DSDPFunctionBegin;
  DSDPEventLogBegin(sdpdotevent);
  info=DSDPVMatScaleDiagonal(X,0.5); DSDPCHKERR(info);
  info=DSDPVMatGetSize(X, &n); DSDPCHKERR(info);
  info=DSDPVMatGetArray(X, &x, &nn); DSDPCHKERR(info);
  for (ii=0;ii<ADATA->nnzmats; ii++){  /* Matrix Entries */
    vari=ADATA->nzmat[ii];
    info=DSDPVecGetElement(Alpha,vari,&aalpha);DSDPCHKVARERR(vari,info);
    if (aalpha==0.0) continue;
    info=DSDPDataMatDot(ADATA->A[ii],x,nn,n,&sum);DSDPCHKVARERR(vari,info);
    info=DSDPVecAddElement(AX,vari,aa*aalpha*sum*scl);DSDPCHKVARERR(vari,info);
  }
  info=DSDPVMatRestoreArray(X, &x, &nn); DSDPCHKERR(info);
  info=DSDPVMatScaleDiagonal(X,2.0); DSDPCHKERR(info);
  DSDPEventLogEnd(sdpdotevent);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockvAv"
/*!
\fn int int DSDPBlockvAv(DSDPBlockData *ADATA, double aa, DSDPVec Alpha, SDPConeVec V, DSDPVec VAV){

\brief Set VAV[i] to aa * Alpha[i] * V' A[i] V
\param ADATA block of data.
\param aa scalar
\param Alpha scalar.
\param V vecotr
\param VAV the product
*/
int DSDPBlockvAv(DSDPBlockData *ADATA, double aa, DSDPVec Alpha, SDPConeVec V, DSDPVec VAV){

  int    ii,vari,info;
  double sum=0,aalpha=0,scl=ADATA->scl;

  DSDPFunctionBegin;
  DSDPEventLogBegin(sdpvecvecevent);
  if (aa==0){DSDPFunctionReturn(0);}
  for (ii=0;ii<ADATA->nnzmats; ii++){  /* Matrix Entries */
    vari=ADATA->nzmat[ii];
    info=DSDPVecGetElement(Alpha,vari,&aalpha);DSDPCHKVARERR(vari,info);
    if (aalpha==0.0) continue;
    info=DSDPDataMatVecVec(ADATA->A[ii],V,&sum);DSDPCHKVARERR(vari,info);
    info=DSDPVecAddElement(VAV,vari,aa*aalpha*sum*scl);DSDPCHKVARERR(vari,info);
  }
  DSDPEventLogEnd(sdpvecvecevent);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockFactorData"
/*!
\fn int DSDPBlockFactorData(DSDPBlockData *ADATA, DSDPVMat X, SDPConeVec W){

\brief Factor the data matrices
\param ADATA block of data.
\param X work array
\param W Work vector
*/
int DSDPBlockFactorData(DSDPBlockData *ADATA, 
			DSDPVMat X, SDPConeVec W){
  
  int    ii,vari,n,nn,info,*iwork3n,i13,n26;
  double *x,*dwork3n;
  DSDPDataMat AA;

  DSDPFunctionBegin;
  info=DSDPVMatGetSize(X, &n); DSDPCHKERR(info);
  i13=13*n+1;n26=26*n+1;
  DSDPCALLOC2(&dwork3n,double,n26,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&iwork3n,int,i13,&info);DSDPCHKERR(info);
  info=DSDPVMatGetArray(X, &x, &nn); DSDPCHKERR(info);
  for (ii=0;ii<ADATA->nnzmats; ii++){  /* Matrix Entries */
    info=DSDPBlockGetMatrix(ADATA,ii,&vari,0,&AA);DSDPCHKVARERR(vari,info);
    DSDPLogInfo(0,39,"SDP Data Mat Setup: %d\n",vari);
    if (vari==0) continue;
    info=DSDPDataMatFactor(AA,W,x,nn,dwork3n,n26,iwork3n,i13); DSDPCHKVARERR(vari,info);
  }
  info=DSDPVMatRestoreArray(X, &x, &nn); DSDPCHKERR(info);
  DSDPFREE(&dwork3n,&info);DSDPCHKERR(info);
  DSDPFREE(&iwork3n,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPBlockEventZero"
int DSDPBlockEventZero(void){
  DSDPFunctionBegin;
  sdpvecvecevent=0;sdpdotevent=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockEventInitialize"
int DSDPBlockEventInitialize(void){
  DSDPFunctionBegin;
  if (sdpvecvecevent==0){DSDPEventLogRegister("SDP VecMatVec",&sdpvecvecevent);}
  if (sdpdotevent==0){DSDPEventLogRegister("SDP Dot",&sdpdotevent);}
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockDataInitialize"
/*!
\fn int DSDPBlockDataInitialize(DSDPBlockData *ADATA){
\brief Set pointers to null.
\param ADATA block of data.
*/
int DSDPBlockDataInitialize(DSDPBlockData *ADATA){
  DSDPFunctionBegin;
  ADATA->nnzmats=0;
  ADATA->maxnnzmats=0;
  ADATA->nzmat=0;
  ADATA->A=0;
  ADATA->r=1.0;
  ADATA->scl=1.0;
  /*  ADATA->n=0; */
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockTakeDownData"
/*!
\fn int DSDPBlockTakeDownData(DSDPBlockData *ADATA){
\brief Free structures in block of data.
\param ADATA block of data.
*/
int DSDPBlockTakeDownData(DSDPBlockData *ADATA){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPBlockDataDestroy"
/*!
\fn int DSDPBlockDataDestroy(DSDPBlockData *ADATA);
\brief Free the data matrices
\param ADATA block of data.
*/
int DSDPBlockDataDestroy(DSDPBlockData *ADATA){
  int ii,vari,info;
  DSDPFunctionBegin;
  if (!ADATA){DSDPFunctionReturn(0);}
  DSDPLogInfo(0,18,"Destroying All Existing Data Matrices \n");
  for (ii=0; ii<ADATA->nnzmats; ii++){
    vari=ADATA->nzmat[ii];
    info = DSDPDataMatDestroy(&ADATA->A[ii]);DSDPCHKVARERR(vari,info);
    ADATA->nzmat[ii]=0;
  }
  ADATA->nnzmats=0;
  info=DSDPBlockTakeDownData(ADATA);DSDPCHKERR(info);
  DSDPFREE(&ADATA->nzmat,&info);DSDPCHKERR(info);
  DSDPFREE(&ADATA->A,&info);DSDPCHKERR(info);
  info=DSDPBlockDataInitialize(ADATA);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPBlockDataAllocate"
/*!
\fn int DSDPBlockDataAllocate(DSDPBlockData *ADATA, int nnz);
\brief Allocate some structures.
\param ADATA block of data.
\param nnz number of data matrices to allocate space
*/
int DSDPBlockDataAllocate(DSDPBlockData *ADATA, int nnz){
  int j,info,*nzmat;
  DSDPDataMat *A;
  DSDPFunctionBegin;
  if (!ADATA){DSDPFunctionReturn(0);}
  if (nnz<=ADATA->maxnnzmats){DSDPFunctionReturn(0);}
  DSDPLogInfo(0,18,"REALLOCATING SPACE FOR %d SDP BLOCK MATRICES! Previously allocated: %d \n",nnz,ADATA->maxnnzmats);
  DSDPCALLOC2(&A,struct DSDPDataMat_C,nnz,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&nzmat,int,nnz,&info);DSDPCHKERR(info);
  for (j=0;j<nnz;j++){nzmat[j]=0;}
  for (j=0;j<nnz;j++){info = DSDPDataMatInitialize(&A[j]);DSDPCHKERR(info);}
  if (ADATA->maxnnzmats>0){
    for (j=0;j<ADATA->nnzmats;j++){nzmat[j]=ADATA->nzmat[j];}
    for (j=0;j<ADATA->nnzmats;j++){A[j]=ADATA->A[j];}
    DSDPFREE(&ADATA->A,&info);DSDPCHKERR(info);
    DSDPFREE(&ADATA->nzmat,&info);DSDPCHKERR(info);
  } else {
    ADATA->nnzmats=0;
  }
  ADATA->maxnnzmats=nnz;
  ADATA->nzmat=nzmat;
  ADATA->A=A;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockDataMarkNonzeroMatrices"
/*!
\fn int DSDPBlockDataMarkNonzeroMatrices(DSDPBlockData *ADATA, int *annz){
\brief Mark which variable in block have a data matrix.
\param ADATA block of data.
\param annz array of integers to mark.
*/
int DSDPBlockDataMarkNonzeroMatrices(DSDPBlockData *ADATA,int*annz){
  int i;
  DSDPFunctionBegin;
  for (i=0; i<ADATA->nnzmats; i++){
    annz[ADATA->nzmat[i]]++;
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockCountNonzerosMatrices"
/*!
\fn int DSDPBlockCountNonzeroMatrices(DSDPBlockData *ADATA,int*nzmats);

\brief Count how many data matrices are in a block of data.
\param ADATA block of data.
\param nzmats set to number of nonzero matrices.  Numbers from 0 to this number can be used as a matrix id in DSDPBlockGetMatrix()
*/
int DSDPBlockCountNonzeroMatrices(DSDPBlockData *ADATA,int*nzmats){
  DSDPFunctionBegin;
  *nzmats=ADATA->nnzmats;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockDataRank"
int DSDPBlockDataRank(DSDPBlockData *ADATA, int *trank, int n){
  int    ii,vari,info,ri,r2=0;
  DSDPDataMat AA;

  DSDPFunctionBegin;
  for (ii=0;ii<ADATA->nnzmats;ii++){
    info=DSDPBlockGetMatrix(ADATA,ii,&vari,0,&AA);DSDPCHKVARERR(vari,info);
    if (vari==0) continue;
    info=DSDPDataMatGetRank(AA,&ri,n); DSDPCHKVARERR(vari,info);
    r2+=ri;
  }
  *trank=r2;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockGetMatrix"
/*!
\fn int DSDPBlockGetMatrix(DSDPBlockData *ADATA,int id, int *vari, double *scl, DSDPDataMat *A);

\brief Get a data matrix from a block of data.
\param ADATA block of data.
\param id identfier of the matrices, numbered consecutively from 0.
\param vari set to variable number corresponding to A.
\param scl scaling
\param A data matrix.
*/
int DSDPBlockGetMatrix(DSDPBlockData *ADATA,int id, int *vari, double *scl, DSDPDataMat *A){
  DSDPFunctionBegin;
  if (id>=0 && id < ADATA->nnzmats){
    if (vari) *vari=ADATA->nzmat[id];
    if (scl) *scl=ADATA->scl;
    if (A) *A=ADATA->A[id];
  } else {
    DSDPSETERR2(2,"Invalid Matrix request.  0 <= %d < %d\n",id,ADATA->nnzmats);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockDataRowSparsity"
/*!
\fn int DSDPBlockDataRowSparsity(DSDPBlockData *ADATA,int row, int ai[], int rnnz[],int n);
\brief Determine sparsity pattern of data.
\param ADATA block of data.
\param row from 0 to n.
\param ai (input) array of ones and zeros that identify which data matrices to check.
\param rnnz (output) array of length m where nonzeros indicate nonzero data.
\param n dimension of block.
*/
int DSDPBlockDataRowSparsity(DSDPBlockData *ADATA,int row, int ai[], int rnnz[],int n){
  int info,i,vari,rn;
  DSDPFunctionBegin;
  if (ai){
    for (i=0; i<ADATA->nnzmats; i++){
      vari=ADATA->nzmat[i];
      if (ai[vari]==0){continue;}
      info=DSDPDataMatGetRowNonzeros(ADATA->A[i],row, n, rnnz, &rn); DSDPCHKVARERR(vari,info);
    }
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPBlockRemoveDataMatrix"
/*!
\fn int DSDPBlockRemoveDataMatrix(DSDPBlockData *ADATA, int vari);
\brief Remove a data matrix.
\param ADATA block of data.
\param vari variable corresponding the matrix.
*/
int DSDPBlockRemoveDataMatrix(DSDPBlockData *ADATA,int vari){ 
  int info,ii,k;
  DSDPFunctionBegin;
  for (ii=0;ii<ADATA->nnzmats;ii++){
    if (ADATA->nzmat[ii]==vari){
      info=DSDPDataMatDestroy(&ADATA->A[ii]);DSDPCHKVARERR(vari,info);
      info=DSDPSetDataMatZero(&ADATA->A[ii]);DSDPCHKVARERR(vari,info);
      for (k=ii;k<ADATA->nnzmats;k++){
	ADATA->A[k]=ADATA->A[k+1];
	ADATA->nzmat[k]=ADATA->nzmat[k+1];
      }
      ADATA->nnzmats--;
      info=DSDPSetDataMatZero(&ADATA->A[ADATA->nnzmats]);DSDPCHKERR(info);
      DSDPFunctionReturn(0);
    }
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPBlockAddDataMatrix"
/*!
\fn int DSDPBlockAddDataMatrix(DSDPBlockData *ADATA,int vari, struct DSDPDataMat_Ops* dsdpdataops, void* data){ 

\brief Add data matrix into SDP block
\param ADATA block of data.
\param vari the variable to which the matrix corresponds.
\param dsdpdataops function pointers
\param data opaque pointer to matrix.
*/
int DSDPBlockAddDataMatrix(DSDPBlockData *ADATA,int vari, struct DSDPDataMat_Ops* dsdpdataops, void* data){ 
  int info,ii;
  DSDPFunctionBegin;
  if (ADATA->nnzmats>=ADATA->maxnnzmats){
    info=DSDPBlockDataAllocate(ADATA,2*ADATA->maxnnzmats+7);DSDPCHKERR(info);
  }
  ii=ADATA->nnzmats;
  info=DSDPDataMatDestroy(&ADATA->A[ii]);DSDPCHKERR(info);
  info=DSDPDataMatSetData(&ADATA->A[ii], dsdpdataops, data);DSDPCHKVARERR(vari,info);
  ADATA->nzmat[ii]=vari;
  ADATA->nnzmats++;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPBlockSetDataMatrix"
/*!
\fn int DSDPBlockSetDataMatrix(DSDPBlockData *ADATA,int vari, struct DSDPDataMat_Ops* dsdpdataops, void* data){ 

\brief Set data matrix into SDP block
\param ADATA block of data.
\param vari the variable to which the matrix corresponds.
\param dsdpdataops function pointers
\param data opaque pointer to matrix.
*/
int DSDPBlockSetDataMatrix(DSDPBlockData *ADATA,int vari, struct DSDPDataMat_Ops* dsdpdataops, void* data){ 
  int info;
  DSDPFunctionBegin;
  info=DSDPBlockRemoveDataMatrix(ADATA,vari);DSDPCHKERR(info);
  info=DSDPBlockAddDataMatrix(ADATA,vari,dsdpdataops,data);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockNorm2"
int DSDPBlockNorm2(DSDPBlockData *ADATA, int n){
  double fn2,tt=0;
  int    ii,info;
  DSDPFunctionBegin;
  for (ii=0;ii<ADATA->nnzmats;ii++){
    info=DSDPDataMatFNorm2(ADATA->A[ii],n,&fn2); DSDPCHKERR(info);
    tt+=fn2;
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPBlockANorm2"
int DSDPBlockANorm2(DSDPBlockData *ADATA, DSDPVec ANorm2, int n){

  double fn2,scl=ADATA->scl;
  int    ii,vari,info;

  DSDPFunctionBegin;
  info=DSDPBlockNorm2(ADATA,n);DSDPCHKERR(info);
  scl=ADATA->scl;
  for (ii=0;ii<ADATA->nnzmats;ii++){
    vari=ADATA->nzmat[ii];
    info=DSDPDataMatFNorm2(ADATA->A[ii],n,&fn2); DSDPCHKVARERR(vari,info);
    info=DSDPVecAddElement(ANorm2,vari,fn2*scl);DSDPCHKVARERR(vari,info);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPBlockView"
/*!
\fn int DSDPBlockView(DSDPBlockData *ADATA);

\brief Print the structure of the block.
\param ADATA block of data.
*/
int DSDPBlockView(DSDPBlockData *ADATA){
  int    ii,kk;

  DSDPFunctionBegin;
  for (ii=0;ii<ADATA->nnzmats;ii++){
    kk=ADATA->nzmat[ii];
    if (kk==0){ printf("+ C\n");}
    else { printf(" - A[%d] y%d\n",kk,kk);}
  }
  printf(" = S >= 0\n");
  DSDPFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DSDPBlockView2"
/*!
\fn int DSDPBlockView2(DSDPBlockData *ADATA);

\brief Print the data 
\param ADATA block of data.
*/
int DSDPBlockView2(DSDPBlockData *ADATA){
  int    ii,vari,info;

  DSDPFunctionBegin;
  for (ii=0;ii<ADATA->nnzmats;ii++){
    vari=ADATA->nzmat[ii];
    printf("A[%d] y%d \n",vari,vari);
    info=DSDPDataMatView(ADATA->A[ii]); DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPDataMatCheck"
/*!
\fn int DSDPDataMatCheck(DSDPDataMat AA, SDPConeVec W, DSDPIndex IS, DSDPVMat XX);

\brief Check correctness of operations on the data.
\param AA data matrix.
\param W work vector
\param IS work array
\param XX work array
*/
int DSDPDataMatCheck(DSDPDataMat AA, SDPConeVec W, DSDPIndex IS, DSDPVMat XX){

  double *xx,ack,vAv=0,esum=0,sum,eignorm,fnorm22,dnorm,scl=1;
  int    k,n,nn,rank,info;

  DSDPFunctionBegin;
  info=SDPConeVecGetSize(W,&n);DSDPCHKERR(info);

  info=DSDPVMatZeroEntries(XX);DSDPCHKERR(info);
  info=DSDPDataMatGetRank(AA,&rank,n);DSDPCHKERR(info);
  for (k=0; k<rank; k++){
    info=DSDPDataMatGetEig(AA,k,W,IS,&ack); DSDPCHKERR(info);
    info=SDPConeVecDot(W,W,&eignorm);DSDPCHKERR(info);
    info=DSDPVMatAddOuterProduct(XX,scl*ack,W);DSDPCHKERR(info);
    info=DSDPDataMatVecVec(AA,W,&sum);DSDPCHKERR(info);
    vAv+=ack*ack*eignorm*eignorm*scl;
  }
  info=DSDPDataMatFNorm2(AA,n,&fnorm22); DSDPCHKERR(info);
  
  info=DSDPVMatScaleDiagonal(XX,0.5); DSDPCHKERR(info);
  info=DSDPVMatGetArray(XX, &xx, &nn); DSDPCHKERR(info);
  info=DSDPDataMatDot(AA,xx,nn,n,&esum);DSDPCHKERR(info);
  info=DSDPVMatRestoreArray(XX, &xx, &nn); DSDPCHKERR(info);
  info=DSDPVMatScaleDiagonal(XX,2.0); DSDPCHKERR(info);
  
  info=DSDPVMatGetArray(XX, &xx, &nn); DSDPCHKERR(info);
  info=DSDPDataMatAddMultiple(AA,-1.0,xx,nn,n); DSDPCHKERR(info);
  info=DSDPVMatRestoreArray(XX, &xx, &nn); DSDPCHKERR(info);
  if (0==1){info=DSDPVMatView(XX);DSDPCHKERR(info);}
  info=DSDPVMatNormF2(XX,&dnorm); DSDPCHKERR(info);
  printf("  %4.4e, %4.4e  %4.4e\n",esum,vAv,fnorm22);
  printf("  error1: %4.4e, error2: %4.4e,  error3: %4.4e\n",sqrt(dnorm),fabs(esum-vAv),fabs(fnorm22-vAv));
  if (dnorm>1)  printf("Check Add or eigs\n");
  if (fabs(esum-vAv) > 1.0) printf("Check vAv \n");
  if (fabs(fnorm22-vAv) > 1.0) printf("Check fnorm22\n");
  
  DSDPFunctionReturn(0);
}

