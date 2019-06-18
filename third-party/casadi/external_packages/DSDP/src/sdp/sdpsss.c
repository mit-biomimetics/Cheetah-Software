#include "dsdpsdp.h"
#include "dsdpsys.h"
/*!
\file sdpsss.c
\brief Determine the sparsity of each block, and choose suitable dual, Delta S, and X matrix structures.
*/

int DSDPCreateS(DSDPBlockData*,char,int,DSDPVec,DSDPVMat,SDPConeVec, SDPConeVec,DSDPDualMat*, DSDPDualMat*, DSDPDSMat*, void*);

static int DSDPCreateDS(DSDPBlockData*, DSDPVMat,int*,int,int,int,int*,int*,DSDPDSMat*);
static int DSDPCreateDS2(DSDPBlockData*, DSDPVMat,int*,int,int,int,int*,int*,DSDPDSMat*);

static int DSDPCreateS1(DSDPBlockData*,int,DSDPVec,DSDPVMat,SDPConeVec, SDPConeVec,DSDPDualMat*, DSDPDualMat*, DSDPDSMat*, void*);
static int DSDPCreateS2(DSDPBlockData*,int,DSDPVec,DSDPVMat,SDPConeVec, SDPConeVec,DSDPDualMat*, DSDPDualMat*, DSDPDSMat*, void*);

extern int DSDPXMatPCreate(int,struct DSDPVMat_Ops**,void**);
extern int DSDPXMatPCreateWithData(int,double[],int,struct DSDPVMat_Ops**,void**);
extern int DSDPXMatUCreate(int,struct DSDPVMat_Ops**,void**);
extern int DSDPXMatUCreateWithData(int,double[],int,struct DSDPVMat_Ops**,void**);

extern int DSDPLAPACKSUDualMatCreate2(int,struct DSDPDualMat_Ops**,void**, struct DSDPDualMat_Ops**,void**);
extern int DSDPLAPACKSUDualMatCreate2P(int,struct DSDPDualMat_Ops**,void**, struct DSDPDualMat_Ops**,void**);
extern int DSDPLAPACKPUDualMatCreate2(int,struct DSDPDualMat_Ops**,void**, struct DSDPDualMat_Ops**,void**);
extern int DSDPDiagDualMatCreateP(int,struct DSDPDualMat_Ops**,void**, struct DSDPDualMat_Ops**,void**);
extern int DSDPDiagDualMatCreateU(int,struct DSDPDualMat_Ops**,void**, struct DSDPDualMat_Ops**,void**);
extern int DSDPDenseDualMatCreate(int, char,struct DSDPDualMat_Ops**,void**, struct DSDPDualMat_Ops**,void**);
extern int DSDPSparseDualMatCreate(int, int*, int*, int,char,int*,struct DSDPDualMat_Ops**,void**, struct DSDPDualMat_Ops**,void**);

extern int DSDPSparseMatCreatePattern2P(int,int[],int[],int,struct DSDPDSMat_Ops**,void**);
extern int DSDPSparseMatCreatePattern2U(int,int[],int[],int,struct DSDPDSMat_Ops**,void**);

extern int DSDPCreateDiagDSMatP(int,struct DSDPDSMat_Ops**,void**);
extern int DSDPCreateDiagDSMatU(int,struct DSDPDSMat_Ops**,void**);

extern int DSDPCreateDSMatWithArray(int,double[],int,struct DSDPDSMat_Ops**,void**);
extern int DSDPCreateDSMatWithArray2(int,double[],int, struct DSDPDSMat_Ops**,void**);

extern int DSDPCreateURDSMat(int,struct DSDPDSMat_Ops**,void**);

/*
#undef __FUNCT__
#define __FUNCT__ "DSDPSetDualMatrix"
int DSDPSetDualMatrix(SDPCone sdpcone,int (*createdualmatrix)(DSDPBlockData*,DSDPVec,DSDPVMat,SDPConeVec,SDPConeVec,DSDPDualMat*, DSDPDualMat*, DSDPDSMat*, void*),void*ctx){
  DSDPFunctionBegin;
  sdpcone->createdualmatrix=createdualmatrix;
  sdpcone->dualmatctx=ctx;
  DSDPFunctionReturn(0);
}
*/
#undef __FUNCT__
#define __FUNCT__ "CountNonzeros"
static int CountNonzeros(DSDPBlockData *ADATA,int m, int rnnz[], int innz[],int n,int *nnz1, int *nnz2)
{
  int i,j,info,totalnnz1=0,totalnnz2=0;
  
  for (i=0;i<n;i++){
    memset(rnnz,0,n*sizeof(int));
    /* SParsity pattern for DS only requires the constraint matrices A and residual */
    for (j=0;j<m;j++) innz[j]=1;innz[0]=0;
    info=DSDPBlockDataRowSparsity(ADATA,i,innz,rnnz,n);DSDPCHKERR(info);
    for (j=0; j<i; j++){ if (rnnz[j]>0) totalnnz1++;}
    /* Adjacency pattern for S also requires the objective matrix C */
    for (j=0;j<m;j++) innz[j]=0;innz[0]=1;
    info=DSDPBlockDataRowSparsity(ADATA,i,innz,rnnz,n);DSDPCHKERR(info);
    for (j=0; j<i; j++){ if (rnnz[j]>0) totalnnz2++;}
  }
  *nnz1=totalnnz1;
  *nnz2=totalnnz2;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "CreateS1b"
static int CreateS1b(DSDPBlockData *ADATA, int innz[], int m, int n, int tnnz[], int rnnz[], int snnz[])
{
  int i,j,info;
  if (ADATA->nnzmats<=0){return 0;};

  memset(rnnz,0,n*sizeof(int));
  for (i=0;i<m;i++) innz[i]=1;
  innz[0]=0;

  /* Check matrices A for nonzero entries. */
  for (i=0;i<n;i++){
    memset(tnnz,0,n*sizeof(int));
    info=DSDPBlockDataRowSparsity(ADATA,i,innz,tnnz,n);DSDPCHKERR(info);
    for (j=0; j<=i; j++){
      if (tnnz[j]>0){ *snnz=j; snnz++; rnnz[i]++;}
    }
  }
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPCreateDS"
int DSDPCreateDS(DSDPBlockData *ADATA, DSDPVMat T, int iworkm[], int m, int n, int nnz, int rnnz[], int tnnz[], DSDPDSMat *B){
  int i,n1,*cols,allnnz,info;
  double *ss;
  void *dsmat;
  struct DSDPDSMat_Ops* dsops;
  DSDPFunctionBegin;
  DSDPLogInfo(0,19,"DS Matrix has %d nonzeros of %d\n",nnz,n*(n-1)/2);
  if (nnz==0){
    info=DSDPCreateDiagDSMatP(n,&dsops,&dsmat); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Using Diagonal Delta S matrix\n");
  } else if (2*nnz +n < n*n/5){
    info=DSDPVMatGetArray(T,&ss,&n1); DSDPCHKERR(info);
    cols=(int*)ss;
    info = CreateS1b(ADATA, iworkm, m, n, tnnz, rnnz, cols); DSDPCHKERR(info);
    for (allnnz=0,i=0;i<n;i++){allnnz+=rnnz[i];}
    info = DSDPSparseMatCreatePattern2P(n,rnnz,cols,allnnz,&dsops,&dsmat); DSDPCHKERR(info);
    info=DSDPVMatRestoreArray(T,&ss,&n1); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Using Sparse Delta S matrix\n");
  } else {
    info=DSDPVMatGetArray(T,&ss,&n1); DSDPCHKERR(info);
    info=DSDPCreateDSMatWithArray(n,ss,n1,&dsops,&dsmat); DSDPCHKERR(info);
    info=DSDPVMatRestoreArray(T,&ss,&n1); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Using Full Delta S matrix\n");
  }
  info=DSDPDSMatSetData(B,dsops,dsmat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CreateS1c"
static int CreateS1c(DSDPBlockData *ADATA, int innz[], int m, int n, int tnnz[], int rnnz[], int snnz[])
{
  int i,j,info;
  memset(rnnz,0,n*sizeof(int));
  for (i=0;i<m;i++) innz[i]=1;
  /* Check matrix C and A for nonzero entries. */
  for (i=0;i<n;i++){
    memset(tnnz,0,n*sizeof(int));
    info=DSDPBlockDataRowSparsity(ADATA,i,innz,tnnz,n);DSDPCHKERR(info);
    for (j=i+1; j<n; j++){
      if (tnnz[j]>0){ *snnz=j; snnz++; rnnz[i]++;}
    } 
  }
  return 0;
}

static int dsdpuselapack=1;
#undef __FUNCT__  
#define __FUNCT__ "SDPConeUseLAPACKForDualMatrix"
int SDPConeUseLAPACKForDualMatrix(SDPCone sdpcone,int flag){ 
  DSDPFunctionBegin;
  dsdpuselapack = flag;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPCreateS"
static int DSDPCreateS1(DSDPBlockData *ADATA, int trank, DSDPVec WY, DSDPVMat T, SDPConeVec W1, SDPConeVec W2, DSDPDualMat *S, DSDPDualMat *SS, DSDPDSMat *DS, void*ctx){

  int nnz,n,n1,*cols,*rnnz,*tnnz,*iworkm,m,info;
  int dsnnz,snnz,sfnnz;
  double *pss;
  void *smat1,*smat2;
  struct DSDPDualMat_Ops *sops1,*sops2;

  DSDPFunctionBegin;
  info=DSDPVecGetSize(WY,&m);DSDPCHKERR(info); 
  info=DSDPVecGetArray(WY,&pss);DSDPCHKERR(info); 
  iworkm=(int*)pss;

  info=SDPConeVecGetSize(W1,&n);DSDPCHKERR(info); 
  info=SDPConeVecGetArray(W1,&pss);DSDPCHKERR(info); 
  tnnz=(int*)pss;
  info=SDPConeVecGetArray(W2,&pss);DSDPCHKERR(info); 
  rnnz=(int*)pss;

  DSDPLogInfo(0,19,"Compute Sparsity\n");
  info = CountNonzeros(ADATA, m, rnnz, iworkm, n, &dsnnz,&snnz); DSDPCHKERR(info);
  nnz=snnz;
  /* printf("Nonzeros: %d %d of %d\n",dsnnz,snnz,n*(n-1)/2); */
  /* TT and DS could share memory */
  info=DSDPCreateDS(ADATA,T,iworkm,m,n,dsnnz,rnnz,tnnz,DS);DSDPCHKERR(info);

  if (nnz==0){
    info=DSDPDiagDualMatCreateP(n,&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Using Diagonal S matrix\n");
  } else if (2*snnz+n+2<n*n/8){
    info=DSDPVMatGetArray(T,&pss,&n1); DSDPCHKERR(info);
    cols=(int*)pss;
    info = CreateS1c(ADATA, iworkm, m, n, tnnz, rnnz, cols); DSDPCHKERR(info);
    info=DSDPSparseDualMatCreate(n,rnnz,cols,trank,'P',&sfnnz,&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    info=DSDPVMatRestoreArray(T,&pss,&n1); DSDPCHKERR(info);
    /*   printf("NNZ: %d %d of %d\n",2*snnz+n,sfnnz*2+n,n*n); */
    DSDPLogInfo(0,19,"Count %d of %d nonzeros: Using Sparse S matrix\n",nnz,n*(n-1)/2);
    DSDPLogInfo(0,19,"Total rank of block: %d, n= %d\n",trank,n);
  } else if (n>20 && dsdpuselapack){
    info=DSDPLAPACKSUDualMatCreate2P(n,&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Count %d of %d nonzeros: Using Full Dense LAPACK S matrix\n",nnz,n*(n-1)/2);
  } else if (dsdpuselapack){
    info=DSDPLAPACKPUDualMatCreate2(n,&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Count %d of %d nonzeros: Using Packed Dense LAPACK S matrix\n",nnz,n*(n-1)/2);
  } else {
    info=DSDPDenseDualMatCreate(n,'P',&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Count %d of %d nonzeros: Using Dense S matrix\n",nnz,n*(n-1)/2);
  }
  info=DSDPDualMatSetData(S,sops1,smat1); DSDPCHKERR(info);
  info=DSDPDualMatSetData(SS,sops2,smat2); DSDPCHKERR(info);

  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DSDPCreateDS2"
int DSDPCreateDS2(DSDPBlockData *ADATA, DSDPVMat T, int iworkm[], int m, int n, int nnz, int rnnz[], int tnnz[], DSDPDSMat *B){
  int i,n1,*cols,allnnz,info;
  double *ss;
  void *dsmat;
  struct DSDPDSMat_Ops* dsops;
  DSDPFunctionBegin;
  DSDPLogInfo(0,19,"DS Matrix has %d nonzeros of %d\n",nnz,n*(n-1)/2);
  if (nnz==0){
    info=DSDPCreateDiagDSMatU(n,&dsops,&dsmat); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Using Diagonal Delta S matrix\n");
  } else if (2*nnz +n < n*n/4){
    info=DSDPVMatGetArray(T,&ss,&n1); DSDPCHKERR(info);
    cols=(int*)ss;
    info = CreateS1b(ADATA, iworkm, m, n, tnnz, rnnz, cols); DSDPCHKERR(info);
    for (allnnz=0,i=0;i<n;i++){allnnz+=rnnz[i];}
    info = DSDPSparseMatCreatePattern2U(n,rnnz,cols,allnnz,&dsops,&dsmat); DSDPCHKERR(info);
    info=DSDPVMatRestoreArray(T,&ss,&n1); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Using Sparse Delta S matrix\n");
  } else {
    info=DSDPVMatGetArray(T,&ss,&n1); DSDPCHKERR(info);
    info=DSDPCreateDSMatWithArray2(n,ss,n1,&dsops,&dsmat); DSDPCHKERR(info);
    info=DSDPVMatRestoreArray(T,&ss,&n1); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Using Full Delta S matrix\n");
  }
  info=DSDPDSMatSetData(B,dsops,dsmat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateS2"
static int DSDPCreateS2(DSDPBlockData *ADATA, int trank, DSDPVec WY, DSDPVMat T, SDPConeVec W1, SDPConeVec W2, DSDPDualMat *S, DSDPDualMat *SS, DSDPDSMat *DS, void*ctx){

  int nnz,n,n1,*cols,*rnnz,*tnnz,*iworkm,m,info;
  int dsnnz,snnz,sfnnz;
  double *pss;
  void *smat1,*smat2;
  struct DSDPDualMat_Ops *sops1,*sops2;
  DSDPFunctionBegin;

  info=DSDPVecGetSize(WY,&m);DSDPCHKERR(info); 
  info=DSDPVecGetArray(WY,&pss);DSDPCHKERR(info); 
  iworkm=(int*)pss;

  info=SDPConeVecGetSize(W1,&n);DSDPCHKERR(info); 
  info=SDPConeVecGetArray(W1,&pss);DSDPCHKERR(info); 
  tnnz=(int*)pss;
  info=SDPConeVecGetArray(W2,&pss);DSDPCHKERR(info); 
  rnnz=(int*)pss;

  DSDPLogInfo(0,19,"Compute Sparsity\n");
  info = CountNonzeros(ADATA, m, rnnz, iworkm, n, &dsnnz,&snnz); DSDPCHKERR(info);
  nnz=snnz;
  /* printf("Nonzeros: %d %d of %d\n",dsnnz,snnz,n*(n-1)/2); */
  /* TT and DS could share memory */
  info=DSDPCreateDS2(ADATA,T,iworkm,m,n,dsnnz,rnnz,tnnz,DS);DSDPCHKERR(info);

  if (nnz==0){
    info=DSDPDiagDualMatCreateU(n,&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Using Diagonal S matrix\n");
  } else if (2*snnz+n+2<n*n/10){
    info=DSDPVMatGetArray(T,&pss,&n1); DSDPCHKERR(info);
    cols=(int*)pss;
    info = CreateS1c(ADATA, iworkm, m, n, tnnz, rnnz, cols); DSDPCHKERR(info);
    info=DSDPSparseDualMatCreate(n,rnnz,cols,trank,'U',&sfnnz,&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    info=DSDPVMatRestoreArray(T,&pss,&n1); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Count %d of %d nonzeros: Using Sparse S matrix\n",nnz,n*(n-1)/2);
  } else if (dsdpuselapack){
    info=DSDPLAPACKSUDualMatCreate2(n,&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Count %d of %d nonzeros: Using Full Dense LAPACK S matrix\n",nnz,n*(n-1)/2);
  } else {
    info=DSDPDenseDualMatCreate(n,'U',&sops1,&smat1,&sops2,&smat2); DSDPCHKERR(info);
    DSDPLogInfo(0,19,"Count %d of %d nonzeros: Using Packed Dense S matrix\n",nnz,n*(n-1)/2);
  }
  info=DSDPDualMatSetData(S,sops1,smat1); DSDPCHKERR(info);
  info=DSDPDualMatSetData(SS,sops2,smat2); DSDPCHKERR(info);

  DSDPFunctionReturn(0);
}


int DSDPCreateS(DSDPBlockData*,char,int,DSDPVec,DSDPVMat,SDPConeVec,SDPConeVec,DSDPDualMat*, DSDPDualMat*, DSDPDSMat*, void*);

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateS"
/*!
\fn int DSDPCreateS(DSDPBlockData *ADATA, char UPLQ, int trank, DSDPVec WY, DSDPVMat T, SDPConeVec W1, SDPConeVec W2, DSDPDualMat *S, DSDPDualMat *SS, DSDPDSMat *DS, void*ctx);

\brief Create S1, S2, and DS
\param ADATA semidefinite block of data.
\param UPLQ such as packed symmetric or upper full symmetric
\param trank rank of data in block
\param WY Work vector
\param T Work matrix.
\param W1 Work vector
\param W2 Work vector
\param S New dual matrix.
\param SS New dual matrix.
\param DS New DS matrix.
\param ctx pointer to structure.
*/
int DSDPCreateS(DSDPBlockData *ADATA, char UPLQ, int trank, DSDPVec WY, DSDPVMat T, SDPConeVec W1, SDPConeVec W2, DSDPDualMat *S, DSDPDualMat *SS, DSDPDSMat *DS, void*ctx){
  int info;
  DSDPFunctionBegin;
  switch (UPLQ){
  case 'P':
    info=DSDPCreateS1(ADATA,trank,WY,T,W1,W2,S,SS,DS,ctx);DSDPCHKERR(info); 
    break;
  case 'U':
    info=DSDPCreateS2(ADATA,trank,WY,T,W1,W2,S,SS,DS,ctx);DSDPCHKERR(info); 
    break;
  }
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DSDPCreateS"
int DSDPUseDefaultDualMatrix(SDPCone sdpcone){
  DSDPFunctionBegin;
  /*
  int info;
  info=DSDPSetDualMatrix(sdpcone,DSDPCreateS2,0); DSDPCHKERR(info);
  */
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPMakeVMat"
/*!
\fn int DSDPMakeVMat(char UPLQ, int n, DSDPVMat *X);
\brief Allocate V matrix.
\param UPLQ format
\param n dimension
\param X new matrix.
\sa SDPConeSetStorageFormat()
*/
int DSDPMakeVMat(char UPLQ, int n, DSDPVMat *X){
  int info;
  void *xmat;
  struct DSDPVMat_Ops* xops;
  DSDPFunctionBegin;
  switch (UPLQ){
  case 'P':
    info=DSDPXMatPCreate(n,&xops,&xmat);DSDPCHKERR(info); 
    break;
  case 'U':
    info=DSDPXMatUCreate(n,&xops,&xmat);DSDPCHKERR(info); 
    break;
  }
  info=DSDPVMatSetData(X,xops,xmat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPMakeVMatWithArray"
/*!
\fn int DSDPMakeVMatWithArray(char UPLQ, double xx[], int nnz, int n, DSDPVMat *X);
\brief Allocate V matrix using the given array.
\param UPLQ format
\param xx array
\param nnz length of the array
\param n dimension
\param X new matrix.
\sa SDPConeSetStorageFormat()
\sa SDPConeGetXArray()
*/
int DSDPMakeVMatWithArray(char UPLQ, double xx[], int nnz, int n, DSDPVMat *X){
  int info;
  void *xmat;
  struct DSDPVMat_Ops* xops;
  DSDPFunctionBegin;
  switch (UPLQ){
  case 'P':
    info=DSDPXMatPCreateWithData(n,xx,nnz,&xops,&xmat);DSDPCHKERR(info); 
    break;
  case 'U':
    info=DSDPXMatUCreateWithData(n,xx,nnz,&xops,&xmat);DSDPCHKERR(info); 
    break;
  }
  info=DSDPVMatSetData(X,xops,xmat); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
