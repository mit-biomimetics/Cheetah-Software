#include "numchol.h"
#include "dsdpschurmat_impl.h"
#include "dsdpsys.h"
#include "dsdpvec.h"
#include "dsdp5.h"

/*!
\file cholmat.c 
\brief Sparse Cholesky for Schur complement matrix
 */
int DSDPSparsityInSchurMat(DSDP, int, int *, int);
int DSDPSetSchurMatOps(DSDP,struct DSDPSchurMat_Ops*,void*);

typedef struct {
  chfac *M;
  int m;
  int is_dense;
  int *rnnz;
  int *colnnz;
  int nnz;
  DSDPVec    D1;
  DSDP     dsdp;
} MCholSolverALL;


static int dsdpuselapack=1;
#undef __FUNCT__  
#define __FUNCT__ "DSDPUseLAPACKForSchur"
int DSDPUseLAPACKForSchur(DSDP dsdp,int flag){ 
  DSDPFunctionBegin;
  dsdpuselapack = flag;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Taddline"
static int Taddline(void* M, int row, double dd, double v[], int m){
  int info;
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  DSDPFunctionBegin;
  info=MatAddColumn4(AMA->M,dd,v,row);  DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Tadddiagonal"
static int Tadddiagonal(void* M, int row, double v){
  int info;
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  DSDPFunctionBegin;
  info=MatAddDiagonalElement(AMA->M,row,v);  DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Tassemble"
static int Tassemble(void*M){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPCheckForSparsity"
static int DSDPCheckForSparsity( DSDP dsdp, int m, int *rnnz, int tnnz[], int *totalnnz){
  int info,i,j,tottalnnz=0;
  DSDPFunctionBegin;
  memset(rnnz,0,(m+1)*sizeof(int));
  for (i=0;i<m;i++){
    info=DSDPSparsityInSchurMat(dsdp,i,tnnz,m); DSDPCHKERR(info);
    for (j=i+1;j<m;j++){ if (tnnz[j]>0) rnnz[i+1]++;} 
  }

  for (i=0;i<m;i++){tottalnnz+=rnnz[i+1];}
  *totalnnz=tottalnnz;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateM"
static int DSDPCreateM(MCholSolverALL *ABA, chfac **M, int rrnnz[],int tnnz[], int totalnnz){

  int *snnz,*rnnz;
  int i,j,tt,info;
  int col,*perm,k;
  chfac * sfptr;
  int n=ABA->m,m=ABA->m;
  DSDP dsdp=ABA->dsdp;
 
  DSDPFunctionBegin;

  DSDPCALLOC2(&snnz,int,(totalnnz+1),&info); DSDPCHKERR(info);
  DSDPCALLOC2(&rnnz,int,(m+1),&info); DSDPCHKERR(info);
  for (i=0;i<=m;i++){ rnnz[i]=rrnnz[i];}
  tt=0;
  for (i=0;i<m;i++){
    info=DSDPSparsityInSchurMat(dsdp,i,tnnz,m); DSDPCHKERR(info);
    for (j=i+1;j<m;j++){ if (tnnz[j]>0) {snnz[tt]=j; tt++;} } 
  }

  printf("Trying Sparse M: Total nonzeros: %d of %d \n",totalnnz,m*(m-1)/2 );
  /* Create sparse dual matrix structure */
  SymbProc(rnnz+1,snnz,n,&sfptr);
  ABA->is_dense=0;
  ABA->M=sfptr;
  ABA->nnz=totalnnz;
  ABA->rnnz=rnnz;
  ABA->colnnz=snnz;
  *M=sfptr;

  for (i=0;i<m;i++){
    rnnz[i+1]+=rnnz[i];
  }

  perm=sfptr->invp;
  for (i=m-1;i>=0;i--){
    for (j=rnnz[i+1]-1;j>=rnnz[i];j--){
      col=snnz[j];
      tt=perm[col];
      if (perm[i] > tt ){
	for (k=j;k<rnnz[col]-1;k++){ snnz[k]=snnz[k+1];}
	for (k=i;k<col;k++){ rnnz[k+1]--;}
	snnz[rnnz[col]]=i;
      }
    }
  }
  
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPLinearSolverPrepare"
static int DSDPLinearSolverPrepare(void* ctx,int*flag){

  cfc_sta Cfact;
  chfac *sf;
  MCholSolverALL *AMA = (MCholSolverALL *)ctx; 

  DSDPFunctionBegin;
  *flag=0;
  sf=AMA->M;  
  /*
  Cfact=(cfc_sta)ChlFact(sf,sf->iw,sf->rw,FALSE);
  */
  Cfact=(cfc_sta)ChlFact(sf,sf->iw,sf->rw,TRUE);
  if (CfcOk!=Cfact ){ *flag=1; /*  printf("Not Good factorization \n"); */ }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPLinearSolve"
static int DSDPLinearSolve2(void* ctx, double dd[], double xx[], int n){

  int i,info;
  double *bb;
  MCholSolverALL *AMA = (MCholSolverALL *)ctx; 

  DSDPFunctionBegin;
  info=DSDPVecGetArray(AMA->D1, &bb);DSDPCHKERR(info);
  for (i=0;i<n;i++){ bb[i]=dd[i];}
  ChlSolve(AMA->M, bb, xx);
  info=DSDPVecRestoreArray(AMA->D1, &bb);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPGramMatRowNonzeros"
static int DSDPGramMatRowNonzeros(void*M, int row, double cols[], int *ncols, int nrows){
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  int i;
  DSDPFunctionBegin;
  if (AMA->is_dense){
    *ncols = nrows-row;
    for (i=row;i<nrows;i++){
      cols[i]=1.0;
    }
  } else {
    *ncols = AMA->rnnz[row+1] - AMA->rnnz[row]+1;
    cols[row]=1.0;
    for (i=AMA->rnnz[row]; i<AMA->rnnz[row+1]; i++){
      cols[AMA->colnnz[i]]=1.0;
    }
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Tzero"
static int Tzero(void*M){
  int info;
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  DSDPFunctionBegin;
  info=MatZeroEntries4(AMA->M); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Tdestroy"
static int Tdestroy(void*M){
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  int info;
  DSDPFunctionBegin;
  CfcFree(&AMA->M);
  info = DSDPVecDestroy(&AMA->D1); DSDPCHKERR(info);
  if (AMA->is_dense==0 && AMA->rnnz ){
    DSDPFREE(&AMA->rnnz,&info);DSDPCHKERR(info);
    DSDPFREE(&AMA->colnnz,&info);DSDPCHKERR(info);
  }
  DSDPFREE(&AMA,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int Tsetup(void*M, int m){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

static int TTTMatMult(void*M,double x[],double y[],int n){
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  int info;
  DSDPFunctionBegin;
  info=MatMult4(AMA->M,x,y,n); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int TTTMatShiftDiagonal(void *M, double d){
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  int info;
  DSDPFunctionBegin;
  info=Mat4DiagonalShift(AMA->M,d); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int TTTMatAddDiagonal(void *M, double diag[], int m){
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  int info;
  DSDPFunctionBegin;
  info=Mat4AddDiagonal(AMA->M,diag,m); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int TTTMatView(void *M){
  MCholSolverALL *AMA = (MCholSolverALL *)M; 
  int info;
  DSDPFunctionBegin;
  info=Mat4View(AMA->M); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


/*
static int TTTMatGetDiagonal(void *M, double d[], int n){
  chfac*A = (chfac*)M;
  int info;
  DSDPFunctionBegin;
  info=Mat4GetDiagonal(A,d,n); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
*/
static const char* tmatname="SPARSE PSD";

static int TMatOpsInit(struct  DSDPSchurMat_Ops *mops){  
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatOpsInitialize(mops); DSDPCHKERR(info);
  mops->matrownonzeros=DSDPGramMatRowNonzeros;
  mops->mataddrow=Taddline;
  mops->matadddiagonal=TTTMatAddDiagonal;
  mops->mataddelement=Tadddiagonal;
  mops->matshiftdiagonal=TTTMatShiftDiagonal;
  mops->matassemble=Tassemble;
  mops->matscaledmultiply=TTTMatMult;
  mops->matfactor=DSDPLinearSolverPrepare;
  mops->matsolve=DSDPLinearSolve2;
  mops->matdestroy=Tdestroy;
  mops->matzero=Tzero;
  mops->matsetup=Tsetup;
  mops->matview=TTTMatView;
  mops->id=5;
  mops->matname=tmatname;
  DSDPFunctionReturn(0);
}

static struct  DSDPSchurMat_Ops dsdpmatops;

int DSDPGetDiagSchurMat(int, struct  DSDPSchurMat_Ops **,void **);
int DSDPGetLAPACKSUSchurOps(int,struct DSDPSchurMat_Ops**,void**);
int DSDPGetLAPACKPUSchurOps(int,struct DSDPSchurMat_Ops**,void**);

static int DSDPCreateM(MCholSolverALL*,chfac**,int[],int[],int);
static int DSDPCreateSchurMatrix(void*,int);

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateSchurMatrix"
static int DSDPCreateSchurMatrix(void *ctx, int m){

  int info;
  int *rnnz,*tnnz,totalnnz;
  int gotit=0;
  DSDP dsdp=(DSDP)ctx;
  chfac *sf;
  MCholSolverALL *AMA;
  void *tdata;
  struct  DSDPSchurMat_Ops *tmatops;

  DSDPFunctionBegin;
  if (m<=1){
    info=DSDPGetDiagSchurMat(m,&tmatops,&tdata); DSDPCHKERR(info);
    info=DSDPSetSchurMatOps(dsdp,tmatops,tdata); DSDPCHKERR(info);
    return 0;
  }

  DSDPCALLOC2(&rnnz,int,(m+1),&info); DSDPCHKERR(info);
  DSDPCALLOC2(&tnnz,int,(m+1),&info); DSDPCHKERR(info);

  info=DSDPCheckForSparsity(dsdp,m,rnnz,tnnz,&totalnnz); DSDPCHKERR(info);

  if (totalnnz*2+m > m*m*0.1 && dsdpuselapack) {
    gotit=1;
    info=DSDPGetLAPACKSUSchurOps(m,&tmatops,&tdata); 
    /* DSDPCHKERR(info); */ if (info) {gotit=0;printf("Try packed format\n"); }
    DSDPLogInfo(0,8,"Creating Dense Full LAPACK Schur Matrix\n");
    info=DSDPSetSchurMatOps(dsdp,tmatops,tdata); DSDPCHKERR(info);
  } 

  if ( 0 && totalnnz*2+m > m*m*0.1 && dsdpuselapack) {

    info=DSDPGetLAPACKPUSchurOps(m,&tmatops,&tdata); DSDPCHKERR(info);
    DSDPLogInfo(0,8,"Creating Dense Packed LAPACK Schur Matrix\n");
    info=DSDPSetSchurMatOps(dsdp,tmatops,tdata); DSDPCHKERR(info);
    gotit=1;

  } 
  if (gotit==0){

    DSDPCALLOC1(&AMA,MCholSolverALL,&info);DSDPCHKERR(info);
    AMA->dsdp=dsdp; AMA->m=m;
    info=DSDPVecCreateSeq(m,&AMA->D1); DSDPCHKERR(info);
    if (totalnnz*2+m > m*m * 0.11 ){
      info=MchlSetup2(m,&sf); DSDPCHKERR(info);
      AMA->M=sf;      AMA->is_dense=1;
      AMA->rnnz=0;    AMA->colnnz=0;
      DSDPLogInfo(0,8,"Creating Dense Full non LAPACK Schur Matrix\n");
    } else {
      info=DSDPCreateM(AMA,&sf,rnnz,tnnz,totalnnz); DSDPCHKERR(info);
      DSDPLogInfo(0,8,"Creating Sparse Schur Matrix\n");
    }
    AMA->M=sf;  
    info=TMatOpsInit(&dsdpmatops); DSDPCHKERR(info);
    info=DSDPSetSchurMatOps(dsdp,&dsdpmatops,(void*)AMA); DSDPCHKERR(info);
  }
  DSDPFREE(&tnnz,&info);DSDPCHKERR(info);
  DSDPFREE(&rnnz,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static struct  DSDPSchurMat_Ops dsdpmatops000;

#undef __FUNCT__
#define __FUNCT__ "DSDPSetDefaultSchurMatrixStructure"
int DSDPSetDefaultSchurMatrixStructure(DSDP dsdp){
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatOpsInitialize(&dsdpmatops000); DSDPCHKERR(info);
  dsdpmatops000.matsetup=DSDPCreateSchurMatrix;
  info=DSDPSetSchurMatOps(dsdp,&dsdpmatops000,(void*)dsdp);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
 
