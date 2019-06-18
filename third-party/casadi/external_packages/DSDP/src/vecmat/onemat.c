#include "dsdpdatamat_impl.h"
#include "dsdpsys.h"
/*!
\file onemat.c
\brief DSDPDataMat object such that all elements are the same value.
*/

typedef struct {
  double cnst;
  char UPLQ;
  int n;
} cmat;

static int ConstMatDestroy(void*);
static int ConstMatView(void*);
static int ConstMatVecVec(void*, double[], int, double *);
static int ConstMatDot(void*, double[],int,int,double *);
static int ConstMatGetRank(void*, int*, int);
static int ConstMatFactor(void*);
static int ConstMatGetEig(void*, int, double*, double[], int,int[],int*);
static int ConstMatRowNnz(void*, int, int[], int*, int);
static int ConstMatAddRowMultiple(void*, int, double, double[], int);
static int ConstMatAddMultiple(void*, double, double[], int,int);
static int ConstMatTest(void*);

static struct  DSDPDataMat_Ops constantmatops;
static int ConstMatOpsInitialize(struct  DSDPDataMat_Ops*);

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetConstantMat"
int DSDPGetConstantMat(int n, double value, char UPLQ, struct  DSDPDataMat_Ops**mops, void**mmat){ 
  int info;
  cmat*AA;
  DSDPFunctionBegin;
  AA=(cmat*) malloc(1*sizeof(cmat));
  if (AA==NULL) return 1;
  AA->cnst=value;
  AA->n=n;
  AA->UPLQ=UPLQ;
  info=ConstMatOpsInitialize(&constantmatops); if(info){return 1;}
  if (mops){*mops=&constantmatops;}
  if (mmat){*mmat=(void*)AA;}
  DSDPFunctionReturn(0);
}


static int ConstMatDot(void* A, double x[], int nn, int n, double *v){

  cmat* AA = (cmat*)A;
  double dtmp=0.0;
  int i,j;

  for (i=0;i<n;i++){
    for (j=0;j<=i;j++){
      dtmp+= (x[j]);
    }
    if (AA->UPLQ=='U'){
      x=x+n;
    } else {
      x=x+i+1;
    }
  }

  *v=2*dtmp*AA->cnst;
  return 0;
}

static int ConstMatVecVec(void* A, double x[], int n, double *v){

  cmat* AA = (cmat*)A;
  double dtmp=0.0;
  int i;

  for (i=0; i<n; i++){
    dtmp+=x[i];
  }
  *v=dtmp*dtmp*AA->cnst;
  return 0;
}

static int ConstMatAddMultiple(void*A, double dd, double vv[], int nn, int n){
  cmat* AA = (cmat*)A;
  int i,j;
  double ddd=dd*AA->cnst;
  for (i=0;i<n;i++){
    for (j=0;j<i;j++){
      (vv[j])+=ddd;
    }
    vv[i]+=ddd;
    if (AA->UPLQ=='U'){
      vv=vv+n;
    } else {
      vv=vv+i+1;
    }
  }
  return 0;
}

static int ConstMatAddRowMultiple(void*A, int nrow, double dd, double row[], int n){
  cmat* AA = (cmat*)A;
  int i;
  double ddd=dd*AA->cnst;
  for (i=0;i<n;i++){
    row[i] += ddd;
  }
  row[nrow] -= ddd;
 return 0;
}



static int ConstMatFactor(void*A){
  return 0;
}

static int ConstMatGetRank(void *A, int*rank, int n){
  *rank=1;
  return 0;
}

static int ConstMatGetEig(void*A, int neig, double *eig, double v[], int n, int  indx[], int*nind){
  cmat* AA = (cmat*)A;
  int i;
  if (neig!=0) return 1;
  if (neig==0){
    for (i=0;i<n;i++){ v[i]=1.0; indx[i]=i;}
    *eig=AA->cnst; *nind=n;
  } else {  /* Or return an error */
    for (i=0;i<n;i++){ v[i]=0.0; }
    *eig=0; *nind=0;
  }
  return 0;
}


static int ConstMatRowNnz(void*A, int row, int nz[], int *nnz, int n){
  int i;
  for (i=0;i<n;i++){ nz[i]++; }
  *nnz=n;
  return 0;
}

static int ConstMatFNorm2(void*AA, int n, double *fnorm2){
  cmat* A = (cmat*)AA;
  *fnorm2=A->cnst*A->cnst*n*n;
  return 0;
}

static int ConstMatCountNonzeros(void*A, int *nnz, int n){
  *nnz=n*n;
  *nnz=1;
  *nnz=n;
  return 0;
}


static int ConstMatView(void* AA){
  cmat* A = (cmat*)AA;
  printf("Every element of the matrix is the same: %10.8e\n",A->cnst);
  return 0;
}

static int ConstMatTest(void* AA){
  return 0;
}


static int ConstMatDestroy(void* A){
  if (A) free(A);
  return 0;
}

static const char *datamatname="ALL ELEMENTS THE SAME";
static int ConstMatOpsInitialize(struct  DSDPDataMat_Ops* cmatops){
  int info;
  if (cmatops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(cmatops); DSDPCHKERR(info);
  cmatops->matfactor1=ConstMatFactor;
  cmatops->matgetrank=ConstMatGetRank;
  cmatops->matgeteig=ConstMatGetEig;
  cmatops->matvecvec=ConstMatVecVec;
  cmatops->matdot=ConstMatDot;
  cmatops->mataddrowmultiple=ConstMatAddRowMultiple;
  cmatops->mataddallmultiple=ConstMatAddMultiple;
  cmatops->matdestroy=ConstMatDestroy;
  cmatops->mattest=ConstMatTest;
  cmatops->matview=ConstMatView;
  cmatops->matrownz=ConstMatRowNnz;
  cmatops->matfnorm2=ConstMatFNorm2;
  cmatops->matnnz=ConstMatCountNonzeros;
  cmatops->id=14;
  cmatops->matname=datamatname;
  return 0;
}

