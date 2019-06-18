#include "dsdpdatamat_impl.h"
#include "dsdpsys.h"
/*!
\file zeromat.c
\brief DSDPDataMat object with all zero entries.
*/
static int ZDestroy(void*);
static int ZView(void*);
static int ZVecVec(void*, double[], int, double *);
static int ZDot(void*, double[], int, int,double *);
static int ZGetRank(void*, int*,int);
static int ZFactor(void*);
static int ZGetEig(void*, int, double*, double[], int,int[],int*);
static int ZAddRowMultiple(void*, int, double, double[], int);
static int ZAddMultiple(void*, double, double[], int,int);
static int ZRowNnz(void*, int, int[], int*, int);

static struct  DSDPDataMat_Ops zeromatops;
static int ZeroMatopsInitialize(struct DSDPDataMat_Ops*);

int DSDPGetZeroDataMatOps(struct DSDPDataMat_Ops** zops){
  int info;
  info=ZeroMatopsInitialize(&zeromatops); if (info){return info;}
  if (zops){*zops=&zeromatops;}
  return info;
}

static int ZFactor(void *A){
  return 0;
}

static int ZGetRank(void*A,int*rank,int n){
  *rank=0;
  return 0;
}

static int ZGetEig(void*A,int neig, double *eig, double v[], int n,int indx[],int*nind){
  *eig=0.0;
  *nind=0;
  return 0;
}

static int ZDot(void*A, double x[], int nn, int n, double *sum){
  *sum=0.0;
  return 0;
}

static int ZVecVec(void*A, double x[], int n, double *sum){
  *sum=0.0;
  return 0;
}

static int ZAddMultiple(void*A, double dd, double row[], int nn, int n){
  return 0;
}

static int ZAddRowMultiple(void*A, int nrow, double dd, double row[], int n){
  return 0;
}

static int ZRowNnz(void*A, int row, int nz[], int *nnz, int n){
  *nnz=0;
  return 0;
}

static int ZDestroy(void*A){
  return 0;
}

static int ZNorm2(void*A,int n,double *v){
  *v=0;
  return 0;
}

static int ZView(void*A){
  printf("All zeros\n");
  return 0;
}

static const char* datamatname="MATRIX OF ZEROS";

static int ZeroMatopsInitialize(struct  DSDPDataMat_Ops* sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(sops); if (info){ return info;}
  sops->matfactor1=ZFactor;
  sops->matgetrank=ZGetRank;
  sops->matgeteig=ZGetEig;
  sops->matvecvec=ZVecVec;
  sops->matdot=ZDot;
  sops->matfnorm2=ZNorm2;
  sops->matrownz=ZRowNnz;
  sops->mataddrowmultiple=ZAddRowMultiple;
  sops->mataddallmultiple=ZAddMultiple;
  sops->matdestroy=ZDestroy;
  sops->matview=ZView;
  sops->id=10; 
  sops->matname=datamatname;
  return 0;
}



