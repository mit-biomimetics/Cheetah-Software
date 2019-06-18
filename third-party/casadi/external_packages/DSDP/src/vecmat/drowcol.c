
typedef struct {
  int    rc;
  const double *val;
  int    n;
  double x,y;
} rcmat;

#include "dsdpdatamat_impl.h"
#include "dsdpsys.h"
/*!
\file drowcol.c
\brief DSDPDataMat object such that A(i,j) is nonzero only if
i or j equals some integer k. Not completed.
*/

static int RCMatDestroy(void*);
static int RCMatView(void*);
static int RCMatVecVec(void*, double[], int, double *);
static int RCMatDot(void*, double[], int, int, double *);
static int RCMatGetRank(void*, int*, int);
static int RCMatFactor(void*);
static int RCMatGetEig(void*, int, double*, double[], int, int[], int*);
static int RCMatAddRowMultiple(void*, int, double, double[], int);
static int RCMatAddMultiple(void*, double, double[], int, int);
static int RCMatGetRowNnz(void*, int, int[], int*, int);

static int RCMatCreate(int n,int rowcol, const double v[], rcmat**M){
  rcmat *M16;
  M16=(rcmat*) malloc(1*sizeof(rcmat));
  M16->val=v;
  M16->rc=rowcol;
  M16->n=n;
  *M=M16;
  return 0;
}

static int RCMatDestroy(void* AA){
  free(AA);
  return 0;
}
static int RCMatVecVec(void* A, double x[], int n, double *v){
  rcmat*RC=(rcmat*)A;
  int i;
  double vv=0;
  const double *val=RC->val;
  for (i=0;i<n;i++){ vv+=val[i]*x[i];}
  *v=2*vv*x[RC->rc];
  return 0;
}
static int RCMatDot(void* A, double x[], int nn, int n1, double *v){
  rcmat*RC=(rcmat*)A;
  int i,k=0,rc=RC->rc,n=RC->n;
  double vv=0;
  const double *val=RC->val;
  k=rc*(rc+1)/2;
  for (i=0;i<=rc;i++,k++){ vv+=v[k]*val[i]; }
  for (i=rc+1;i<n;i++,k+=i){ vv+=v[k+rc]*val[i];}
  *v=2*vv;
  return 0;
}
static int RCMatView(void* A){
  rcmat*RC=(rcmat*)A;
  int i;
  printf("Row Col %d\n",RC->rc);
  for (i=0;i<RC->n;i++){
    printf("%4.4e ",RC->val[i]);
  }
  printf("\n");
  return 0;
}

static int RCMatFactor(void* A){
  /* Must compute x,y such that eigs are orthogonal */
  rcmat*RC=(rcmat*)A;
  int i,rc=RC->rc;
  double vnorm2=0;
  const double *val=RC->val;
  for (i=0;i<RC->n;i++){ vnorm2+=val[i]*val[i];} 
  vnorm2=sqrt(vnorm2);
  if (val[rc]>=0){
    RC->y=-sqrt(vnorm2);
    RC->x=sqrt(2*val[rc]+RC->y*RC->y);
  } else {
    RC->x=sqrt(vnorm2);
    RC->y=-sqrt(-2*val[rc]+RC->x*RC->x);
  }
  return 0;
}
static int RCMatGetRank(void *A, int*rank, int n){
  *rank=2;
  return 0;
}

static int RCMatGetEig(void*A, int neig, double *eig, double v[], int n,int spind[], int *nind){
  rcmat*RC=(rcmat*)A;
  int i,rc=RC->rc;
  double x=RC->x,y=RC->y,xmy=x-y;
  const double *val=RC->val;
  if (neig==0){
    for (i=0;i<n;i++){spind[i]=i;}
    for (i=0;i<n;i++){v[i]=val[i]/xmy;}
    v[rc]=RC->x;
    *nind=n;
    *eig=1.0;
  } else if (neig==1){
    for (i=0;i<n;i++){spind[i]=i;}
    for (i=0;i<n;i++){v[i]=val[i]/xmy;}
    v[rc]=RC->y;
    *nind=n;
    *eig=-1.0;
  } else { *eig=0;*nind=0;}
  return 0;
}

static int RCMatGetRowNnz(void*A, int nrow, int nz[], int *nnzz, int n){
  rcmat*RC=(rcmat*)A;
  int i;
  *nnzz=1;
  if (nrow==RC->rc){for (i=0;i<n;i++){nz[i]++;}*nnzz=n;}
  nz[nrow]++;
  return 0;
}

static int RCMatAddRowMultiple(void*A, int nrow, double dd, double rrow[], int n){
  rcmat*RC=(rcmat*)A;
  int i;
  if (nrow==RC->rc){ 
    for (i=0;i<n;i++){rrow[i]+=dd*RC->val[i];}
  }
  rrow[nrow]+=dd*RC->val[nrow];
  return 0;
}
static int RCMatAddMultiple(void*A, double dd, double vv[], int nn, int n1){
  rcmat*RC=(rcmat*)A;
  int i,rc=RC->rc,k=rc*(rc+1)/2,n=RC->n;
  const double *val=RC->val;
  for (i=0;i<=rc;i++,k++){
    vv[k]+=dd*val[i];
  }
  for (i=rc+1;i<n;i++,k+=i){
    vv[k+rc]+=dd*val[i];
  }
  return 0;
}
static int RCMatFNorm(void*A, int n, double *fnorm){
  rcmat*RC=(rcmat*)A;
  int i,rc=RC->rc;
  double ff=0;
  const double *val=RC->val;
  for (i=0;i<n;i++){
    ff+=val[i]*val[i];
  }
  ff*=2;
  ff+=2*val[rc]*val[rc];
  *fnorm=ff;
  return 0;
}
static int RCMatCountNonzeros(void*A, int *nnz, int n){
  *nnz=2*n-1;
  return 0;
}
static struct  DSDPDataMat_Ops rcmatops;
static const char *datamatname="One Row and Column matrix";
static int RCMatOperationsInitialize(struct  DSDPDataMat_Ops* rcmatoperator){
  int info;
  if (rcmatoperator==NULL) return 0;
  info=DSDPDataMatOpsInitialize(rcmatoperator); if (info){ return info;}
  rcmatoperator->matfactor1=RCMatFactor;
  rcmatoperator->matgetrank=RCMatGetRank;
  rcmatoperator->matgeteig=RCMatGetEig;
  rcmatoperator->matvecvec=RCMatVecVec;
  rcmatoperator->matrownz=RCMatGetRowNnz;
  rcmatoperator->matdot=RCMatDot;
  rcmatoperator->matfnorm2=RCMatFNorm;
  rcmatoperator->matnnz=RCMatCountNonzeros;
  rcmatoperator->mataddrowmultiple=RCMatAddRowMultiple;
  rcmatoperator->mataddallmultiple=RCMatAddMultiple;
  rcmatoperator->matdestroy=RCMatDestroy;
  rcmatoperator->matview=RCMatView;
  rcmatoperator->matname=datamatname;
  rcmatoperator->id=27;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPGetRCMat"
int DSDPGetRCMat(int n,const double val[],int rowcol, struct  DSDPDataMat_Ops**sops, void**smat){
  rcmat *AA;
  int info;
  DSDPFunctionBegin;
  info=RCMatCreate(n,rowcol,val,&AA);DSDPCHKERR(info);
  info=RCMatOperationsInitialize(&rcmatops); DSDPCHKERR(info);
  if (sops){*sops=&rcmatops;}
  if (smat){*smat=(void*)AA;}
  DSDPFunctionReturn(0);
}

