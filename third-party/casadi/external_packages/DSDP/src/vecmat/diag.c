
#include "dsdpschurmat_impl.h"
#include "dsdpdualmat_impl.h"
#include "dsdpdsmat_impl.h"
#include "dsdpsys.h"
/*!
\file diag.c
\brief DSDPDualMat, DSDPDSMat, and the DSDPSchurMat implentations for
diagonal matrices.
*/

typedef struct {
  int    n;
  double *val;
  int   owndata;
} diagmat;

static int DiagMatCreate(int,diagmat**);
static int DiagMatMult(void*,double[],double[],int);
static int DiagMatGetSize(void*, int *);
static int DiagMatAddRow2(void*, int, double, double[], int);
static int DiagMatDestroy(void*);
static int DiagMatView(void*);
static int DiagMatLogDeterminant(void*, double *);

/* static int DiagMatScale(double *, diagmat *); */

static int DiagMatCreate(int n,diagmat**M){
  int info;
  diagmat *M7;

  DSDPCALLOC1(&M7,diagmat,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&M7->val,double,n,&info);DSDPCHKERR(info);
  if (n>0 && M7->val==NULL) return 1;
  M7->owndata=1; M7->n=n;
  *M=M7;
  return 0;
}

static int DiagMatDestroy(void* AA){
  int info;
  diagmat* A=(diagmat*) AA;
  if (A->owndata && A->val){ DSDPFREE(&A->val,&info);DSDPCHKERR(info);}
  DSDPFREE(&A,&info);DSDPCHKERR(info);
  return 0;
}


static int DiagMatMult(void* AA, double x[], double y[], int n){

  diagmat* A=(diagmat*) AA;
  double *vv=A->val;
  int i;

  if (A->n != n) return 1;
  if (x==0 && n>0) return 3;
  if (y==0 && n>0) return 3;

  for (i=0; i<n; i++){
    y[i]=x[i]*vv[i];
  }  
  return 0;
}


static int DiagMatGetSize(void *AA, int *n){
  diagmat* A=(diagmat*) AA;
  *n=A->n;
  return 0;
}

static int DiagMatView(void* AA){
  diagmat* A=(diagmat*) AA;
  int i;
  for (i=0;i<A->n; i++){
    printf(" Row: %d, Column: %d, Value: %8.4e \n",i,i,A->val[i]);
  }
  return 0;
}

static int DiagMatAddRow2(void* AA, int nrow, double dd, double row[], int n){
  diagmat* A=(diagmat*) AA;
  A->val[nrow] += dd*row[nrow];
  return 0;
}


static int DiagMatAddElement(void*A, int nrow, double dd){
  diagmat* AA = (diagmat*)A;
  AA->val[nrow] += dd;
  return 0;
}

static int DiagMatZeros(void*A){
  diagmat* AA = (diagmat*)A;
  int n=AA->n;
  memset(AA->val,0,n*sizeof(double));
  return 0;
}

static int DiagMatSolve(void* A, double b[], double x[],int n){
  diagmat* AA = (diagmat*)A;
  double *v=AA->val;
  int i;
  for (i=0;i<n;i++){
    x[i]=b[i]/v[i];
  }
  return 0;
}

static int DiagMatSolve2(void* A, int indx[], int nindx, double b[], double x[],int n){
  diagmat* AA = (diagmat*)A;
  double *v=AA->val;
  int i,j;
  memset((void*)x,0,n*sizeof(double));
  for (j=0;j<nindx;j++){
    i=indx[j];
    x[i]=b[i]/v[i];
  }
  return 0;
}

static int DiagMatCholeskySolveBackward(void* A, double b[], double x[],int n){
  int i;
  for (i=0;i<n;i++){
    x[i]=b[i];
  }
  return 0;
}

static int DiagMatInvert(void *A){
  return 0;
}

static int DiagMatCholeskyFactor(void*A,int *flag){
  diagmat* AA = (diagmat*)A;
  double *v=AA->val;
  int i,n=AA->n;
  *flag=0;
  for (i=0;i<n;i++){
    if (v[i]<=0){ *flag=i+1; break;}
  }
  return 0;
}

static int DiagMatLogDeterminant(void*A, double *dd){
  diagmat* AA = (diagmat*)A;
  double d=0,*val=AA->val;
  int i,n=AA->n;
  for (i=0;i<n;i++){
    if (val[i]<=0) return 1;
    d+=log(val[i]);
  }
  *dd=d;
  return 0;
}

static int DiagMatTakeUREntriesP(void*A, double dd[], int nn, int n){
  diagmat* AA = (diagmat*)A;
  int i,ii;
  double *val=AA->val;
  for (i=0;i<n;i++){
    ii=(i+1)*(i+2)/2-1;
    val[i]=dd[ii];
  }
  return 0;
}
static int DiagMatTakeUREntriesU(void*A, double dd[], int nn, int n){
  diagmat* AA = (diagmat*)A;
  int i;
  double *val=AA->val;
  for (i=0;i<n;i++){
    val[i]=dd[i*n+i];
  }
  return 0;
}

static int DiagMatInverseAddP(void*A, double alpha, double dd[], int nn, int n){
  diagmat* AA = (diagmat*)A;
  int i,ii;
  double *val=AA->val;
  for (i=0;i<n;i++){
    ii=(i+1)*(i+2)/2-1;
    dd[ii]+=alpha/val[i];
  }
  return 0;
}
static int DiagMatInverseAddU(void*A, double alpha, double dd[], int nn, int n){
  diagmat* AA = (diagmat*)A;
  int i;
  double *val=AA->val;
  for (i=0;i<n;i++){
    dd[i*n+i]+=alpha/val[i];
  }
  return 0;
}

static int DiagMatFull(void*A, int* dfull){
  *dfull=1;
  return 0;
}

static struct  DSDPDualMat_Ops sdmatopsp;
static struct  DSDPDualMat_Ops sdmatopsu;
static const char* diagmatname="DIAGONAL";

static int DiagDualOpsInitializeP(struct  DSDPDualMat_Ops* sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDualMatOpsInitialize(sops); DSDPCHKERR(info);
  sops->matcholesky=DiagMatCholeskyFactor;
  sops->matsolveforward=DiagMatSolve;
  sops->matsolvebackward=DiagMatCholeskySolveBackward;
  sops->matinvert=DiagMatInvert;
  sops->matinverseadd=DiagMatInverseAddP;
  sops->matinversemultiply=DiagMatSolve2;
  sops->matseturmat=DiagMatTakeUREntriesP;
  sops->matfull=DiagMatFull;
  sops->matdestroy=DiagMatDestroy;
  sops->matgetsize=DiagMatGetSize;
  sops->matview=DiagMatView;
  sops->matlogdet=DiagMatLogDeterminant;
  sops->id=9;
  sops->matname=diagmatname;
  return 0;
}
static int DiagDualOpsInitializeU(struct  DSDPDualMat_Ops* sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDualMatOpsInitialize(sops); DSDPCHKERR(info);
  sops->matcholesky=DiagMatCholeskyFactor;
  sops->matsolveforward=DiagMatSolve;
  sops->matsolvebackward=DiagMatCholeskySolveBackward;
  sops->matinvert=DiagMatInvert;
  sops->matinversemultiply=DiagMatSolve2;
  sops->matseturmat=DiagMatTakeUREntriesU;
  sops->matfull=DiagMatFull;
  sops->matinverseadd=DiagMatInverseAddU;
  sops->matdestroy=DiagMatDestroy;
  sops->matgetsize=DiagMatGetSize;
  sops->matview=DiagMatView;
  sops->matlogdet=DiagMatLogDeterminant;
  sops->id=9;
  sops->matname=diagmatname;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDiagDualMatCreateP"
int DSDPDiagDualMatCreateP(int n, 
			   struct  DSDPDualMat_Ops **sops1, void**smat1,
			   struct  DSDPDualMat_Ops **sops2, void**smat2){
  diagmat *AA;
  int info;
  DSDPFunctionBegin;

  info=DiagMatCreate(n,&AA); DSDPCHKERR(info);
  info=DiagDualOpsInitializeP(&sdmatopsp); DSDPCHKERR(info);
  *sops1=&sdmatopsp;
  *smat1=(void*)AA;

  info=DiagMatCreate(n,&AA); DSDPCHKERR(info);
  info=DiagDualOpsInitializeP(&sdmatopsp); DSDPCHKERR(info);
  *sops2=&sdmatopsp;
  *smat2=(void*)AA;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDiagDualMatCreateU"
int DSDPDiagDualMatCreateU(int n, 
			   struct  DSDPDualMat_Ops **sops1, void**smat1,
			   struct  DSDPDualMat_Ops **sops2, void**smat2){
  diagmat *AA;
  int info;
  DSDPFunctionBegin;
  info=DiagMatCreate(n,&AA); DSDPCHKERR(info);
  info=DiagDualOpsInitializeU(&sdmatopsu); DSDPCHKERR(info);
  *sops1=&sdmatopsu;
  *smat1=(void*)AA;
  info=DiagMatCreate(n,&AA); DSDPCHKERR(info);
  info=DiagDualOpsInitializeU(&sdmatopsu); DSDPCHKERR(info);
  *sops2=&sdmatopsu;
  *smat2=(void*)AA;
  DSDPFunctionReturn(0);
}



static int DiagMatVecVec(void*A,double x[], int n, double *vv){
  diagmat* AA = (diagmat*)A;
  double *v=AA->val,vAv=0;
  int i;
  for (i=0;i<n;i++){
    vAv+=x[i]*x[i]*v[i];
  }
  *vv=vAv;
  return 0;
}

static int DDiagDSMatOpsInitP(struct  DSDPDSMat_Ops *ddiagops){
  int info;
  if (ddiagops==NULL) return 0;
  info=DSDPDSMatOpsInitialize(ddiagops);DSDPCHKERR(info);
  ddiagops->matseturmat=DiagMatTakeUREntriesP;
  ddiagops->matview=DiagMatView;
  ddiagops->matgetsize=DiagMatGetSize;
  ddiagops->matmult=DiagMatMult;
  ddiagops->matvecvec=DiagMatVecVec;
  ddiagops->matzeroentries=DiagMatZeros;
  ddiagops->matdestroy=DiagMatDestroy;
  ddiagops->id=9;
  ddiagops->matname=diagmatname;
  DSDPFunctionReturn(0);
}
static int DDiagDSMatOpsInitU(struct  DSDPDSMat_Ops *ddiagops){
  int info;
  if (ddiagops==NULL) return 0;
  info=DSDPDSMatOpsInitialize(ddiagops);DSDPCHKERR(info);
  ddiagops->matseturmat=DiagMatTakeUREntriesU;
  ddiagops->matview=DiagMatView;
  ddiagops->matgetsize=DiagMatGetSize;
  ddiagops->matmult=DiagMatMult;
  ddiagops->matvecvec=DiagMatVecVec;
  ddiagops->matzeroentries=DiagMatZeros;
  ddiagops->matdestroy=DiagMatDestroy;
  ddiagops->id=9;
  ddiagops->matname=diagmatname;
  DSDPFunctionReturn(0);
}

static struct  DSDPDSMat_Ops dsdiagmatopsp;
static struct  DSDPDSMat_Ops dsdiagmatopsu;

#undef __FUNCT__
#define __FUNCT__ "DSDPDiagDSMatP"
int DSDPCreateDiagDSMatP(int n,struct  DSDPDSMat_Ops* *dsmatops, void**dsmat){

  int info=0;
  diagmat *AA;

  DSDPFunctionBegin;
  info=DiagMatCreate(n,&AA); DSDPCHKERR(info);
  info=DDiagDSMatOpsInitP(&dsdiagmatopsp); DSDPCHKERR(info);
  *dsmatops=&dsdiagmatopsp;
  *dsmat=(void*)AA;
  DSDPFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DSDPDiagDSMatU"
int DSDPCreateDiagDSMatU(int n,struct  DSDPDSMat_Ops* *dsmatops, void**dsmat){

  int info=0;
  diagmat *AA;

  DSDPFunctionBegin;
  info=DiagMatCreate(n,&AA); DSDPCHKERR(info);
  info=DDiagDSMatOpsInitU(&dsdiagmatopsu); DSDPCHKERR(info);
  *dsmatops=&dsdiagmatopsu;
  *dsmat=(void*)AA;
  DSDPFunctionReturn(0);
}

static int DiagRowNonzeros(void*M, int row, double cols[], int *ncols,int nrows){
  DSDPFunctionBegin;
  *ncols = 1;
  cols[row]=1.0;
  DSDPFunctionReturn(0);
}

static int DiagAddDiag(void*M, double diag[], int m){
  diagmat* AA = (diagmat*)M;
  double *v=AA->val;
  int i;
  DSDPFunctionBegin;
  for (i=0;i<m;i++){
    v[i]+=diag[i];
  }
  DSDPFunctionReturn(0);
}

static int DiagMultiply(void*M, double xin[], double xout[], int m){
  diagmat* AA = (diagmat*)M;
  double *v=AA->val;
  int i;
  DSDPFunctionBegin;
  for (i=0;i<m;i++){
    xout[i]+=v[i]*xin[i];
  }
  DSDPFunctionReturn(0);
}

static int DiagShiftDiag(void*M, double dd){
  diagmat* AA = (diagmat*)M;
  double *v=AA->val;
  int i,m=AA->n;
  DSDPFunctionBegin;
  for (i=0;i<m;i++){
    v[i]+=dd;
  }
  DSDPFunctionReturn(0);
}

static int DiagAddElement(void*M, int ii, double dd){
  diagmat* AA = (diagmat*)M;
  DSDPFunctionBegin;
  AA->val[ii]+=dd;
  DSDPFunctionReturn(0);
}

static int DiagMatOnProcessor(void*A,int row,int*flag){
  *flag=1;
  return 0;
}

static int DiagAssemble(void*M){
  return 0;
}

static struct  DSDPSchurMat_Ops dsdpdiagschurops;

#undef __FUNCT__
#define __FUNCT__ "DSDPDiagSchurOps"
static int DiagSchurOps(struct  DSDPSchurMat_Ops *sops){
  int info;
  DSDPFunctionBegin;
  if (!sops) return 0;
  info=DSDPSchurMatOpsInitialize(sops); DSDPCHKERR(info);
  sops->matzero=DiagMatZeros;
  sops->mataddrow=DiagMatAddRow2;
  sops->mataddelement=DiagMatAddElement;
  sops->matdestroy=DiagMatDestroy;
  sops->matfactor=DiagMatCholeskyFactor;
  sops->matsolve=DiagMatSolve;
  sops->matadddiagonal=DiagAddDiag;
  sops->matshiftdiagonal=DiagShiftDiag;
  sops->mataddelement=DiagAddElement;
  sops->matscaledmultiply=DiagMultiply;
  sops->matassemble=DiagAssemble;
  sops->pmatonprocessor=DiagMatOnProcessor;
  sops->matrownonzeros=DiagRowNonzeros;
  sops->id=9;
  sops->matname=diagmatname;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPGetDiagSchurMat"
int DSDPGetDiagSchurMat(int n, struct  DSDPSchurMat_Ops **sops, void **data){
  int info=0;
  diagmat *AA;
  DSDPFunctionBegin;
  info=DiagMatCreate(n,&AA); DSDPCHKERR(info);
  info=DiagSchurOps(&dsdpdiagschurops); DSDPCHKERR(info);
  if (sops){*sops=&dsdpdiagschurops;}
  if (data){*data=(void*)AA;}
  DSDPFunctionReturn(0);
}
