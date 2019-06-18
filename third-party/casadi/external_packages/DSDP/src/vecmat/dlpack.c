#include "dsdpsys.h"
#include "dsdpvec.h"
#include "dsdplapack.h"

/*! \file dlpack.c
\brief DSDPDataMat, DSDPDualMat, DSDPDSMat, DSDPSchurMat, DSDPXMat,
objects implemented in dense upper packed symmetric format.
*/

typedef struct{
  char UPLO;
  double *val;
  double *v2;
  double *sscale;
  int  scaleit;
  int n;
  int owndata;
} dtpumat;

static const char* lapackname="DENSE,SYMMETRIC,PACKED STORAGE";

int DTPUMatEigs(void*AA,double W[],double IIWORK[], int nn1, double *mineig){
  dtpumat* AAA=(dtpumat*) AA;
  ffinteger info,INFO=0,M,N=AAA->n;
  ffinteger IL=1,IU=1,LDZ=1,IFAIL;
  ffinteger *IWORK=(ffinteger*)IIWORK;
  double *AP=AAA->val,ABSTOL=1e-13;
  double Z=0,VL=-1e10,VU=1;
  double *WORK;
  char UPLO=AAA->UPLO,JOBZ='N',RANGE='I';

  DSDPCALLOC2(&WORK,double,7*N,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&IWORK,ffinteger,5*N,&info);DSDPCHKERR(info);
  dspevx(&JOBZ,&RANGE,&UPLO,&N,AP,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,&Z,&LDZ,WORK,IWORK,&IFAIL,&INFO);

  /*
  DSDPCALLOC2(&WORK,double,2*N,&info);
  LWORK=2*N;
  dspevd(&JOBZ,&UPLO,&N,AP,W,&Z,&LDZ,WORK,&LWORK,IWORK,&LIWORK,&INFO);
  */
  *mineig=W[0];
  DSDPFREE(&WORK,&info);DSDPCHKERR(info);
  DSDPFREE(&IWORK,&info);DSDPCHKERR(info);
  return INFO;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACKROUTINE"
static void dtpuscalevec(double alpha, double v1[], double v2[], double v3[], int n){
  int i;
  for (i=0;i<n;i++){ 
    v3[i] = (alpha*v1[i]*v2[i]);
  }
}

static void dtpuscalemat(double vv[], double ss[], int n){
  int i;
  for (i=0;i<n;i++,vv+=i){ 
    dtpuscalevec(ss[i],vv,ss,vv,i+1);
  }
}

static int DTPUMatCreateWData(int n,double nz[], int nnz, dtpumat**M){
  int i,nn=(n*n+n)/2,info;
  double dtmp;
  dtpumat*M23;
  if (nnz<nn){DSDPSETERR1(2,"Array must have length of : %d \n",nn);}
  for (i=0;i<nnz;i++) dtmp=nz[i];
  DSDPCALLOC1(&M23,dtpumat,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&M23->sscale,double,n,&info);DSDPCHKERR(info);
  M23->owndata=0; M23->val=nz; M23->n=n; M23->UPLO='U';
  for (i=0;i<n;i++)M23->sscale[i]=1.0;
  M23->scaleit=0;
  *M=M23;
  return 0;
}



#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACK ROUTINE"


static int DTPUMatMult(void* AA, double x[], double y[], int n){
  dtpumat* A=(dtpumat*) AA;
  ffinteger ione=1,N=n;
  double BETA=0.0,ALPHA=1.0;
  double *AP=A->val,*Y=y,*X=x;
  char UPLO=A->UPLO;
  
  if (A->n != n) return 1;
  if (x==0 && n>0) return 3;
  dspmv(&UPLO,&N,&ALPHA,AP,X,&ione,&BETA,Y,&ione);
  return 0;
}

static int DTPUMatSolve(void* AA, double b[], double x[], int n){
  dtpumat* A=(dtpumat*) AA;
  ffinteger INFO,NRHS=1,LDB=A->n,N=A->n;
  double *AP=A->val,*ss=A->sscale;
  char UPLO=A->UPLO;

  if (N>0) LDB=N;  
  dtpuscalevec(1.0,ss,b,x,n);
  dpptrs(&UPLO, &N, &NRHS, AP, x, &LDB, &INFO );
  dtpuscalevec(1.0,ss,x,x,n);
  return INFO;
}

static int DTPUMatCholeskyFactor(void* AA, int *flag){
  dtpumat* A=(dtpumat*) AA;
  int i;
  ffinteger INFO,LDA=1,N=A->n;
  double *AP=A->val,*ss=A->sscale,*v;
  char UPLO=A->UPLO;

  if (N<=0) LDA=1;
  else LDA=N;
  if (A->scaleit){
    for (v=AP,i=0;i<N;i++){ ss[i]=*v;v+=(i+2);}
    for (i=0;i<N;i++){ ss[i]=1.0/sqrt(fabs(ss[i])+1.0e-8); }
    dtpuscalemat(AP,ss,N);
  }
  dpptrf(&UPLO, &N, AP, &INFO );
  *flag=INFO;
  return 0;
}

static int DTPUMatShiftDiagonal(void* AA, double shift){
  dtpumat* A=(dtpumat*) AA;
  int i,n=A->n;
  double *v=A->val;
  for (i=0; i<n; i++){
    *v+=shift;
    v+=i+2;
  }
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DTPUMatAssemble"
static int DTPUMatAssemble(void*M){
  int info;
  double shift=1.0e-15;
  DSDPFunctionBegin;
  info= DTPUMatShiftDiagonal(M, shift); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int DTPUMatRowNonzeros(void*M, int row, double cols[], int *ncols,int nrows){
  int i;
  DSDPFunctionBegin;
  *ncols = row+1;
  for (i=0;i<=row;i++){
    cols[i]=1.0;
  }
  for (i=row+1;i<nrows;i++){
    cols[i]=0.0;
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DTPUMatDiag"
static int DTPUMatDiag(void*M, int row, double dd){
  int ii;
  dtpumat*  ABA=(dtpumat*)M;
  DSDPFunctionBegin;
  ii=row*(row+1)/2+row;
  ABA->val[ii] +=dd;
  DSDPFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DTPUMatDiag2"
static int DTPUMatDiag2(void*M, double diag[], int m){
  int row,ii;
  dtpumat*  ABA=(dtpumat*)M;
  DSDPFunctionBegin;
  for (row=0;row<m;row++){
    ii=row*(row+1)/2+row;
    ABA->val[ii] +=diag[row];
  }
  DSDPFunctionReturn(0);
}

static int DTPUMatAddRow(void* AA, int nrow, double dd, double row[], int n){
  dtpumat* A=(dtpumat*) AA;
  ffinteger ione=1,nn,nnn;
  double *vv=A->val;

  nnn=nrow*(nrow+1)/2;
  nn=nrow+1;
  daxpy(&nn,&dd,row,&ione,vv+nnn,&ione);

  return 0;
}

static int DTPUMatZero(void* AA){
  dtpumat* A=(dtpumat*) AA;
  int mn=A->n*(A->n+1)/2;
  double *vv=A->val;
  memset((void*)vv,0,mn*sizeof(double));
  return 0;
}

static int DTPUMatGetSize(void *AA, int *n){
  dtpumat* A=(dtpumat*) AA;
  *n=A->n;
  return 0;
}

static int DTPUMatDestroy(void* AA){
  int info;
  dtpumat* A=(dtpumat*) AA;
  if (A && A->owndata){DSDPFREE(&A->val,&info);DSDPCHKERR(info);}
  if (A && A->sscale) {DSDPFREE(&A->sscale,&info);DSDPCHKERR(info);}
  if (A) {DSDPFREE(&A,&info);DSDPCHKERR(info);}
  return 0;
}

static int DTPUMatView(void* AA){
  dtpumat* M=(dtpumat*) AA;
  int i,j,kk=0;
  double *val=M->val;
  for (i=0; i<M->n; i++){
    for (j=0; j<=i; j++){
      printf(" %9.2e",val[kk]);
      kk++;
    }
    printf("\n");
  }
  return 0;
}


#include "dsdpschurmat_impl.h"
#include "dsdpsys.h"
static struct  DSDPSchurMat_Ops dsdpmmatops;

static int DSDPInitSchurOps(struct  DSDPSchurMat_Ops* mops){ 
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatOpsInitialize(mops);DSDPCHKERR(info);
  mops->matrownonzeros=DTPUMatRowNonzeros;
  mops->matscaledmultiply=DTPUMatMult;
  mops->mataddrow=DTPUMatAddRow;
  mops->mataddelement=DTPUMatDiag;
  mops->matadddiagonal=DTPUMatDiag2;
  mops->matshiftdiagonal=DTPUMatShiftDiagonal;
  mops->matassemble=DTPUMatAssemble;
  mops->matfactor=DTPUMatCholeskyFactor;
  mops->matsolve=DTPUMatSolve;
  mops->matdestroy=DTPUMatDestroy;
  mops->matzero=DTPUMatZero;
  mops->matview=DTPUMatView;
  mops->id=1;
  mops->matname=lapackname;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPGetLAPACKPUSchurOps"
int DSDPGetLAPACKPUSchurOps(int n,struct DSDPSchurMat_Ops** sops, void** mdata){
  int info,nn=n*(n+1)/2;
  double *vv;
  dtpumat*AA;
  DSDPFunctionBegin;
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DTPUMatCreateWData(n,vv,nn,&AA); DSDPCHKERR(info);
  AA->owndata=1;
  AA->scaleit=1;
  info=DSDPInitSchurOps(&dsdpmmatops); DSDPCHKERR(info);
  *sops=&dsdpmmatops;
  *mdata=(void*)AA;
  DSDPFunctionReturn(0);
}


static void daddrow(double *v, double alpha, int i, double row[], int n){
  double *s1;
  ffinteger j,nn=n,ione=1;
  nn=i+1; s1=v+i*(i+1)/2;
  daxpy(&nn,&alpha,s1,&ione,row,&ione);
  for (j=i+1;j<n;j++){
    s1+=j;
    row[j]+=alpha*s1[i];
  }
  return;
}

static int DTPUMatInverseMult(void* AA, int indx[], int nind, double x[], double y[], int n){
  dtpumat* A=(dtpumat*) AA;
  ffinteger ione=1,N=n;
  double BETA=0.0,ALPHA=1.0;
  double *AP=A->v2,*Y=y,*X=x;
  int i,ii;
  char UPLO=A->UPLO;
  
  if (A->n != n) return 1;
  if (x==0 && n>0) return 3;
  
  if (nind<n/4 ){
    memset((void*)y,0,n*sizeof(double));    
    for (ii=0;ii<nind;ii++){
      i=indx[ii];  ALPHA=x[i];
      daddrow(AP,ALPHA,i,y,n);
    }
  } else {
    ALPHA=1.0;
    dspmv(&UPLO,&N,&ALPHA,AP,X,&ione,&BETA,Y,&ione);
  }
  return 0;
}


static int DTPUMatCholeskyBackward(void* AA, double b[], double x[], int n){
  dtpumat* M=(dtpumat*) AA;
  ffinteger N=M->n,INCX=1;
  double *AP=M->val,*ss=M->sscale;
  char UPLO=M->UPLO,TRANS='N',DIAG='N';
  memcpy(x,b,N*sizeof(double));
  dtpsv(&UPLO,&TRANS,&DIAG, &N, AP, x, &INCX );
  dtpuscalevec(1.0,ss,x,x,n);
  return 0;
}


static int DTPUMatCholeskyForward(void* AA, double b[], double x[], int n){
  dtpumat* M=(dtpumat*) AA;
  ffinteger N=M->n,INCX=1;
  double *AP=M->val,*ss=M->sscale;
  char UPLO=M->UPLO,TRANS='T',DIAG='N';
  dtpuscalevec(1.0,ss,b,x,n);
  dtpsv(&UPLO,&TRANS,&DIAG, &N, AP, x, &INCX);
  return 0;
}

static int DenseSymPSDCholeskyForwardMultiply(void* AA, double x[], double y[], int n){
  dtpumat* A=(dtpumat*) AA;
  int i,j,k=0;
  ffinteger N=A->n;
  double *AP=A->val,*ss=A->sscale;
  
  if (x==0 && N>0) return 3;
  for (i=0;i<N;i++){
    for (j=0;j<=i;j++){
      y[i]+=AP[k]*x[j];
      k++;
    }
  }
  for (i=0;i<N;i++){ y[i]=y[i]/ss[i];}
  return 0;
}

static int DTPUMatLogDet(void* AA, double *dd){
  dtpumat* A=(dtpumat*) AA;
  int i,n=A->n;
  double d=0,*v=A->val,*ss=A->sscale;
  for (i=0; i<n; i++){
    if (*v<=0) return 1;
    d+=2*log(*v/ss[i]);
    v+=i+2;
  }
  *dd=d;
  return 0;
}

static int DTPUMatInvert(void* AA){
  dtpumat* A=(dtpumat*) AA;
  ffinteger INFO,N=A->n,nn=N*(N+1)/2;
  double *v=A->val,*AP=A->v2,*ss=A->sscale;
  char UPLO=A->UPLO;
  memcpy((void*)AP,(void*)v,nn*sizeof(double));
  dpptri(&UPLO, &N, AP, &INFO );
  if (INFO){
    INFO=DTPUMatShiftDiagonal(AA,1e-11);
    memcpy((void*)AP,(void*)v,nn*sizeof(double));
    dpptri(&UPLO, &N, AP, &INFO );
  }
  if (A->scaleit){
    dtpuscalemat(AP,ss,N);
  }
  return INFO;
}

static int DTPUMatInverseAdd(void* AA, double alpha, double y[], int nn, int n){
  dtpumat* A=(dtpumat*) AA;
  ffinteger N=nn,ione=1;
  double *v2=A->v2;
  daxpy(&N,&alpha,v2,&ione,y,&ione);
  return 0;
}


static int DTPUMatScaleDiagonal(void* AA, double dd){
  dtpumat* A=(dtpumat*) AA;
  int i,n=A->n;
  double *v=A->val;
  for (i=0; i<n; i++){
    *v*=dd;
    v+=i+2;    
  }
  return 0;
}

static int DTPUMatOuterProduct(void* AA, double alpha, double x[], int n){
  dtpumat* A=(dtpumat*) AA;
  ffinteger ione=1,N=n;
  double *v=A->val;
  char UPLO=A->UPLO;
  dspr(&UPLO,&N,&alpha,x,&ione,v);
  return 0;
}


static int DenseSymPSDNormF2(void* AA, int n, double *dddot){
  dtpumat* A=(dtpumat*) AA;
  ffinteger ione=1,nn=A->n*(A->n+1)/2;
  double dd,tt=sqrt(0.5),*val=A->val;
  int info;
  info=DTPUMatScaleDiagonal(AA,tt);
  dd=dnrm2(&nn,val,&ione);
  info=DTPUMatScaleDiagonal(AA,1.0/tt);
  *dddot=dd*dd*2;
  return 0;
}


/*
static int DTPUMatFNorm2(void* AA, double *mnorm){
  dtpumat* A=(dtpumat*) AA;
  ffinteger ione=1,n;
  double vv=0,*AP=A->val;
  n=A->n*(A->n+1)/2;
  vv=dnrm2(&n,AP,&ione);
  *mnorm=2.0*vv;
  return 1;
}
*/

#include "dsdpdualmat_impl.h"
#include "dsdpdatamat_impl.h"
#include "dsdpxmat_impl.h"
#include "dsdpdsmat_impl.h"



static int DTPUMatFull(void*A, int*full){
  *full=1;
  return 0;
}


static int DTPUMatGetDenseArray(void* A, double *v[], int*n){
  dtpumat*  ABA=(dtpumat*)A;
  *v=ABA->val;
  *n=(ABA->n)*(ABA->n+1)/2;
  return 0;
}

static int DTPUMatRestoreDenseArray(void* A, double *v[], int *n){
  *v=0;*n=0;
  return 0;
}

static int DDenseSetXMat(void*A, double v[], int nn, int n){
  double *vv;
  dtpumat*  ABA=(dtpumat*)A;
  vv=ABA->val;
  if (v!=vv){
    memcpy((void*)vv,(void*)v,nn*sizeof(double));  
  }
  return 0;
}

static int DDenseVecVec(void* A, double x[], int n, double *v){
  dtpumat*  ABA=(dtpumat*)A;
  int i,j,k=0;
  double dd=0,*val=ABA->val;
  *v=0.0;
  for (i=0; i<n; i++){
    for (j=0;j<i;j++){
      dd+=2*x[i]*x[j]*val[k];
      k++;
    }
    dd+=x[i]*x[i]*val[k];
    k++;
  }
  *v=dd;
  return 0;
}

static struct  DSDPDSMat_Ops tdsdensematops;
static int DSDPDSDenseInitializeOps(struct  DSDPDSMat_Ops* densematops){
  int info;
  if (!densematops) return 0;
  info=DSDPDSMatOpsInitialize(densematops); DSDPCHKERR(info);
  densematops->matseturmat=DDenseSetXMat;
  densematops->matview=DTPUMatView;
  densematops->matdestroy=DTPUMatDestroy;
  densematops->matgetsize=DTPUMatGetSize;
  densematops->matzeroentries=DTPUMatZero;
  densematops->matmult=DTPUMatMult;
  densematops->matvecvec=DDenseVecVec;
  densematops->id=1;
  densematops->matname=lapackname;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateDSMatWithArray"
int DSDPCreateDSMatWithArray(int n,double vv[], int nnz, struct  DSDPDSMat_Ops* *dsmatops, void**dsmat){
  int info;
  dtpumat*AA;
  DSDPFunctionBegin;
  info=DTPUMatCreateWData(n,vv,nnz,&AA); DSDPCHKERR(info);
  AA->owndata=0;
  info=DSDPDSDenseInitializeOps(&tdsdensematops); DSDPCHKERR(info);
  *dsmatops=&tdsdensematops;
  *dsmat=(void*)AA;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPCreateDSMat"
int DSDPCreateDSMat(int n,struct  DSDPDSMat_Ops* *dsmatops, void**dsmat){
  int info,nn=n*(n+1)/2;
  double *vv;
  dtpumat*AA;
  DSDPFunctionBegin;  
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DTPUMatCreateWData(n,vv,nn,&AA); DSDPCHKERR(info);
  info=DSDPDSDenseInitializeOps(&tdsdensematops); DSDPCHKERR(info);
  *dsmatops=&tdsdensematops;
  *dsmat=(void*)AA;
  AA->owndata=1;
  DSDPFunctionReturn(0);
}

static struct  DSDPVMat_Ops turdensematops;

static int DSDPDenseXInitializeOps(struct  DSDPVMat_Ops* densematops){
  int info;
  if (!densematops) return 0;
  info=DSDPVMatOpsInitialize(densematops); DSDPCHKERR(info);
  densematops->matview=DTPUMatView;
  densematops->matscalediagonal=DTPUMatScaleDiagonal;
  densematops->matshiftdiagonal=DTPUMatShiftDiagonal;
  densematops->mataddouterproduct=DTPUMatOuterProduct;
  densematops->matdestroy=DTPUMatDestroy;
  densematops->matfnorm2=DenseSymPSDNormF2;
  densematops->matgetsize=DTPUMatGetSize;
  densematops->matzeroentries=DTPUMatZero;
  densematops->matgeturarray=DTPUMatGetDenseArray;
  densematops->matrestoreurarray=DTPUMatRestoreDenseArray;
  densematops->matmineig=DTPUMatEigs;
  densematops->matmult=DTPUMatMult;
  densematops->id=1;
  densematops->matname=lapackname;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPXMatCreate"
int DSDPXMatPCreate(int n,struct  DSDPVMat_Ops* *xops, void * *xmat){
  int info,nn=n*(n+1)/2;
  double *vv;
  dtpumat*AA;
  DSDPFunctionBegin;
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DTPUMatCreateWData(n,vv,nn,&AA); DSDPCHKERR(info);
  AA->owndata=1;
  info=DSDPDenseXInitializeOps(&turdensematops); DSDPCHKERR(info);
  *xops=&turdensematops;
  *xmat=(void*)AA;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPXMatCreate"
int DSDPXMatPCreateWithData(int n,double nz[],int nnz,struct  DSDPVMat_Ops* *xops, void * *xmat){
  int i,info;
  double dtmp;
  dtpumat*AA;
  DSDPFunctionBegin;
  for (i=0;i<n*(n+1)/2;i++) dtmp=nz[i];
  info=DTPUMatCreateWData(n,nz,nnz,&AA); DSDPCHKERR(info);
  info=DSDPDenseXInitializeOps(&turdensematops); DSDPCHKERR(info);
  *xops=&turdensematops;
  *xmat=(void*)AA;
  DSDPFunctionReturn(0);
}


static struct  DSDPDualMat_Ops sdmatops;
static int SDualOpsInitialize(struct  DSDPDualMat_Ops* sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDualMatOpsInitialize(sops);DSDPCHKERR(info);
  sops->matseturmat=DDenseSetXMat;
  sops->matcholesky=DTPUMatCholeskyFactor;
  sops->matsolveforward=DTPUMatCholeskyForward;
  sops->matsolvebackward=DTPUMatCholeskyBackward;
  sops->matinvert=DTPUMatInvert;
  sops->matinverseadd=DTPUMatInverseAdd;
  sops->matinversemultiply=DTPUMatInverseMult;
  sops->matforwardmultiply=DenseSymPSDCholeskyForwardMultiply;
  sops->matfull=DTPUMatFull;
  sops->matdestroy=DTPUMatDestroy;
  sops->matgetsize=DTPUMatGetSize;
  sops->matview=DTPUMatView;
  sops->matlogdet=DTPUMatLogDet;
  sops->matname=lapackname;
  sops->id=1;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACKDualMatCreate"
int DSDPLAPACKPUDualMatCreate(int n,struct  DSDPDualMat_Ops **sops, void**smat){
  int info,nn=n*(n+1)/2;
  double *vv;
  dtpumat*AA;
  DSDPFunctionBegin;
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DTPUMatCreateWData(n,vv,nn,&AA); DSDPCHKERR(info);
  AA->owndata=1;;
  AA->scaleit=1;
  info=SDualOpsInitialize(&sdmatops);DSDPCHKERR(info);
  *sops=&sdmatops;
  *smat=(void*)AA;
  DSDPFunctionReturn(0);
}

static int switchptr(void* SD,void *SP){
  dtpumat *s1,*s2;
  s1=(dtpumat*)(SD);
  s2=(dtpumat*)(SP);
  s1->v2=s2->val;
  s2->v2=s1->val;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACKDualMatCreate2"
int DSDPLAPACKPUDualMatCreate2(int n,
			       struct  DSDPDualMat_Ops **sops1, void**smat1,
			       struct  DSDPDualMat_Ops **sops2, void**smat2){
  int info;
  DSDPFunctionBegin;
  info=DSDPLAPACKPUDualMatCreate(n,sops1,smat1);DSDPCHKERR(info);
  info=DSDPLAPACKPUDualMatCreate(n,sops2,smat2);DSDPCHKERR(info);
  info=switchptr(*smat1,*smat2);
  DSDPFunctionReturn(0);
}


typedef struct {
  int    neigs;
  double *eigval;
  double *an;
} Eigen;

typedef struct {
  dtpumat* AA;
  double alpha;
  Eigen    Eig;
} dvechmat;

#undef __FUNCT__  
#define __FUNCT__ "CreateDvechmatWData"
static int CreateDvechmatWdata(int n, double alpha, double vv[], dvechmat **A){
  int info,nn=(n*n+n)/2;
  dvechmat* V;
  DSDPCALLOC1(&V,dvechmat,&info);DSDPCHKERR(info);
  info=DTPUMatCreateWData(n,vv,nn,&V->AA); DSDPCHKERR(info);
  V->Eig.neigs=-1;
  V->Eig.eigval=0;
  V->Eig.an=0;
  V->alpha=alpha;
  *A=V;
  return 0;
}


static int DvechmatGetRowNnz(void* AA, int trow, int nz[], int *nnzz,int n){
  int k;
  *nnzz=n;
  for (k=0;k<n;k++) nz[k]++;
  return 0;
}

static int DTPUMatGetRowAdd(void* AA, int nrow, double ytmp, double row[], int n){
  dtpumat* A=(dtpumat*) AA;
  ffinteger i,k,nnn=n;
  double *v=A->val;

  nnn=nrow*(nrow+1)/2;
  for (i=0;i<nrow;i++){
    row[i]+=ytmp*v[nnn+i];
  }
  row[nrow]+=ytmp*v[nnn+nrow];
  for (i=nrow+1;i<n;i++){
    k=i*(i+1)/2+nrow;
    row[i]+=ytmp*v[k];
  }
  return 0;
}

static int DvechmatGetRowAdd(void* AA, int trow, double scl, double r[], int m){
  int info;
  dvechmat* A=(dvechmat*)AA;
  info=DTPUMatGetRowAdd((void*)A->AA ,trow,scl*A->alpha,r,m);
  return 0;
}

static int DvechmatAddMultiple(void* AA, double alpha, double r[], int nnn, int n){
  dvechmat* A=(dvechmat*)AA;
  ffinteger nn=nnn, ione=1;
  double *val=A->AA->val;
  alpha*=A->alpha;
  daxpy(&nn,&alpha,val,&ione,r,&ione);
  return 0;
}


static int DvechEigVecVec(void*, double[], int, double*);
static int DvechmatVecVec(void* AA, double x[], int n, double *v){
  dvechmat* A=(dvechmat*)AA;
  int i,j,k=0;
  double dd=0,*val=A->AA->val;
  *v=0.0;
  if (A->Eig.neigs<n/5){
    i=DvechEigVecVec(AA,x,n,&dd);
    *v=dd*A->alpha;
  } else {
    for (i=0; i<n; i++){
      for (j=0;j<i;j++){
	dd+=2*x[i]*x[j]*val[k];
	k++;
      }
      dd+=x[i]*x[i]*val[k];
      k++;
    }
    *v=dd*A->alpha;
  }
  return 0;
}


static int DvechmatFNorm2(void* AA, int n, double *v){
  dvechmat* A=(dvechmat*)AA;
  long int i,j,k=0;
  double dd=0,*x=A->AA->val;
  for (i=0; i<n; i++){
    for (j=0;j<i;j++){
      dd+=2*x[k]*x[k];
      k++;
    }
    dd+=x[k]*x[k];
    k++;
  }
  *v=dd*A->alpha*A->alpha;
  return 0;
}


static int DvechmatCountNonzeros(void* AA, int *nnz, int n){
  *nnz=n*(n+1)/2;
  return 0;
}


static int DvechmatDot(void* AA, double x[], int nn, int n, double *v){
  dvechmat* A=(dvechmat*)AA;
  ffinteger ione=1,nnn=nn;
  double dd,*val=A->AA->val;
  dd=ddot(&nnn,val,&ione,x,&ione);
  *v=2*dd*A->alpha;
  return 0;
}

/*
static int DvechmatNormF2(void* AA, int n, double *v){
  dvechmat* A=(dvechmat*)AA;
  return(DTPUMatNormF2((void*)(A->AA), n,v));
}
*/
#undef __FUNCT__  
#define __FUNCT__ "DvechmatDestroy"
static int DvechmatDestroy(void* AA){
  dvechmat* A=(dvechmat*)AA;
  int info;
  info=DTPUMatDestroy((void*)(A->AA));
  if (A->Eig.an){DSDPFREE(&A->Eig.an,&info);DSDPCHKERR(info);}
  if (A->Eig.eigval){DSDPFREE(&A->Eig.eigval,&info);DSDPCHKERR(info);}
  A->Eig.neigs=-1;
  DSDPFREE(&A,&info);DSDPCHKERR(info);
  return 0;
}


static int DvechmatView(void* AA){
  dvechmat* A=(dvechmat*)AA;
  dtpumat* M=A->AA;
  int i,j,kk=0;
  double *val=M->val;
  for (i=0; i<M->n; i++){
    for (j=0; j<=i; j++){
      printf(" %4.2e",A->alpha*val[kk]);
      kk++;
    }
    printf(" \n");
  }
  return 0;
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPCreateDvechmatEigs"
static int CreateEigenLocker(Eigen *E,int neigs, int n){
  int info;
  DSDPCALLOC2(&E->eigval,double,neigs,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&E->an,double,n*neigs,&info);DSDPCHKERR(info);
  E->neigs=neigs;
  return 0;
}


static int EigMatSetEig(Eigen* A,int row, double eigv, double v[], int n){
  double *an=A->an;
  A->eigval[row]=eigv;
  memcpy((void*)(an+n*row),(void*)v,n*sizeof(double));
  return 0;
}


static int EigMatGetEig(Eigen* A,int row, double *eigenvalue, double eigenvector[], int n){
  double* an=A->an;
  *eigenvalue=A->eigval[row];
  memcpy((void*)eigenvector,(void*)(an+n*row),n*sizeof(double));
  return 0;
}


static int DvechmatComputeEigs(dvechmat*,double[],int,double[],int,double[],int,int[],int);

static int DvechmatFactor(void*AA, double dmatp[], int nn0, double dwork[], int n, double ddwork[], int n1, int iptr[], int n2){

  int info;
  dvechmat*  A=(dvechmat*)AA;
  if (A->Eig.neigs>=0) return 0;
  info=DvechmatComputeEigs(A,dmatp,nn0,dwork,n,ddwork,n1,iptr,n2);DSDPCHKERR(info);
  return 0;
}

static int DvechmatGetRank(void *AA,int *rank, int n){
  dvechmat*  A=(dvechmat*)AA;
  if (A->Eig.neigs>=0){
    *rank=A->Eig.neigs;
  } else {
    DSDPSETERR(1,"Vech Matrix not factored yet\n");
  }
  return 0;
}

static int DvechmatGetEig(void* AA, int rank, double *eigenvalue, double vv[], int n, int indz[], int *nind){
  dvechmat*  A=(dvechmat*)AA;
  int i,info;
  double dd;
  if (A->Eig.neigs>0){
    info=EigMatGetEig(&A->Eig,rank,&dd,vv,n);DSDPCHKERR(info);
    *nind=n;
    *eigenvalue=dd*A->alpha;
    for (i=0;i<n;i++){ indz[i]=i;}
  } else {
    DSDPSETERR(1,"Vech Matrix not factored yet\n");
  }
  return 0;  
}

static int DvechEigVecVec(void* AA, double v[], int n, double *vv){
  dvechmat*  A=(dvechmat*)AA;
  int i,rank,neigs;
  double *an,dd,ddd=0,*eigval;
  if (A->Eig.neigs>=0){
    an=A->Eig.an;
    neigs=A->Eig.neigs;
    eigval=A->Eig.eigval;
    for (rank=0;rank<neigs;rank++){
      for (dd=0,i=0;i<n;i++){
	dd+=v[i]*an[i];
      }
      an+=n;
      ddd+=dd*dd*eigval[rank];
    }
    *vv=ddd*A->alpha;
  } else {
    DSDPSETERR(1,"Vech Matrix not factored yet\n");
  }
  return 0;  
}


static struct  DSDPDataMat_Ops dvechmatops;
static const char *datamatname="DENSE VECH MATRIX";

static int DvechmatOpsInitialize(struct  DSDPDataMat_Ops *sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(sops); DSDPCHKERR(info);
  sops->matvecvec=DvechmatVecVec;
  sops->matdot=DvechmatDot;
  sops->mataddrowmultiple=DvechmatGetRowAdd;
  sops->mataddallmultiple=DvechmatAddMultiple;
  sops->matview=DvechmatView;
  sops->matdestroy=DvechmatDestroy;
  sops->matfactor2=DvechmatFactor;
  sops->matgetrank=DvechmatGetRank;
  sops->matgeteig=DvechmatGetEig;
  sops->matrownz=DvechmatGetRowNnz;
  sops->matfnorm2=DvechmatFNorm2;
  sops->matnnz=DvechmatCountNonzeros;
  sops->id=1;
  sops->matname=datamatname;
  return 0;
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDmat"
int DSDPGetDMat(int n,double alpha, double *val, struct  DSDPDataMat_Ops**sops, void**smat){ 
  int info,k;
  double dtmp;
  dvechmat* A;
  DSDPFunctionBegin;

  for (k=0;k<n*(n+1)/2;++k) dtmp=val[k];
  info=CreateDvechmatWdata(n,alpha,val,&A); DSDPCHKERR(info);
  A->Eig.neigs=-1;
  info=DvechmatOpsInitialize(&dvechmatops); DSDPCHKERR(info);
  if (sops){*sops=&dvechmatops;}
  if (smat){*smat=(void*)A;}
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DvechmatComputeEigs"
static int DvechmatComputeEigs(dvechmat* AA,double DD[], int nn0, double W[], int n, double WORK[], int n1, int iiptr[], int n2){

  int i,j,k,neigs,info;
  long int *i2darray=(long int*)DD;
  int ownarray1=0,ownarray2=0,ownarray3=0;
  double *val=AA->AA->val;
  double *dmatarray=0,*dworkarray=0,eps=1.0e-12;
  int nn1=0,nn2=0;
  
  /* create a dense array in which to put numbers */
  if (n*n>nn1){
    DSDPCALLOC2(&dmatarray,double,(n*n),&info); DSDPCHKERR(info);
    ownarray1=1;
  }
  memset((void*)dmatarray,0,n*n*sizeof(double));

  if (n*n>nn2){
    DSDPCALLOC2(&dworkarray,double,(n*n),&info); DSDPCHKERR(info);
    ownarray2=1;
  }

  if (n*n*sizeof(long int)>nn0*sizeof(double)){
    DSDPCALLOC2(&i2darray,long int,(n*n),&info); DSDPCHKERR(info);
    ownarray3=1;
  }


  for (k=0,i=0; i<n; i++){
    for (j=0; j<=i; j++){
      dmatarray[i*n+j] += val[k];
      if (i!=j){
	dmatarray[j*n+i] += val[k];      
      }
      k++;
    }
  }
  /* Call LAPACK to compute the eigenvalues */
  info=DSDPGetEigs(dmatarray,n,dworkarray,n*n,i2darray,n*n,
		   W,n,WORK,n1,iiptr+3*n,n2-3*n); DSDPCHKERR(info);
  
  /* Count the nonzero eigenvalues */
  for (neigs=0,i=0;i<n;i++){
    if (fabs(W[i])> eps ){ neigs++;}
  }

  info=CreateEigenLocker(&AA->Eig,neigs,n);DSDPCHKERR(info);
  
  /* Copy into structure */
  for (neigs=0,i=0; i<n; i++){
    if (fabs(W[i]) > eps){
      info=EigMatSetEig(&AA->Eig,neigs,W[i],dmatarray+n*i,n);DSDPCHKERR(info);
      neigs++;
    }
  }
  
  if (ownarray1){ DSDPFREE(&dmatarray,&info);DSDPCHKERR(info);}
  if (ownarray2){ DSDPFREE(&dworkarray,&info);DSDPCHKERR(info);}
  if (ownarray3){ DSDPFREE(&i2darray,&info);DSDPCHKERR(info);}
  return 0;
}

