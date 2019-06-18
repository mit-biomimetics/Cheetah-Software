#include "dsdpsys.h"
#include "dsdpvec.h"
#include "dsdplapack.h"
#include "dsdpdatamat_impl.h"

/*! \file dufull.c
\brief DSDPDataMat, DSDPDualMat, DSDPDSMat, DSDPSchurMat, DSDPXMat
objects implemented in symmetric upper full symmetric format.
*/

typedef enum {
  Init=0,
  Assemble=1,
  Factored=2,   /* fail to allocate required space */
  Inverted=3,    /* indefinity is detected          */
  ISymmetric=4
} MatStatus;

typedef struct{
  char UPLO;
  int LDA;
  double *val,*v2;
  double *sscale;
  double *workn;
  int  scaleit;
  int n;
  int owndata;
  MatStatus status;
} dtrumat;

static int DTRUMatView(void*);


#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACKROUTINE"
static void dtruscalevec(double alpha, double v1[], double v2[], double v3[], int n){
  int i;
  for (i=0;i<n;i++){ 
    v3[i] = (alpha*v1[i]*v2[i]);
  }
  return;
}

static void dsydadd(double x[], double s[], double y[],int n){
  int i;
  for (i=0;i<n;i++){ 
    y[i] += x[i]*(1/(s[i]*s[i])-2);
  }
  return;
}


static void dtruscalemat(double vv[], double ss[], int n,int LDA){
  int i;
  for (i=0;i<n;i++,vv+=LDA){ 
    dtruscalevec(ss[i],vv,ss,vv,i+1);
  }
  return;
}

static void DTRUMatOwnData(void* A, int owndata){
  dtrumat* AA=(dtrumat*)A;
  AA->owndata=owndata;
  return;
}

static int SUMatGetLDA(int n){
  int nlda=n;
  if (n>8 && nlda%2==1){ nlda++;}
  if (n>100){
    while (nlda%8!=0){ nlda++;}
  }
  /*
  printf("LDA: %d %d %d \n",n,nlda,(int)sizeof(double));
  */
  return (nlda);
}

static int DTRUMatCreateWData(int n,int LDA,double nz[], int nnz, dtrumat**M){
  int i,info;
  dtrumat*M23;
  if (nnz<n*n){DSDPSETERR1(2,"Array must have length of : %d \n",n*n);}
  DSDPCALLOC1(&M23,dtrumat,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&M23->sscale,double,n,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&M23->workn,double,n,&info);DSDPCHKERR(info);
  M23->owndata=0; M23->val=nz; M23->n=n; M23->UPLO='U';M23->LDA=n;
  M23->status=Init;
  for (i=0;i<n;i++)M23->sscale[i]=1.0;
  M23->scaleit=1;
  M23->LDA=LDA;
  if (n<=0){M23->LDA=1;}
  *M=M23;
  return 0;
}



#undef __FUNCT__
#define __FUNCT__ "DSDPGetEigs"
int DSDPGetEigs(double A[],int n, double AA[], int nn0, long int IA[], int nn1, 
		double W[],int n2,
		double WORK[],int nd, int LLWORK[], int ni){
  ffinteger N=n,LDA=n,LDZ=n,LWORK=nd,INFO=0;
  char UPLO='U',JOBZ='V',RANGE='A';
  /* Faster, but returns more error codes. ie. INFO>0 sometimes*/

  LDA=DSDPMax(1,n); 
  LDZ=DSDPMax(1,n); 
  if ( n2/2.5 > n || (ni<10*n+1) || (nd<26*n+1) || (nn0 < n*LDA) || (nn1<LDZ*n) ){
    /*
    printf("n: %d, ni: %d, nd: %d\n",n,ni/n,nd/n);
    printf("SLOW VERSION\n");
    */
    dsyev(&JOBZ,&UPLO,&N,A,&LDA,W,WORK,&LWORK,&INFO);
  } else {

    int i;
    ffinteger M,IL=1,IU=n,*ISUPPZ=(ffinteger*)IA;
    ffinteger *IWORK=(ffinteger*)LLWORK,LIWORK=(ffinteger)ni;
    double *Z=AA,VL=-1e10,VU=1e10,ABSTOL=0;
    /*   ABSTOL=dlamch_("Safe minimum" ); */
    if (0==1){
      dtrumat*M;
      DTRUMatCreateWData(n,n,A,n*n,&M);
      DTRUMatView((void*)M);
    }
    /*
      printf("N: %d N2: %d , ",n,n2);
    */
    /*
    LWORK=26*n; LIWORK=10*n;
    */
    /*
    printf("JOBZ: %c, RANGE: %c, UPLO: %c, N: %d LDA: %d \n",
	   JOBZ,RANGE,UPLO, N,LDA);
    printf("VL: %4.4e, VU: %4.4e, IL: %d, IU: %d, ABSTOL: %4.4e, LDZ: %d\n",
	   VL,VU,IL,IU,ABSTOL,LDZ);
    printf("LWORK: %d, LIWORK: %d\n",LWORK,LIWORK);
    */

      dsyevr(&JOBZ,&RANGE,&UPLO,&N,A,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z,&LDZ,ISUPPZ,WORK,&LWORK,IWORK,&LIWORK,&INFO);
    for (i=0;i<N*N;i++){A[i]=Z[i];}    

  }
  return INFO;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPGetEigsSTEP"
int DSDPGetEigsSTEP(double A[],int n, double AA[], int nn0, long int IA[], int nn1, 
		    double W[],int n2,
		    double WORK[],int nd, int LLWORK[], int ni){
  int info;
  info=DSDPGetEigs(A,n,AA,nn0,IA,nn1,W,n2,WORK,nd,LLWORK,ni);
  return info;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPGetEigs2"
int DSDPGetEigs2(double A[],int n, double AA[], int nn0, long int IA[], int nn1, 
		double W[],int n2,
		double WORK[],int nd, int LLWORK[], int ni){
  ffinteger N=n,LDA=n,LDZ=n,LWORK=nd,INFO=0;
  char UPLO='U',JOBZ='V';
  /* Works and uses less memory, but slower by a factor of about 2 or 3 */
  LDA=DSDPMax(1,n); 
  LDZ=DSDPMax(1,n); 
  if (0==1){
      dtrumat*M;
      DTRUMatCreateWData(n,n,A,n*n,&M);
      DTRUMatView((void*)M);
  }
  dsyev(&JOBZ,&UPLO,&N,A,&LDA,W,WORK,&LWORK,&INFO);
  return INFO;
}


#undef __FUNCT__
#define __FUNCT__ "DSDP FULL SYMMETRIC U LAPACK ROUTINE"

static int DTRUMatMult(void* AA, double x[], double y[], int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger N=n,LDA=A->LDA,INCX=1,INCY=1;
  double BETA=0.0,ALPHA=1.0;
  double *AP=A->val,*Y=y,*X=x;
  char UPLO=A->UPLO,TRANS='N';
  
  if (A->n != n) return 1;
  if (x==0 && n>0) return 3;
  if (0==1){
    dgemv(&TRANS,&N,&N,&ALPHA,AP,&LDA,X,&INCX,&BETA,Y,&INCY);
  } else {
    dsymv(&UPLO,&N,&ALPHA,AP,&LDA,X,&INCX,&BETA,Y,&INCY);
  }
  return 0;
}


static int DTRUMatMultR(void* AA, double x[], double y[], int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger N=n,LDA=A->LDA,INCX=1,INCY=1;
  double ALPHA=1.0;
  double *AP=A->val,*Y=y,*X=x,*ss=A->sscale,*W=A->workn;
  char UPLO=A->UPLO,TRANS='N',DIAG='U';
  
  UPLO='L';
  if (A->n != n) return 1;
  if (x==0 && n>0) return 3;
  /*  dsymv(&UPLO,&N,&ALPHA,AP,&LDA,X,&INCX,&BETA,Y,&INCY); */
  
  memset(y,0,n*sizeof(double));
  
  memcpy(W,X,n*sizeof(double));
  TRANS='N'; UPLO='L';
  dtrmv(&UPLO,&TRANS,&DIAG,&N,AP,&LDA,W,&INCY);
  daxpy(&N,&ALPHA,W,&INCX,Y,&INCY);
  
  memcpy(W,X,n*sizeof(double));
  TRANS='T'; UPLO='L';
  dtrmv(&UPLO,&TRANS,&DIAG,&N,AP,&LDA,W,&INCY);
  daxpy(&N,&ALPHA,W,&INCX,Y,&INCY);
  
  dsydadd(x,ss,y,n);
  return 0;
}


static void DTRUMatScale(void*AA){
  dtrumat* A=(dtrumat*) AA;
  int i,N=A->n,LDA=A->LDA;
  double *ss=A->sscale,*v=A->val;
  for (i=0;i<N;i++){ ss[i]=*v;v+=(LDA+1);}
  for (i=0;i<N;i++){ if (ss[i]!=0.0){ss[i]=1.0/sqrt(fabs(ss[i]));} else {ss[i]=1.0; }}
  dtruscalemat(A->val,ss,N,LDA);
}

static int DTRUMatCholeskyFactor(void* AA, int *flag){
  dtrumat* A=(dtrumat*) AA;
  ffinteger INFO,LDA=A->LDA,N=A->n;
  double *AP=A->val;
  char UPLO=A->UPLO;

  if (A->scaleit){ DTRUMatScale(AA);}
  dpotrf(&UPLO, &N, AP, &LDA, &INFO );
  *flag=INFO;
  A->status=Factored;
  return 0;
}


int dtrsm2(char*,char*,char*,char*,ffinteger*,ffinteger*,double*,double*,ffinteger*,double*,ffinteger*); 
 
static int DTRUMatSolve(void* AA, double b[], double x[],int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger ierr,INFO=0,NRHS=1,LDA=A->LDA,LDB=A->LDA,N=A->n;
  double *AP=A->val,*ss=A->sscale,ONE=1.0;
  char SIDE='L',UPLO=A->UPLO,TRANSA='N',DIAG='N';

  dtruscalevec(1.0,ss,b,x,n);

  if (0==1){
    dpotrs(&UPLO, &N, &NRHS, AP, &LDA, x, &LDB, &INFO );
  } else {
    TRANSA='T';
    ierr=dtrsm2(&SIDE,&UPLO,&TRANSA,&DIAG,&N, &NRHS, &ONE, AP, &LDA, x, &LDB );
    TRANSA='N';
    ierr=dtrsm2(&SIDE,&UPLO,&TRANSA,&DIAG,&N, &NRHS, &ONE, AP, &LDA, x, &LDB );
  }

  dtruscalevec(1.0,ss,x,x,n);
  return INFO;
}


static int DTRUMatShiftDiagonal(void* AA, double shift){
  dtrumat* A=(dtrumat*) AA;
  int i,n=A->n, LDA=A->LDA;
  double *v=A->val;
  for (i=0; i<n; i++){
    v[i*LDA+i]+=shift;
  }
  return 0;
}

static int DTRUMatAddRow(void* AA, int nrow, double dd, double row[], int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger ione=1,LDA=A->LDA,nn,INCX=1,INCY=A->LDA;
  double *vv=A->val;

  nn=nrow;
  daxpy(&nn,&dd,row,&INCX,vv+nrow,&INCY);
  nn=nrow+1;
  daxpy(&nn,&dd,row,&ione,vv+nrow*LDA,&ione);

  return 0;
}

static int DTRUMatZero(void* AA){
  dtrumat* A=(dtrumat*) AA;
  int mn=A->n*(A->LDA);
  double *vv=A->val;
  memset((void*)vv,0,mn*sizeof(double));
  A->status=Assemble;
  return 0;
}

static int DTRUMatGetSize(void *AA, int *n){
  dtrumat* A=(dtrumat*) AA;
  *n=A->n;
  return 0;
}

static int DTRUMatDestroy(void* AA){
  int info;
  dtrumat* A=(dtrumat*) AA;
  if (A && A->owndata){DSDPFREE(&A->val,&info);DSDPCHKERR(info);}
  if (A && A->sscale) {DSDPFREE(&A->sscale,&info);DSDPCHKERR(info);}
  if (A && A->workn) {DSDPFREE(&A->workn,&info);DSDPCHKERR(info);}
  if (A) {DSDPFREE(&A,&info);DSDPCHKERR(info);}
  return 0;
}

static int DTRUMatView(void* AA){
  dtrumat* M=(dtrumat*) AA;
  int i,j;
  double *val=M->val;
  ffinteger LDA=M->LDA;
  for (i=0; i<M->n; i++){
    for (j=0; j<=i; j++){
      printf(" %9.2e",val[i*LDA+j]);
    }
    for (j=i+1; j<M->LDA; j++){
      printf(" %9.1e",val[i*LDA+j]);
    }
    printf("\n");
  }
  return 0;
}

static int DTRUMatView2(void* AA){
  dtrumat* M=(dtrumat*) AA;
  int i,j;
  double *val=M->v2;
  ffinteger LDA=M->LDA;
  for (i=0; i<M->n; i++){
    for (j=0; j<=i; j++){
      printf(" %9.2e",val[i*LDA+j]);
    }
    for (j=i+1; j<M->LDA; j++){
      printf(" %9.2e",val[i*LDA+j]);
    }
    printf("\n");
  }
  return 0;
}


#include "dsdpschurmat_impl.h"
#include "dsdpdualmat_impl.h"
#include "dsdpdatamat_impl.h"
#include "dsdpxmat_impl.h"
#include "dsdpdsmat_impl.h"
#include "dsdpsys.h"

#undef __FUNCT__
#define __FUNCT__ "Tassemble"
static int DTRUMatAssemble(void*M){
  int info;
  double shift=1.0e-15;
  DSDPFunctionBegin;
  info= DTRUMatShiftDiagonal(M, shift); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int DTRUMatRowNonzeros(void*M, int row, double cols[], int *ncols,int nrows){
  int i;
  DSDPFunctionBegin;
  *ncols = row+1;
  for (i=0;i<=row;i++){
    cols[i]=1.0;
  }
  memset((void*)(cols+row+1),0,(nrows-row-1)*sizeof(int));
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TAddDiag"
static int DTRUMatAddDiag(void*M, int row, double dd){
  int ii;
  dtrumat*  ABA=(dtrumat*)M;
  ffinteger LDA=ABA->LDA;
  DSDPFunctionBegin;
  ii=row*LDA+row;
  ABA->val[ii] +=dd;
  DSDPFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "TAddDiag2"
static int DTRUMatAddDiag2(void*M, double diag[], int m){
  int row,ii;
  dtrumat*  ABA=(dtrumat*)M;
  ffinteger LDA=ABA->LDA;
  DSDPFunctionBegin;
  for (row=0;row<m;row++){
    ii=row*LDA+row;
    ABA->val[ii] +=diag[row];
  }
  DSDPFunctionReturn(0);
}
static struct  DSDPSchurMat_Ops dsdpmmatops;
static const char* lapackname="DENSE,SYMMETRIC U STORAGE";

static int DSDPInitSchurOps(struct  DSDPSchurMat_Ops* mops){ 
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatOpsInitialize(mops);DSDPCHKERR(info);
  mops->matrownonzeros=DTRUMatRowNonzeros;
  mops->matscaledmultiply=DTRUMatMult;
  mops->matmultr=DTRUMatMultR;
  mops->mataddrow=DTRUMatAddRow;
  mops->mataddelement=DTRUMatAddDiag;
  mops->matadddiagonal=DTRUMatAddDiag2;
  mops->matshiftdiagonal=DTRUMatShiftDiagonal;
  mops->matassemble=DTRUMatAssemble;
  mops->matfactor=DTRUMatCholeskyFactor;
  mops->matsolve=DTRUMatSolve;
  mops->matdestroy=DTRUMatDestroy;
  mops->matzero=DTRUMatZero;
  mops->matview=DTRUMatView;
  mops->id=1;
  mops->matname=lapackname;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPGetLAPACKSUSchurOps"
int DSDPGetLAPACKSUSchurOps(int n,struct DSDPSchurMat_Ops** sops, void** mdata){
  int info,nn,LDA;
  double *vv;
  dtrumat *AA;
  DSDPFunctionBegin;

  LDA=SUMatGetLDA(n);
  nn=n*LDA;
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DTRUMatCreateWData(n,LDA,vv,nn,&AA); DSDPCHKERR(info);
  AA->owndata=1;
  info=DSDPInitSchurOps(&dsdpmmatops); DSDPCHKERR(info);
  *sops=&dsdpmmatops;
  *mdata=(void*)AA;
  DSDPFunctionReturn(0);
}


static int DTRUMatCholeskyBackward(void* AA, double b[], double x[], int n){
  dtrumat* M=(dtrumat*) AA;
  ffinteger N=M->n,INCX=1,LDA=M->LDA;
  double *AP=M->val,*ss=M->sscale;
  char UPLO=M->UPLO,TRANS='N',DIAG='N';
  memcpy(x,b,N*sizeof(double));
  dtrsv(&UPLO,&TRANS,&DIAG, &N, AP, &LDA, x, &INCX );
  dtruscalevec(1.0,ss,x,x,n);
  return 0;
}


static int DTRUMatCholeskyForward(void* AA, double b[], double x[], int n){
  dtrumat* M=(dtrumat*) AA;
  ffinteger N=M->n,INCX=1,LDA=M->LDA;
  double *AP=M->val,*ss=M->sscale;
  char UPLO=M->UPLO,TRANS='T',DIAG='N';
  dtruscalevec(1.0,ss,b,x,n);
  dtrsv(&UPLO,&TRANS,&DIAG, &N, AP, &LDA, x, &INCX );
  return 0;
}

static int DTRUMatLogDet(void* AA, double *dd){
  dtrumat* A=(dtrumat*) AA;
  int i,n=A->n,LDA=A->LDA;
  double d=0,*v=A->val,*ss=A->sscale;
  for (i=0; i<n; i++){
    if (*v<=0) return 1;
    d+=2*log(*v/ss[i]);
    v+=LDA+1;
  }
  *dd=d;
  return 0;
}


static int DTRUMatCholeskyForwardMultiply(void* AA, double x[], double y[], int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger i,j,N=A->n,LDA=A->LDA;
  double *AP=A->val,*ss=A->sscale;
  /*  char UPLO=A->UPLO,TRANS='N',DIAG='N'; */
  
  if (x==0 && N>0) return 3;
  /*
  memcpy(y,x,N*sizeof(double));
  dtrmv(&UPLO,&TRANS,&DIAG,&N,AP,&LDA,Y,&INCY);
  */
  for (i=0;i<N;i++)y[i]=0;
  for (i=0;i<N;i++){
    for (j=0;j<=i;j++){
      y[i]+=AP[j]*x[j];
    }
    AP=AP+LDA;
  }

  for (i=0;i<N;i++){ y[i]=y[i]/ss[i];}
  return 0;
}

static int DTRUMatCholeskyBackwardMultiply(void* AA, double x[], double y[], int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger i,j,N=A->n,LDA=A->LDA;
  double *AP=A->val,*ss=A->sscale;
  /*  char UPLO=A->UPLO,TRANS='N',DIAG='N'; */
  
  if (x==0 && N>0) return 3;
  /*
  memcpy(y,x,N*sizeof(double));
  dtrmv(&UPLO,&TRANS,&DIAG,&N,AP,&LDA,Y,&INCY);
  */
  for (i=0;i<N;i++)y[i]=0;
  for (i=0;i<N;i++){
    for (j=0;j<=i;j++){
      y[j]+=AP[j]*x[i]/ss[i];
    }
    AP=AP+LDA;
  }

  for (i=0;i<-N;i++){ y[i]=y[i]/ss[i];}
  return 0;
}

static int DTRUMatInvert(void* AA){
  dtrumat* A=(dtrumat*) AA;
  ffinteger INFO,LDA=A->LDA,N=A->n,nn=LDA*N;
  double *v=A->val,*AP=A->v2,*ss=A->sscale;
  char UPLO=A->UPLO;
  memcpy((void*)AP,(void*)v,nn*sizeof(double));
  dpotri(&UPLO, &N, AP, &LDA, &INFO );
  if (INFO){
    INFO=DTRUMatShiftDiagonal(AA,1e-11); INFO=0;
    memcpy((void*)AP,(void*)v,nn*sizeof(double));
    dpotri(&UPLO, &N, AP, &LDA, &INFO );
  }
  if (A->scaleit){
    dtruscalemat(AP,ss,N,LDA);
  }
  A->status=Inverted;
  return INFO;

}

static void DTRUSymmetrize(dtrumat* A){
  int     i,j,n=A->n,row,LDA=A->LDA;
  double *v2=A->v2,*r1=A->v2,*r2=A->v2+LDA,*c1;
  for (i=0;i<n/2;i++){
    row=2*i;
    r1=v2+LDA*(row);
    r2=v2+LDA*(row+1);
    c1=v2+LDA*(row+1);
    r1[row+1]=c1[row];
    c1+=LDA;
    for (j=row+2;j<n;j++){
      r1[j]=c1[row];
      r2[j]=c1[row+1];
      c1+=LDA;
    }
    r1+=LDA*2;
    r2+=LDA*2;
  }
 
  for (row=2*n/2;row<n;row++){
    r1=v2+LDA*(row);
    c1=v2+LDA*(row+1);
    for (j=row+1;j<n;j++){
      r1[j]=c1[row];
      c1+=LDA;
    }
    r1+=LDA;
  }
  A->status=ISymmetric;
  return;
}

static int DTRUMatInverseAdd(void* AA, double alpha, double y[], int nn, int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger i,LDA=A->LDA,N,ione=1;
  double *v2=A->v2;
  for (i=0;i<n;i++){
    N=i+1;
    daxpy(&N,&alpha,v2,&ione,y,&ione);
    v2+=LDA; y+=n;
  }
  return 0;
}

static int DTRUMatInverseAddP(void* AA, double alpha, double y[], int nn, int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger N,LDA=A->LDA,i,ione=1;
  double *v2=A->v2;
  for (i=0;i<n;i++){
    N=i+1;
    daxpy(&N,&alpha,v2,&ione,y,&ione);
    v2+=LDA; y+=i+1;
  }
  return 0;
}

static void daddrow(double *v, double alpha, int i, double row[], int n){
  double *s1;
  ffinteger j,nn=n,ione=1;
  if (alpha==0){return;}
  nn=i+1; s1=v+i*n;
  daxpy(&nn,&alpha,s1,&ione,row,&ione);
  s1=v+i*n+n;
  for (j=i+1;j<n;j++){ row[j]+=alpha*s1[i]; s1+=n; }
  return;
}
/*
static void printrow(double r[], int n){int i;
    for (i=0;i<n;i++){printf(" %4.2e",r[i]);} printf("\n"); }
*/
static int DTRUMatInverseMultiply(void* AA, int indx[], int nind, double x[], double y[],int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger nn=n,LDA=A->LDA,N=A->n,INCX=1,INCY=1;
  double *AP=A->v2,*s1=A->v2,*s2,*X=x,*Y=y,ALPHA=1.0,BETA=0.0;
  char UPLO=A->UPLO,TRANS='N';
  int i,ii,usefull=1;
  
  if (usefull){
    if (A->status==Inverted){
      DTRUSymmetrize(A);
    }
    if (nind < n/4){
      memset((void*)y,0,n*sizeof(double));    
      for (ii=0;ii<nind;ii++){
	i=indx[ii]; nn=n; ALPHA=x[i];s2=s1+i*LDA;
	daxpy(&nn,&ALPHA,s2,&INCY,y,&INCX);
      }
    } else{
      ALPHA=1.0;
      dgemv(&TRANS,&N,&N,&ALPHA,AP,&LDA,X,&INCX,&BETA,Y,&INCY);
    }
    
  } else {
    if (nind<n/4 ){
      memset((void*)y,0,n*sizeof(double));    
      for (ii=0;ii<nind;ii++){
	i=indx[ii];  ALPHA=x[i];
	daddrow(s1,ALPHA,i,y,n);
      } 
    } else {
      ALPHA=1.0;
      dsymv(&UPLO,&N,&ALPHA,AP,&LDA,X,&INCX,&BETA,Y,&INCY);
    }
  }
  return 0;
}

static int DTRUMatSetXMat(void*A, double v[], int nn, int n){
  dtrumat*  ABA=(dtrumat*)A;
  int i,LDA=ABA->LDA;
  double *vv=ABA->val;
  if (vv!=v){
    for (i=0;i<n;i++){
      memcpy((void*)vv,(void*)v,(i+1)*sizeof(double));  
      vv+=LDA; v+=n;
    }
  }
  ABA->status=Assemble;
  return 0;
}
static int DTRUMatSetXMatP(void*A, double v[], int nn, int n){
  dtrumat*  ABA=(dtrumat*)A;
  int i,LDA=ABA->LDA;
  double *vv=ABA->val;
  if (vv!=v){
    for (i=0;i<n;i++){
      memcpy((void*)vv,(void*)v,(i+1)*sizeof(double));  
      v+=(i+1); vv+=LDA;
    }
  }
  ABA->status=Assemble;
  return 0;
}

static int DTRUMatFull(void*A, int*full){
  *full=1;
  return 0;
}

static int DTRUMatGetArray(void*A,double **v,int *n){
  dtrumat*  ABA=(dtrumat*)A;
  *n=ABA->n*ABA->LDA;
  *v=ABA->val;
  return 0;
}

static struct  DSDPDualMat_Ops sdmatops;
static int SDualOpsInitialize(struct  DSDPDualMat_Ops* sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDualMatOpsInitialize(sops);DSDPCHKERR(info);
  sops->matseturmat=DTRUMatSetXMat;
  sops->matgetarray=DTRUMatGetArray;
  sops->matcholesky=DTRUMatCholeskyFactor;
  sops->matsolveforward=DTRUMatCholeskyForward;
  sops->matsolvebackward=DTRUMatCholeskyBackward;
  sops->matinvert=DTRUMatInvert;
  sops->matinverseadd=DTRUMatInverseAdd;
  sops->matinversemultiply=DTRUMatInverseMultiply;
  sops->matforwardmultiply=DTRUMatCholeskyForwardMultiply;
  sops->matbackwardmultiply=DTRUMatCholeskyBackwardMultiply;
  sops->matfull=DTRUMatFull;
  sops->matdestroy=DTRUMatDestroy;
  sops->matgetsize=DTRUMatGetSize;
  sops->matview=DTRUMatView;
  sops->matlogdet=DTRUMatLogDet;
  sops->matname=lapackname;
  sops->id=1;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACKSUDualMatCreate"
static int DSDPLAPACKSUDualMatCreate(int n,struct  DSDPDualMat_Ops **sops, void**smat){
  dtrumat *AA;
  int info,nn,LDA=n;
  double *vv;
  DSDPFunctionBegin;
  LDA=SUMatGetLDA(n);
  nn=n*LDA;
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DTRUMatCreateWData(n,LDA,vv,nn,&AA); DSDPCHKERR(info);
  AA->owndata=1;
  info=SDualOpsInitialize(&sdmatops);DSDPCHKERR(info);
  *sops=&sdmatops;
  *smat=(void*)AA;
  DSDPFunctionReturn(0);
}


static int switchptr(void *SD,void *SP){
  dtrumat *s1,*s2;
  s1=(dtrumat*)(SD);
  s2=(dtrumat*)(SP);
  s1->v2=s2->val;
  s2->v2=s1->val;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACKSUDualMatCreate2"
int DSDPLAPACKSUDualMatCreate2(int n, 
			       struct  DSDPDualMat_Ops **sops1, void**smat1,
			       struct  DSDPDualMat_Ops **sops2, void**smat2){
  int info;
  DSDPFunctionBegin;
  info=DSDPLAPACKSUDualMatCreate(n,sops1,smat1);DSDPCHKERR(info);
  info=DSDPLAPACKSUDualMatCreate(n,sops2,smat2);DSDPCHKERR(info);
  info=switchptr(*smat1,*smat2);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static struct  DSDPDualMat_Ops sdmatopsp;
static int SDualOpsInitializeP(struct  DSDPDualMat_Ops* sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDualMatOpsInitialize(sops);DSDPCHKERR(info);
  sops->matseturmat=DTRUMatSetXMatP;
  sops->matgetarray=DTRUMatGetArray;
  sops->matcholesky=DTRUMatCholeskyFactor;
  sops->matsolveforward=DTRUMatCholeskyForward;
  sops->matsolvebackward=DTRUMatCholeskyBackward;
  sops->matinvert=DTRUMatInvert;
  sops->matinverseadd=DTRUMatInverseAddP;
  sops->matinversemultiply=DTRUMatInverseMultiply;
  sops->matforwardmultiply=DTRUMatCholeskyForwardMultiply;
  sops->matbackwardmultiply=DTRUMatCholeskyBackwardMultiply;
  sops->matfull=DTRUMatFull;
  sops->matdestroy=DTRUMatDestroy;
  sops->matgetsize=DTRUMatGetSize;
  sops->matview=DTRUMatView;
  sops->matlogdet=DTRUMatLogDet;
  sops->matname=lapackname;
  sops->id=1;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACKSUDualMatCreate"
static int DSDPLAPACKSUDualMatCreateP(int n,  struct  DSDPDualMat_Ops **sops, void**smat){
  dtrumat *AA;
  int info,nn,LDA;
  double *vv;
  DSDPFunctionBegin;
  LDA=SUMatGetLDA(n);
  nn=LDA*n;
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DTRUMatCreateWData(n,LDA,vv,nn,&AA); DSDPCHKERR(info);
  AA->owndata=1;
  info=SDualOpsInitializeP(&sdmatopsp);DSDPCHKERR(info);
  *sops=&sdmatopsp;
  *smat=(void*)AA;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPLAPACKSUDualMatCreate2P"
int DSDPLAPACKSUDualMatCreate2P(int n, 
			       struct  DSDPDualMat_Ops* *sops1, void**smat1,
			       struct  DSDPDualMat_Ops* *sops2, void**smat2){
  int info;
  DSDPFunctionBegin;
  info=DSDPLAPACKSUDualMatCreateP(n,sops1,smat1);
  info=DSDPLAPACKSUDualMatCreateP(n,sops2,smat2);
  info=switchptr(*smat1,*smat2);
  DSDPFunctionReturn(0);
}

static int DTRUMatScaleDiagonal(void* AA, double dd){
  dtrumat* A=(dtrumat*) AA;
  ffinteger LDA=A->LDA;
  int i,n=A->n;
  double *v=A->val;
  for (i=0; i<n; i++){
    *v*=dd;
    v+=LDA+1;    
  }
  return 0;
}

static int DTRUMatOuterProduct(void* AA, double alpha, double x[], int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger ione=1,N=n,LDA=A->LDA;
  double *v=A->val;
  char UPLO=A->UPLO;
  dsyr(&UPLO,&N,&alpha,x,&ione,v,&LDA);
  return 0;
}

static int DenseSymPSDNormF2(void* AA, int n, double *dddot){
  dtrumat* A=(dtrumat*) AA;
  ffinteger ione=1,nn=A->n*A->n;
  double dd,tt=sqrt(0.5),*val=A->val;
  int info;
  info=DTRUMatScaleDiagonal(AA,tt);
  dd=dnrm2(&nn,val,&ione);
  info=DTRUMatScaleDiagonal(AA,1.0/tt);
  *dddot=dd*dd*2;
  return 0;
}

static int DTRUMatGetDenseArray(void* A, double *v[], int*n){
  dtrumat*  ABA=(dtrumat*)A;
  *v=ABA->val;
  *n=ABA->n*ABA->LDA;
  return 0;
}

static int DTRUMatRestoreDenseArray(void* A, double *v[], int *n){
  *v=0;*n=0;
  return 0;
}

static int DDenseSetXMat(void*A, double v[], int nn, int n){
  dtrumat*  ABA=(dtrumat*)A;
  int i,LDA=ABA->LDA;
  double *vv=ABA->val;
  if (v!=vv){
    for (i=0;i<n;i++){
      memcpy((void*)vv,(void*)v,nn*sizeof(double));  
      v+=n;vv+=LDA;
    }
  }
  ABA->status=Assemble;
  return 0;
}

static int DDenseVecVec(void* A, double x[], int n, double *v){
  dtrumat*  ABA=(dtrumat*)A;
  int i,j,k=0,LDA=ABA->LDA;
  double dd=0,*val=ABA->val;
  *v=0.0;
  for (i=0; i<n; i++){
    for (j=0;j<i;j++){
      dd+=2*x[i]*x[j]*val[j];
      k++;
    }
    dd+=x[i]*x[i]*val[i];
    k+=LDA;
  }
  *v=dd;
  return 0;
}



static int DTRUMatEigs(void*AA,double W[],double IIWORK[], int nn1, double *mineig){
  dtrumat* AAA=(dtrumat*) AA;
  ffinteger info,INFO=0,M,N=AAA->n;
  ffinteger IL=1,IU=1,LDA=AAA->LDA,LDZ=LDA,LWORK,IFAIL;
  ffinteger *IWORK=(ffinteger*)IIWORK;
  double *AP=AAA->val;
  double Z=0,VL=-1e10,VU=1e10,*WORK,ABSTOL=1e-13;
  char UPLO=AAA->UPLO,JOBZ='N',RANGE='I';
  DSDPCALLOC2(&WORK,double,8*N,&info);
  DSDPCALLOC2(&IWORK,ffinteger,5*N,&info);
  LWORK=8*N;
  dsyevx(&JOBZ,&RANGE,&UPLO,&N,AP,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,&Z,&LDZ,WORK,&LWORK,IWORK,&IFAIL,&INFO);
  /*
  ffinteger LIWORK=nn1;
  dsyevd(&JOBZ,&UPLO,&N,AP,&LDA,W,WORK,&LWORK,IWORK,&LIWORK,&INFO);
  */
  DSDPFREE(&WORK,&info);
  DSDPFREE(&IWORK,&info);
  *mineig=W[0];
  return INFO;
}


static struct  DSDPVMat_Ops turdensematops;

static int DSDPDenseXInitializeOps(struct  DSDPVMat_Ops* densematops){
  int info;
  if (!densematops) return 0;
  info=DSDPVMatOpsInitialize(densematops); DSDPCHKERR(info);
  densematops->matview=DTRUMatView;
  densematops->matscalediagonal=DTRUMatScaleDiagonal;
  densematops->matshiftdiagonal=DTRUMatShiftDiagonal;
  densematops->mataddouterproduct=DTRUMatOuterProduct;
  densematops->matmult=DTRUMatMult;
  densematops->matdestroy=DTRUMatDestroy;
  densematops->matfnorm2=DenseSymPSDNormF2;
  densematops->matgetsize=DTRUMatGetSize;
  densematops->matzeroentries=DTRUMatZero;
  densematops->matgeturarray=DTRUMatGetDenseArray;
  densematops->matrestoreurarray=DTRUMatRestoreDenseArray;
  densematops->matmineig=DTRUMatEigs;
  densematops->id=1;
  densematops->matname=lapackname;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPXMatUCreateWithData"
int DSDPXMatUCreateWithData(int n,double nz[],int nnz,struct  DSDPVMat_Ops* *xops, void * *xmat){
  int i,info;
  double dtmp;
  dtrumat*AA;
  DSDPFunctionBegin;
  if (nnz<n*n){DSDPSETERR1(2,"Array must have length of : %d \n",n*n);}
  for (i=0;i<n*n;i++) dtmp=nz[i];
  info=DTRUMatCreateWData(n,n,nz,nnz,&AA); DSDPCHKERR(info);
  AA->owndata=0;
  info=DSDPDenseXInitializeOps(&turdensematops); DSDPCHKERR(info);
  *xops=&turdensematops;
  *xmat=(void*)AA;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPXMatUCreate"
int DSDPXMatUCreate(int n,struct  DSDPVMat_Ops* *xops, void * *xmat){
  int info,nn=n*n;
  double *vv;
  DSDPFunctionBegin;
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DSDPXMatUCreateWithData(n,vv,nn,xops,xmat);DSDPCHKERR(info);
  DTRUMatOwnData((dtrumat*)(*xmat),1);
  DSDPFunctionReturn(0);
}

static struct  DSDPDSMat_Ops tdsdensematops;
static int DSDPDSDenseInitializeOps(struct  DSDPDSMat_Ops* densematops){
  int info;
  if (!densematops) return 0;
  info=DSDPDSMatOpsInitialize(densematops); DSDPCHKERR(info);
  densematops->matseturmat=DDenseSetXMat;
  densematops->matview=DTRUMatView;
  densematops->matdestroy=DTRUMatDestroy;
  densematops->matgetsize=DTRUMatGetSize;
  densematops->matzeroentries=DTRUMatZero;
  densematops->matmult=DTRUMatMult;
  densematops->matvecvec=DDenseVecVec;
  densematops->id=1;
  densematops->matname=lapackname;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateDSMatWithArray2"
int DSDPCreateDSMatWithArray2(int n,double vv[],int nnz,struct  DSDPDSMat_Ops* *dsmatops, void**dsmat){
  int info;
  dtrumat*AA;
  DSDPFunctionBegin;
  info=DTRUMatCreateWData(n,n,vv,nnz,&AA); DSDPCHKERR(info);
  AA->owndata=0;
  info=DSDPDSDenseInitializeOps(&tdsdensematops); DSDPCHKERR(info);
  *dsmatops=&tdsdensematops;
  *dsmat=(void*)AA;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateXDSMat2"
int DSDPCreateXDSMat2(int n,struct  DSDPDSMat_Ops* *dsmatops, void**dsmat){
  int info,nn=n*n;
  double *vv;
  DSDPFunctionBegin;  
  DSDPCALLOC2(&vv,double,nn,&info);DSDPCHKERR(info);
  info=DSDPCreateDSMatWithArray2(n,vv,nn,dsmatops,dsmat);DSDPCHKERR(info);
  DTRUMatOwnData((dtrumat*)(*dsmat),1);
  DSDPFunctionReturn(0);
}


typedef struct {
  int    neigs;
  double *eigval;
  double *an;
} Eigen;

typedef struct {
  dtrumat* AA;
  Eigen   *Eig;
} dvecumat;

#undef __FUNCT__  
#define __FUNCT__ "CreateDvecumatWData"
static int CreateDvecumatWdata(int n, double vv[], dvecumat **A){
  int info,nnz=n*n;
  dvecumat* V;
  DSDPCALLOC1(&V,dvecumat,&info);DSDPCHKERR(info);
  info=DTRUMatCreateWData(n,n,vv,nnz,&V->AA); DSDPCHKERR(info);
  V->Eig=0;
  *A=V;
  return 0;
}


static int DvecumatGetRowNnz(void* AA, int trow, int nz[], int *nnzz,int n){
  int k;
  *nnzz=n;
  for (k=0;k<n;k++) nz[k]++;
  return 0;
}

static int DTRUMatGetRowAdd(void* AA, int nrow, double ytmp, double row[], int n){
  dtrumat* A=(dtrumat*) AA;
  ffinteger i,nnn=n;
  double *v=A->val;

  nnn=nrow*n;
  for (i=0;i<=nrow;i++){
    row[i]+=ytmp*v[nnn+i];
  }
  for (i=nrow+1;i<n;i++){
    nnn+=nrow;
    row[i]+=ytmp*v[nrow];
  }
  return 0;
}

static int DvecumatGetRowAdd(void* AA, int trow, double scl, double r[], int m){
  int info;
  dvecumat* A=(dvecumat*)AA;
  info=DTRUMatGetRowAdd((void*)A->AA ,trow,scl,r,m);
  return 0;
}

static int DvecumatAddMultiple(void* AA, double alpha, double r[], int nnn, int n){
  dvecumat* A=(dvecumat*)AA;
  ffinteger nn=nnn, ione=1;
  double *val=A->AA->val;
  daxpy(&nn,&alpha,val,&ione,r,&ione);
  return 0;
}


static int DvecuEigVecVec(void*, double[], int, double*);
static int DvecumatVecVec(void* AA, double x[], int n, double *v){
  dvecumat* A=(dvecumat*)AA;
  int i,j,k=0,LDA=A->AA->LDA;
  double dd=0,*val=A->AA->val;
  *v=0.0;
  if (A->Eig && A->Eig->neigs<n/5){
    i=DvecuEigVecVec(AA,x,n,v);
  } else {
    for (i=0; i<n; i++){
      for (j=0;j<i;j++){
	dd+=2*x[i]*x[j]*val[j];
      }
      dd+=x[i]*x[i]*val[i];
      k+=LDA;
    }
    *v=dd;
  }
  return 0;
}


static int DvecumatFNorm2(void* AA, int n, double *v){
  dvecumat* A=(dvecumat*)AA;
  long int i,j,k=0,LDA=A->AA->LDA;
  double dd=0,*x=A->AA->val;
  for (i=0; i<n; i++){
    for (j=0;j<i;j++){
      dd+=2*x[j]*x[j];
    }
    dd+=x[i]*x[i];
    k+=LDA;
  }
  *v=dd;
  return 0;
}


static int DvecumatCountNonzeros(void* AA, int *nnz, int n){
  *nnz=n*(n+1)/2;
  return 0;
}


static int DvecumatDot(void* AA, double x[], int nn, int n, double *v){
  dvecumat* A=(dvecumat*)AA;
  double d1,dd=0,*v1=x,*v2=A->AA->val;
  ffinteger i,n2,ione=1,LDA=A->AA->LDA;

  for (i=0;i<n;i++){
    n2=i+1; 
    d1=ddot(&n2,v1,&ione,v2,&ione);
    v1+=n; v2+=LDA;
    dd+=d1;
  }
  *v=2*dd;
  return 0;
}

/*
static int DvecumatNormF2(void* AA, int n, double *v){
  dvecumat* A=(dvecumat*)AA;
  return(DTRUMatNormF2((void*)(A->AA), n,v));
}
*/
#undef __FUNCT__  
#define __FUNCT__ "DvecumatDestroy"
static int DvecumatDestroy(void* AA){
  dvecumat* A=(dvecumat*)AA;
  int info;
  info=DTRUMatDestroy((void*)(A->AA));
  if (A->Eig){
    DSDPFREE(&A->Eig->an,&info);DSDPCHKERR(info);
    DSDPFREE(&A->Eig->eigval,&info);DSDPCHKERR(info);
  }
  DSDPFREE(&A->Eig,&info);DSDPCHKERR(info);
  DSDPFREE(&A,&info);DSDPCHKERR(info);
  return 0;
}


static int DvecumatView(void* AA){
  dvecumat* A=(dvecumat*)AA;
  dtrumat* M=A->AA;
  int i,j,LDA=M->LDA;
  double *val=M->val;
  for (i=0; i<M->n; i++){
    for (j=0; j<M->n; j++){
      printf(" %4.2e",val[j]);
    }
    val+=LDA;
  }
  return 0;
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPCreateDvecumatEigs"
static int CreateEigenLocker(Eigen **EE,int neigs, int n){
  int info;
  Eigen *E;

  DSDPCALLOC1(&E,Eigen,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&E->eigval,double,neigs,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&E->an,double,n*neigs,&info);DSDPCHKERR(info);
  E->neigs=neigs;
  *EE=E;
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


static int DvecumatComputeEigs(dvecumat*,double[],int,double[],int,double[],int,int[],int);

static int DvecumatFactor(void*AA, double dmatp[], int nn0, double dwork[], int n, double ddwork[], int n1, int iptr[], int n2){
  
  int info;
  dvecumat*  A=(dvecumat*)AA;
  if (A->Eig) return 0;
  info=DvecumatComputeEigs(A,dmatp,nn0,dwork,n,ddwork,n1,iptr,n2);DSDPCHKERR(info);
  return 0;
}

static int DvecumatGetRank(void *AA,int *rank, int n){
  dvecumat*  A=(dvecumat*)AA;
  if (A->Eig){
    *rank=A->Eig->neigs;
  } else {
    DSDPSETERR(1,"Vecu Matrix not factored yet\n");
  }
  return 0;
}

static int DvecumatGetEig(void* AA, int rank, double *eigenvalue, double vv[], int n, int indz[], int *nind){
  dvecumat*  A=(dvecumat*)AA;
  int i,info;
  if (A->Eig){
    info=EigMatGetEig(A->Eig,rank,eigenvalue,vv,n);DSDPCHKERR(info);
    *nind=n;
    for (i=0;i<n;i++){ indz[i]=i;}
  } else {
    DSDPSETERR(1,"Vecu Matrix not factored yet\n");
  }
  return 0;  
}

static int DvecuEigVecVec(void* AA, double v[], int n, double *vv){
  dvecumat*  A=(dvecumat*)AA;
  int i,rank,neigs;
  double *an,dd,ddd=0,*eigval;
  if (A->Eig){
    an=A->Eig->an;
    neigs=A->Eig->neigs;
    eigval=A->Eig->eigval;
    for (rank=0;rank<neigs;rank++){
      for (dd=0,i=0;i<n;i++){
	dd+=v[i]*an[i];
      }
      an+=n;
      ddd+=dd*dd*eigval[rank];
    }
    *vv=ddd;
  } else {
    DSDPSETERR(1,"Vecu Matrix not factored yet\n");
  }
  return 0;  
}


static struct  DSDPDataMat_Ops dvecumatops;
static const char *datamatname="STANDARD VECU MATRIX";

static int DvecumatOpsInitialize(struct  DSDPDataMat_Ops *sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(sops); DSDPCHKERR(info);
  sops->matvecvec=DvecumatVecVec;
  sops->matdot=DvecumatDot;
  sops->mataddrowmultiple=DvecumatGetRowAdd;
  sops->mataddallmultiple=DvecumatAddMultiple;
  sops->matview=DvecumatView;
  sops->matdestroy=DvecumatDestroy;
  sops->matfactor2=DvecumatFactor;
  sops->matgetrank=DvecumatGetRank;
  sops->matgeteig=DvecumatGetEig;
  sops->matrownz=DvecumatGetRowNnz;
  sops->matfnorm2=DvecumatFNorm2;
  sops->matnnz=DvecumatCountNonzeros;
  sops->id=1;
  sops->matname=datamatname;
  return 0;
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDUmat"
int DSDPGetDUMat(int n,double *val, struct  DSDPDataMat_Ops**sops, void**smat){ 
  int info,k;
  double dtmp;
  dvecumat* A;
  DSDPFunctionBegin;

  for (k=0;k<n*n;++k) dtmp=val[k];
  info=CreateDvecumatWdata(n,val,&A); DSDPCHKERR(info);
  A->Eig=0;
  info=DvecumatOpsInitialize(&dvecumatops); DSDPCHKERR(info);
  if (sops){*sops=&dvecumatops;}
  if (smat){*smat=(void*)A;}
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DvecumatComputeEigs"
static int DvecumatComputeEigs(dvecumat* AA,double DD[], int nn0, double W[], int n, double WORK[], int n1, int iiptr[], int n2){

  int i,neigs,info;
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
  memcpy((void*)dmatarray,(void*)val,n*n*sizeof(double));

  if (n*n>nn2){
    DSDPCALLOC2(&dworkarray,double,(n*n),&info); DSDPCHKERR(info);
    ownarray2=1;
  }

  if (n*n*sizeof(long int)>nn0*sizeof(double)){
    DSDPCALLOC2(&i2darray,long int,(n*n),&info); DSDPCHKERR(info);
    ownarray3=1;
  }


  /* Call LAPACK to compute the eigenvalues */
  info=DSDPGetEigs(dmatarray,n,dworkarray,n*n,i2darray,n*n,
		   W,n,WORK,n1,iiptr,n2);
  if (info){
    memcpy((void*)dmatarray,(void*)val,n*n*sizeof(double));
    info=DSDPGetEigs2(dmatarray,n,dworkarray,n*n,i2darray,n*n,
		      W,n,WORK,n1,iiptr+3*n,n2-3*n); DSDPCHKERR(info);
  } 

  /* Count the nonzero eigenvalues */
  for (neigs=0,i=0;i<n;i++){
    if (fabs(W[i])> eps ){ neigs++;}
  }

  info=CreateEigenLocker(&AA->Eig,neigs,n);DSDPCHKERR(info);
  
  /* Copy into structure */
  for (neigs=0,i=0; i<n; i++){
    if (fabs(W[i]) > eps){
      info=EigMatSetEig(AA->Eig,neigs,W[i],dmatarray+n*i,n);DSDPCHKERR(info);
      neigs++;
    }
  }
  
  if (ownarray1){ DSDPFREE(&dmatarray,&info);DSDPCHKERR(info);}
  if (ownarray2){ DSDPFREE(&dworkarray,&info);DSDPCHKERR(info);}
  if (ownarray3){ DSDPFREE(&i2darray,&info);DSDPCHKERR(info);}
  return 0;
}

