#include "dsdpdatamat_impl.h"
#include "dsdpsys.h"
/*!
\file rmmat.c
\brief DSDPDataMat object of rank one outer product.
*/

typedef struct {
  double ev;
  const double *spval;
  const int    *spai;
  int    nnz;
  int    n;
  int    ishift;
  char   UPLQ;
} r1mat;

static int R1MatDestroy(void*);
static int R1MatView(void*);
static int R1MatVecVec(void*, double[], int, double *);
static int R1MatDotP(void*, double[],int,int,double *);
static int R1MatDotU(void*, double[],int,int,double *);
static int R1MatGetRank(void*, int*, int);
static int R1MatFactor(void*);
static int R1MatGetEig(void*, int, double*, double[], int,int[],int*);
static int R1MatRowNnz(void*, int, int[], int*, int);
static int R1MatAddRowMultiple(void*, int, double, double[], int);
static int R1MatAddMultipleP(void*, double, double[], int,int);
static int R1MatAddMultipleU(void*, double, double[], int,int);

static struct  DSDPDataMat_Ops r1matopsP;
static struct  DSDPDataMat_Ops r1matopsU;
static int R1MatOpsInitializeP(struct  DSDPDataMat_Ops*);
static int R1MatOpsInitializeU(struct  DSDPDataMat_Ops*);


#undef __FUNCT__  
#define __FUNCT__ "DSDPGetR1Mat"
int DSDPGetR1Mat(int n, double ev, int ishift, const int spai[], const double spval[], int nnz, char UPLQ, void**mmat){
  int i;
  r1mat*AA;
  DSDPFunctionBegin;
  for (i=0;i<nnz;i++){
    if (spai[i]-ishift<0 || spai[i]-ishift >=n){
      printf("Invalid entry: Entry %d . Is %d <= %d < %d?\n",i,ishift,spai[i],n+ishift);
      return 1;
    }
  }
  AA=(r1mat*) malloc(1*sizeof(r1mat));
  if (AA==NULL) return 1;
  AA->n=n;
  AA->UPLQ=UPLQ;
  AA->spval=spval;
  AA->spai=spai;
  AA->nnz=nnz;
  AA->ev=ev;
  AA->ishift=ishift;
  if (mmat){*mmat=(void*)AA;}
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPGetR1PMat"
/*!
int DSDPGetR1PMat(int n, double ev, int ishift, const int spai[], const double spval[], int nnz, struct  DSDPDataMat_Ops**mops, void**mmat)

\brief Create a rank one matrix usuable by DSDP in packed symmetric format.
\param n number of rows and columns of the matrix
\param ev multiple of the outer product.
\param ishift index of first element in vector.
\param spai array of indices for vector.
\param spval array of vector values.
\param nnz size of arrays.
\param mops address of a pointer to a table of function pointers
\param mmat address of a pointer to an opaque data type.
*/
int DSDPGetR1PMat(int n, double ev, int ishift, const int spai[], const double spval[], int nnz, struct  DSDPDataMat_Ops**mops, void**mmat){
  int info;
  DSDPFunctionBegin;
  info=DSDPGetR1Mat(n,ev,ishift,spai,spval,nnz,'P',mmat);
  info=R1MatOpsInitializeP(&r1matopsP); if(info){return 1;}
  if (mops){*mops=&r1matopsP;}
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetR1UMat"
/*!
int DSDPGetR1UMat(int n, double ev, int ishift, const int spai[], const double spval[], int nnz, struct  DSDPDataMat_Ops**mops, void**mmat)

\brief Create a rank one matrix usuable by DSDP in full symmetric format.
\param n number of rows and columns of the matrix
\param ev multiple of the outer product.
\param ishift index of first element in vector.
\param spai array of indices for vector.
\param spval array of vector values.
\param nnz size of arrays.
\param mops address of a pointer to a table of function pointers
\param mmat address of a pointer to an opaque data type.
*/
int DSDPGetR1UMat(int n, double ev, int ishift, const int spai[], const double spval[], int nnz, struct  DSDPDataMat_Ops**mops, void**mmat){
  int info;
  DSDPFunctionBegin;
  info=DSDPGetR1Mat(n,ev,ishift,spai,spval,nnz,'U',mmat);
  info=R1MatOpsInitializeU(&r1matopsU); if(info){return 1;}
  if (mops){*mops=&r1matopsU;}
  DSDPFunctionReturn(0);
}

static int R1MatDotP(void* A, double x[], int nn, int n, double *v){
  r1mat* AA = (r1mat*)A;
  int i,i2,i3,j,j2;
  int nnz=AA->nnz,ishift=AA->ishift;
  const int *ai=AA->spai;
  double dtmp=0.0,d3;
  const double *val=AA->spval;
  for (i=0;i<nnz;i++){
    d3=val[i];
    i2=ai[i]-ishift;
    i3=(i2+1)*i2/2;
    for (j=0;j<nnz;j++){
      j2=ai[j]-ishift;
      if (j2<=i2){
	dtmp+=2*x[i3+j2]*d3*val[j];
      } 
    }
  }
  *v=dtmp*AA->ev;
  return 0;
}

static int R1MatDotU(void* A, double x[], int nn, int n, double *v){
  r1mat* AA = (r1mat*)A;
  int i,i2,i3,j,j2;
  int nnz=AA->nnz,ishift=AA->ishift;
  const int *ai=AA->spai;
  const double *val=AA->spval;
  double dtmp=0.0,d3;

  for (i=0;i<nnz;i++){
    d3=val[i];
    i2=ai[i]-ishift;
    i3=i2*n;
    for (j=0;j<nnz;j++){
      j2=ai[j]-ishift;
      if (j2<=i2){
	dtmp+=2*x[i3+j2]*d3*val[j];
      } 
    }
  }
  *v=dtmp*AA->ev;
  return 0;
}

static int R1MatVecVec(void* A, double x[], int n, double *v){

  r1mat* AA = (r1mat*)A;
  double dtmp=0.0;
  const double *val=AA->spval;
  int i,ishift=AA->ishift,nnz=AA->nnz;
  const int *ai=AA->spai;
  for (i=0; i<nnz; i++){
    dtmp+=val[i] * x[ai[i]-ishift];
  }
  *v=dtmp*dtmp*AA->ev;
  return 0;
}

static int R1MatAddMultipleP(void*A, double dd, double vv[], int nn, int n){
  r1mat* AA = (r1mat*)A;
  int i,i2,i3,j,j2;
  int nnz=AA->nnz,ishift=AA->ishift;
  const int *ai=AA->spai;
  const double *val=AA->spval;
  double d3,ddd=dd*AA->ev;
  for (i=0;i<nnz;i++){
    d3=ddd*val[i];
    i2=ai[i]-ishift;
    i3=(i2+1)*i2/2;
    for (j=0;j<nnz;j++){
      j2=ai[j]-ishift;
      if (j2<=i2){
	vv[i3+j2]+=d3*val[j];
      } 
    }
  }
  return 0;
}
static int R1MatAddMultipleU(void*A, double dd, double vv[], int nn, int n){
  r1mat* AA = (r1mat*)A;
  int i,i2,i3,j,j2;
  int nnz=AA->nnz,ishift=AA->ishift;
  const int *ai=AA->spai;
  const double *val=AA->spval;
  double d3,ddd=dd*AA->ev;
  for (i=0;i<nnz;i++){
    d3=ddd*val[i];
    i2=ai[i]-ishift;
    i3=i2*n;      
    for (j=0;j<nnz;j++){
      j2=ai[j]-ishift;
      if (j2<=i2){
	vv[i3+j2]+=d3*val[j];
      } 
    }
  }
  return 0;
}

static int R1MatAddRowMultiple(void*A, int nrow, double dd, double row[], int n){
  r1mat* AA = (r1mat*)A;
  int nnz=AA->nnz,ishift=AA->ishift;
  const int *ai=AA->spai;
  const double *val=AA->spval;
  double ddd=dd*AA->ev;
  int i,j;
  for (i=0;i<nnz;i++){ 
    if (ai[i]-ishift==nrow){
      for (j=0;j<nnz;j++){
	row[ai[j]-ishift]+= ddd*val[i]*val[j];
      }
    }
  }
  return 0;
}


static int R1MatFactor(void*A){
  return 0;
}


static int R1MatGetRank(void *A, int*rank, int n){
  *rank=1;
  return 0;
}

static int R1MatGetEig(void*A, int neig, double *eig, double v[], int n, int  indx[], int*nind){
  r1mat* AA = (r1mat*)A;
  int i,aii,ishift=AA->ishift,nnz=AA->nnz;
  const int *ai=AA->spai;
  const double *val=AA->spval;
  for (i=0;i<n;i++){ v[i]=0.0; }
  *eig=0; *nind=0;
  if (neig==0){
    for (i=0;i<nnz;i++){ 
      aii=ai[i]-ishift;
      v[aii]=val[i]; 
      indx[i]=aii;
    }
    *eig=AA->ev; *nind=AA->nnz;
  }
  return 0;
}


static int R1MatRowNnz(void*A, int row, int nz[], int *rnnz, int n){
  r1mat* AA = (r1mat*)A;
  int i,j;
  int nnz=AA->nnz,ishift=AA->ishift;
  const int *ai=AA->spai;
  *rnnz=0;
  for (i=0;i<nnz;i++){ 
    if (ai[i]-ishift==row){
      for (j=0;j<nnz;j++){
	nz[ai[j]-ishift]++;
      }
    }
    *rnnz=nnz;
  }
  return 0;
}

static int R1MatFNorm2(void*A, int n, double *fnorm2){
  r1mat* AA = (r1mat*)A;
  double dd=0;
  const double *val=AA->spval;
  int i,nnz=AA->nnz;
  for (i=0;i<nnz;i++){
    dd+=val[i]*val[i];
  }
  *fnorm2=dd*dd*AA->ev*AA->ev;
  return 0;
}

static int R1MatCountNonzeros(void*A, int *nnz, int n){
  r1mat* AA = (r1mat*)A;
  *nnz=AA->nnz*AA->nnz;
  return 0;
}


static int R1MatView(void* A){
  int i;
  r1mat* AA = (r1mat*)A;
  printf("This matrix is %4.8e times the outer product of \n",AA->ev);
  for (i=0;i<AA->nnz;i++){
    printf("%d  %4.8e \n",AA->spai[i],AA->spval[i]);
  }
  return 0;
}


static int R1MatDestroy(void* A){
  if (A) free(A);
  return 0;
}

static const char *datamatname="RANK 1 Outer Product";
static int R1MatOpsInitializeP(struct  DSDPDataMat_Ops* r1matops){
  int info;
  if (r1matops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(r1matops); DSDPCHKERR(info);
  r1matops->matfactor1=R1MatFactor;
  r1matops->matgetrank=R1MatGetRank;
  r1matops->matgeteig=R1MatGetEig;
  r1matops->matvecvec=R1MatVecVec;
  r1matops->matdot=R1MatDotP;
  r1matops->mataddrowmultiple=R1MatAddRowMultiple;
  r1matops->mataddallmultiple=R1MatAddMultipleP;
  r1matops->matdestroy=R1MatDestroy;
  r1matops->matview=R1MatView;
  r1matops->matrownz=R1MatRowNnz;
  r1matops->matfnorm2=R1MatFNorm2;
  r1matops->matnnz=R1MatCountNonzeros;
  r1matops->id=15;
  r1matops->matname=datamatname;
  return 0;
}
static int R1MatOpsInitializeU(struct  DSDPDataMat_Ops* r1matops){
  int info;
  if (r1matops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(r1matops); DSDPCHKERR(info);
  r1matops->matfactor1=R1MatFactor;
  r1matops->matgetrank=R1MatGetRank;
  r1matops->matgeteig=R1MatGetEig;
  r1matops->matvecvec=R1MatVecVec;
  r1matops->matdot=R1MatDotU;
  r1matops->mataddrowmultiple=R1MatAddRowMultiple;
  r1matops->mataddallmultiple=R1MatAddMultipleU;
  r1matops->matdestroy=R1MatDestroy;
  r1matops->matview=R1MatView;
  r1matops->matrownz=R1MatRowNnz;
  r1matops->matfnorm2=R1MatFNorm2;
  r1matops->matnnz=R1MatCountNonzeros;
  r1matops->id=15;
  r1matops->matname=datamatname;
  return 0;
}

