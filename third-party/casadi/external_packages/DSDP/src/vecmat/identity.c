#include "dsdpdatamat_impl.h"
#include "dsdpsys.h"
/*!
\file identity.c
\brief DSDPDataMat object representing a multiple of the identity matrix.
*/

typedef struct {
  int n;
  double dm;
} identitymat;


static int IdentityMatDestroy(void*);
static int IdentityMatView(void*);
static int IdentityMatVecVec(void*, double[], int, double *);
static int IdentityMatDotP(void*, double[], int, int, double*);
static int IdentityMatDotF(void*, double[], int, int, double*);
static int IdentityMatGetRank(void*, int*,int);
static int IdentityMatFactor(void*);
static int IdentityMatGetEig(void*, int, double*, double[], int,int[],int*);
static int IdentityMatAddRowMultiple(void*, int, double, double[], int);
static int IdentityMatAddMultipleP(void*, double, double[], int, int);
static int IdentityMatAddMultipleF(void*, double, double[], int, int);
static int IdentityMatGetRowNnz(void*, int, int[], int*, int);

static struct  DSDPDataMat_Ops identitymatopsp;
static struct  DSDPDataMat_Ops identitymatopsf;
static int IdentitymatOperationsInitializeP(struct  DSDPDataMat_Ops*);
static int IdentitymatOperationsInitializeF(struct  DSDPDataMat_Ops*);


/*!
\fn int DSDPGetIdentityDataMatP(int n, double val, struct  DSDPDataMat_Ops** dops, void** imat)

\brief Create a sparse matrix usuable by DSDP in packed symmetric format.
\param n number of rows and columns of the matrix
\param val multiple of identity matrix.
\param dops address of a pointer to a table of function pointers
\param imat address of a pointer to an opaque data type.
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetIdentityP"
int DSDPGetIdentityDataMatP(int n, double val, struct  DSDPDataMat_Ops** dops, void** imat){
  int info;
  identitymat *AA;
  
  DSDPFunctionBegin;
  AA=(identitymat*) malloc(1*sizeof(identitymat));
  AA->dm=val;
  AA->n=n;
  info=IdentitymatOperationsInitializeP(&identitymatopsp); DSDPCHKERR(info);
  if (dops){*dops=&identitymatopsp;}
  if (imat){*imat=(void*)AA;}
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetIdentityDataMatF(int n, double val, struct  DSDPDataMat_Ops** dops, void** imat)

\brief Create a sparse matrix usuable by DSDP in full symmetric format.
\param n number of rows and columns of the matrix
\param val multiple of identity matrix.
\param dops address of a pointer to a table of function pointers
\param imat address of a pointer to an opaque data type.
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetIdentityF"
int DSDPGetIdentityDataMatF(int n, double val, struct  DSDPDataMat_Ops** dops, void** imat){
  int info;
  identitymat *AA;
  
  DSDPFunctionBegin;
  AA=(identitymat*) malloc(1*sizeof(identitymat));
  AA->dm=val;
  AA->n=n;
  info=IdentitymatOperationsInitializeF(&identitymatopsf); DSDPCHKERR(info);
  if (dops){*dops=&identitymatopsf;}
  if (imat){*imat=(void*)AA;}
  DSDPFunctionReturn(0);
}

static int IdentityMatDestroy(void* AA){
  free(AA);
  return 0;
}


static int IdentityMatVecVec(void* AA, double x[], int n, double *v){
  identitymat* A=(identitymat*)AA;
  int i;
  *v=0;
  for (i=0;i<n;i++){
    *v+=x[i]*x[i];
  }
  *v *= A->dm;
  return 0;
}

static int IdentityMatDotP(void* AA, double x[], int nn, int n, double *v){
  identitymat* A=(identitymat*)AA;
  int i;
  double *xx=x;
  *v=0;
  for (i=0;i<n;i++){
    *v+=*xx;
    xx+=i+2;
  }
  *v *= 2*A->dm;
  return 0;
}

static int IdentityMatDotF(void* AA, double x[], int nn, int n, double *v){
  identitymat* A=(identitymat*)AA;
  int i;
  double *xx=x;
  *v=0;
  for (i=0;i<n;i++){
    *v+=*xx;
    xx+=n+1;
  }
  *v *= 2*A->dm;
  return 0;
}

static int IdentityMatFNorm2(void* AA, int n, double *v){
  identitymat* A=(identitymat*)AA;
  *v=A->n*A->dm*A->dm;
  return 0;
}

static int IdentityMatView(void* AA){
  identitymat* A=(identitymat*)AA;
  printf("Multiple of Identity matrix: All Diagonal elements equal %8.8e \n",A->dm);
  return 0;
}

static int IdentityMatGetRank(void *AA, int*rank, int n){
  identitymat* A=(identitymat*)AA;
  *rank=A->n;;
  return 0;
}

static int IdentityMatFactor(void*A){
  return 0;
}

static int IdentityMatGetEig(void*AA, int neig, double *eig, double v[], int n, int* indx, int *nind){
  identitymat* A = (identitymat*)AA;

  if (neig<0 || neig>= A->n){ *eig=0; return 0;} 
  memset((void*)v,0,(A->n)*sizeof(double)); 
  v[neig]=1.0;
  indx[0]=neig;
  *nind=1;
  *eig=A->dm;  
  return 0;
}

static int IdentityMatGetRowNnz(void*A, int nrow, int nz[], int *nnzz, int n){
  identitymat* AA = (identitymat*)A;
  if (nrow>=0 && nrow < AA->n){
    *nnzz=1;
    nz[nrow]++;
  } else {
    *nnzz=0;    
  }
  return 0;
}

static int IdentityMatCountNonzeros(void*A, int *nnz, int n){
  identitymat* AA = (identitymat*)A;
  *nnz=AA->n;
  return 0;
}

static int IdentityMatAddRowMultiple(void*A, int nrow, double dd, double rrow[], int n){
  identitymat* AA = (identitymat*)A;
  rrow[nrow] += dd*AA->dm;
  return 0;
}

static int IdentityMatAddMultipleP(void*A, double dd, double vv[], int nn, int n){
  identitymat* AA = (identitymat*)A;
  double *v=vv,dm=dd*AA->dm;
  int i;
  for (i=0;i<n;i++){
    *v += dm;
    v+= i+2;
  }
  return 0;
}

static int IdentityMatAddMultipleF(void*A, double dd, double vv[], int nn, int n){
  identitymat* AA = (identitymat*)A;
  double *v=vv,dm=dd*AA->dm;
  int i;
  for (i=0;i<n;i++){
    *v += dm;
    v+= n+1;
  }
  return 0;
}

static const char *datamatname="MULTIPLE OF IDENTITY";

static int IdentitymatOperationsInitializeP(struct  DSDPDataMat_Ops* spdiagops){
  int info;
  if (spdiagops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(spdiagops); if (info){return info;}
  spdiagops->matfactor1=IdentityMatFactor;
  spdiagops->matgetrank=IdentityMatGetRank;
  spdiagops->matgeteig=IdentityMatGetEig;
  spdiagops->matvecvec=IdentityMatVecVec;
  spdiagops->matrownz=IdentityMatGetRowNnz;
  spdiagops->matdot=IdentityMatDotP;
  spdiagops->matfnorm2=IdentityMatFNorm2;
  spdiagops->matnnz=IdentityMatCountNonzeros;
  spdiagops->mataddrowmultiple=IdentityMatAddRowMultiple;
  spdiagops->mataddallmultiple=IdentityMatAddMultipleP;
  spdiagops->matdestroy=IdentityMatDestroy;
  spdiagops->matview=IdentityMatView;
  spdiagops->id=12;
  spdiagops->matname=datamatname;
  return 0;
}

static int IdentitymatOperationsInitializeF(struct  DSDPDataMat_Ops* spdiagops){
  int info;
  if (spdiagops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(spdiagops); if (info){return info;}
  spdiagops->matfactor1=IdentityMatFactor;
  spdiagops->matgetrank=IdentityMatGetRank;
  spdiagops->matgeteig=IdentityMatGetEig;
  spdiagops->matvecvec=IdentityMatVecVec;
  spdiagops->matrownz=IdentityMatGetRowNnz;
  spdiagops->matdot=IdentityMatDotF;
  spdiagops->matfnorm2=IdentityMatFNorm2;
  spdiagops->matnnz=IdentityMatCountNonzeros;
  spdiagops->mataddrowmultiple=IdentityMatAddRowMultiple;
  spdiagops->mataddallmultiple=IdentityMatAddMultipleF;
  spdiagops->matdestroy=IdentityMatDestroy;
  spdiagops->matview=IdentityMatView;
  spdiagops->id=12;
  spdiagops->matname=datamatname;
  return 0;
}
