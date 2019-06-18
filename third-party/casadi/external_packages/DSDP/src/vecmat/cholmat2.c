#include "numchol.h"
#include "dsdpdualmat_impl.h"
#include "dsdpsys.h"
#include "dsdplapack.h"

/*!
\file cholmat2.c 
\brief Sparse Cholesky for dual S matrix
*/

typedef struct{
  chfac* spsym;
  double *sinv;
  char   UPLQ;
  int n;
  int dsinv;
} spmat;


static int SMatDestroy(void*S){
  spmat* SS=(spmat*)S;
  int info;
  CfcFree(&SS->spsym);
  if (SS->dsinv){
    DSDPFREE(&SS->sinv,&info);
  }
  DSDPFREE(&SS,&info);
  return 0;
}

static int SMatGetSize(void *S, int *n){
  spmat* SS=(spmat*)S;
  *n=SS->spsym->nrow;
  return 0;
}

static int SMatView(void* S){
  spmat* SS=(spmat*)S;
  int info;
  info=Mat4View(SS->spsym); DSDPCHKERR(info);
  return 0;
}

static int SMatLogDet(void* S, double *dd){
  spmat* SS=(spmat*)S;
  int info;
  info=Mat4LogDet(SS->spsym,dd); DSDPCHKERR(info);
  return 0;
}



static int SMatSetURMatP(void*S, double v[], int nn, int n){
  spmat* SS=(spmat*)S;
  int k,j,row,info;
  double *rw1,*rw2,*xr=v;
  rw1=SS->spsym->rw;rw2=rw1+n;
  info=MatZeroEntries4(SS->spsym); DSDPCHKERR(info);
  for (k=0;k<n/2;k++){
    row = 2*k;

    xr=v+row*(row+1)/2;
    memcpy((void*)rw1,(void*)(xr),(row+1)*sizeof(double));
    xr+=row+1;
    rw1[row+1]=xr[row];
    memcpy((void*)(rw2),(void*)(xr),(row+2)*sizeof(double));
    xr+=row+2;

    /*  memset((void*)rw,0,(n-row)*sizeof(double)); */
    for (j=row+2;j<n;j++){ 
      rw1[j]=xr[row];
      rw2[j]=xr[row+1];
      xr+=j+1; 
    }

    info=MatSetColumn4(SS->spsym,rw1,row);  DSDPCHKERR(info);
    info=MatSetColumn4(SS->spsym,rw2,row+1);  DSDPCHKERR(info);
  }

  for (row=2*(n/2);row<n;row++){
    /*  memset((void*)rw,0,(n-row)*sizeof(double)); */
    xr=v+row*(row+1)/2;
    memcpy((void*)(rw1),(void*)(xr),(row+1)*sizeof(double));
    xr+=row+1;
    for (j=row+1;j<n;j++){ rw1[j]=xr[row]; xr+=(j+2); }
    info=MatSetColumn4(SS->spsym,rw1,row);  DSDPCHKERR(info);
    xr+=(n-row);
  }
  return 0;
}

static int SMatSetURMatU(void*S, double v[], int nn, int n){
  spmat* SS=(spmat*)S;
  int k,j,row,info;
  double *rw1,*rw2,*xr=v;
  rw1=SS->spsym->rw;rw2=rw1+n;
  info=MatZeroEntries4(SS->spsym); DSDPCHKERR(info);
  for (k=0;k<n/2;k++){
    row = 2*k;

    xr=v+row*n;
    memcpy((void*)rw1,(void*)(xr),(row+1)*sizeof(double));
    xr+=n;
    rw1[row+1]=xr[row];
    memcpy((void*)(rw2),(void*)(xr),(row+2)*sizeof(double));
    xr+=n;
    /*  memset((void*)rw,0,(n-row)*sizeof(double)); */
    for (j=row+2;j<n;j++){ 
      rw1[j]=xr[row];
      rw2[j]=xr[row+1];
      xr+=n; 
    }

    info=MatSetColumn4(SS->spsym,rw1,row);  DSDPCHKERR(info);
    info=MatSetColumn4(SS->spsym,rw2,row+1);  DSDPCHKERR(info);
  }

  for (row=2*(n/2);row<n;row++){
    /*  memset((void*)rw,0,(n-row)*sizeof(double)); */
    xr=v+row*n;
    memcpy((void*)(rw1),(void*)(xr),(row+1)*sizeof(double));
    xr+=n;
    for (j=row+1;j<n;j++){ rw1[j]=xr[row]; xr+=n; }
    info=MatSetColumn4(SS->spsym,rw1,row);  DSDPCHKERR(info);
  }
  return 0;
}

static int SMatSetURMat(void*S, double v[], int nn, int n){
  spmat* SS=(spmat*)S;
  int info;
  if (SS->UPLQ=='P'){
    info=SMatSetURMatP(S,v,nn,n);DSDPCHKERR(info);
  } else if (SS->UPLQ=='U'){
    info=SMatSetURMatU(S,v,nn,n);DSDPCHKERR(info);
  }
  return 0;
}

static int SMatSolve(void *S, int indx[], int nind, double b[], double x[],int n){
  spmat* SS=(spmat*)S;
  int i,ii;
  double alpha,*s1=SS->sinv,*s2;
  ffinteger nn,ione;
  if (SS->sinv && nind < n/4){
    memset((void*)x,0,n*sizeof(double));
    for (i=0;i<nind;i++){
      ii=indx[i];
      ione=1;nn=n;alpha=b[ii];s2=s1+n*ii;
      daxpy(&nn,&alpha,s2,&ione,x,&ione);
    }
  } else {
    memcpy(x,b,n*sizeof(double));
    ChlSolve(SS->spsym, b, x);
  }
  return 0;
}

static int SMatCholeskySolveBackward(void *S, double b[], double x[],int n){
  spmat* SS=(spmat*)S;
  ChlSolveBackward(SS->spsym, b, x);
  return 0;
}

static int SMatCholeskySolveForward(void *S, double b[], double x[], int n){
  spmat* SS=(spmat*)S;
  ChlSolveForward(SS->spsym, b, x);
  return 0;
}

static int SMatFull(void *S, int *full){
  *full=0;
  return 0;
}

static int SMatCholeskyFactor(void *S, int *flag){
  spmat* SS=(spmat*)S;
  int *iw;
  double *rw;
  cfc_sta Cfact;
  iw=SS->spsym->iw;
  rw=SS->spsym->rw;
  Cfact=(cfc_sta)ChlFact(SS->spsym,iw,rw,TRUE);
  if (CfcOk!= Cfact){ 
    *flag=1;
  } else {
    *flag=0;
  }
  return 0;
}

static int SMatInverseAddP(void *S, double alpha, double v[], int nn, int n){
  spmat* SS=(spmat*)S;
  ffinteger ii,ione=1;
  double *x,*b,*ss=SS->sinv;
  int i,j,k=0;

  if (ss){
    for (i=0;i<n;i++){
      v+=i; ii=i+1;
      daxpy(&ii,&alpha,ss,&ione,v,&ione);
      ss+=n;
    }
  } else {
    b=SS->spsym->rw;x=b+n;
    for (i=0;i<n;i++){
      memset((void*)b,0,n*sizeof(double));
      b[i]=alpha;
      ChlSolve(SS->spsym, b, x);
      k=k+i;
      for (j=0;j<=i;j++){
	v[k+j]+=x[j];
      }
    }
  }
  return 0;
}

static int SMatInverseAddU(void *S, double alpha, double v[], int nn, int n){
  spmat* SS=(spmat*)S;
  ffinteger n2=n*n,ione=1;
  double *x,*b,*ss=SS->sinv;
  int i,j,k=0;

  if (ss){
    daxpy(&n2,&alpha,ss,&ione,v,&ione);
  } else {
    b=SS->spsym->rw;x=b+n;
    for (i=0;i<n;i++){
      memset((void*)b,0,n*sizeof(double));
      b[i]=alpha;
      ChlSolve(SS->spsym, b, x);
      k=i*n;
      for (j=0;j<n;j++){
	v[k+j]+=x[j];
      }
    }
  }
  return 0;
}

static int SMatInverseAdd(void *S, double alpha, double v[], int nn, int n){
  spmat* SS=(spmat*)S;
  int info;
  if (SS->UPLQ=='P'){
    info=SMatInverseAddP(S,alpha,v,nn,n);DSDPCHKERR(info);
  } else if (SS->UPLQ=='U'){
    info=SMatInverseAddU(S,alpha,v,nn,n);DSDPCHKERR(info);
  }
  return 0;
}

static int SMatCholeskyForwardMultiply(void *S, double b[], double x[], int n){
  spmat* SS=(spmat*)S;
  GetUhat(SS->spsym,b,x);
  return 0;
}

static int SMatInvert(void *S){
  spmat* SS=(spmat*)S;
  double *x,*b,*v=SS->sinv;
  int i,n=SS->n;

  b=SS->spsym->rw;x=b+n;

  if (v){
    for (i=0;i<n;i++){
      memset((void*)b,0,n*sizeof(double));
      b[i]=1.0;
      ChlSolve(SS->spsym, b, x);
      memcpy((void*)(v+i*n),(void*)x,n*sizeof(double));
    }
  }
  return 0;
}

static struct  DSDPDualMat_Ops sdmatops;
static const char* tmatname="SPARSE PSD";
static int SDualOpsInitialize(struct  DSDPDualMat_Ops* sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDualMatOpsInitialize(sops); DSDPCHKERR(info);
  sops->matcholesky=SMatCholeskyFactor;
  sops->matsolveforward=SMatCholeskySolveForward;
  sops->matsolvebackward=SMatCholeskySolveBackward;
  sops->matinversemultiply=SMatSolve;
  sops->matinvert=SMatInvert;
  sops->matinverseadd=SMatInverseAdd;
  sops->matforwardmultiply=SMatCholeskyForwardMultiply;
  sops->matseturmat=SMatSetURMat;
  sops->matfull=SMatFull;
  sops->matdestroy=SMatDestroy;
  sops->matgetsize=SMatGetSize;
  sops->matview=SMatView;
  sops->matlogdet=SMatLogDet;
  sops->matname=tmatname;
  return 0;
 }

static int dcholmatcreate(int n, char UPLQ, chfac *sp, 
			   struct  DSDPDualMat_Ops **sops, void**smat){
  spmat  *S;
  int info;
  DSDPCALLOC1(&S,spmat,&info);DSDPCHKERR(info);
  S->UPLQ=UPLQ; S->n=n; S->sinv=0; S->dsinv=0; S->spsym=sp;
  info=SDualOpsInitialize(&sdmatops);DSDPCHKERR(info);
  *sops=&sdmatops;
  *smat=(void*)S;
  return 0;
 }

static int dcholmatsinverse(int n, spmat *S1, spmat *S2){
  int info;
  double *ssinv;
  DSDPCALLOC2(&ssinv,double,n*n,&info);
  S1->sinv=ssinv; S2->sinv=ssinv; S2->dsinv=1;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDenseDualMatCreate"
int DSDPDenseDualMatCreate(int n, char UPLQ, 
			   struct  DSDPDualMat_Ops **sops1, void**smat1,
			   struct  DSDPDualMat_Ops **sops2, void**smat2){
  int info=0;
  chfac *sp;

  DSDPFunctionBegin;
  info=MchlSetup2(n,&sp); DSDPCHKERR(info);
  info=dcholmatcreate(n,UPLQ,sp,sops1,smat1);DSDPCHKERR(info);
  info=MchlSetup2(n,&sp); DSDPCHKERR(info);
  info=dcholmatcreate(n,UPLQ,sp,sops1,smat2);DSDPCHKERR(info);
  info=dcholmatsinverse(n,(spmat*)(*smat1),(spmat*)(*smat2));DSDPCHKERR(info);

  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPSparseDualMatCreate"
int DSDPSparseDualMatCreate(int n, int *rnnz, int *snnz, 
			    int trank,char UPLQ,int*nnzz, 
			    struct  DSDPDualMat_Ops **sops1, void**smat1,
			    struct  DSDPDualMat_Ops **sops2, void**smat2){
  int nnz,info=0;
  chfac *sp;

  DSDPFunctionBegin;
  SymbProc(rnnz,snnz,n,&sp); DSDPCHKERR(info);
  info=dcholmatcreate(n,UPLQ,sp,sops1,smat1);DSDPCHKERR(info);
  SymbProc(rnnz,snnz,n,&sp); DSDPCHKERR(info);
  info=dcholmatcreate(n,UPLQ,sp,sops2,smat2);DSDPCHKERR(info);
  nnz=sp->unnz;*nnzz=nnz;
  if (trank>2*n+2){
    info=dcholmatsinverse(n,(spmat*)(*smat1),(spmat*)(*smat2));DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}
