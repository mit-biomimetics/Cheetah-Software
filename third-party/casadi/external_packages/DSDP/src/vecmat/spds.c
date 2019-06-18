#include "dsdpsys.h"
#include "dsdpdsmat_impl.h"

/*! \file spds.c
\brief DSDPDualMat object with sparse data structures. 
*/

typedef struct {
  int    n;
  double *an;
  int    *col;
  int    *nnz;
} spdsmat;

static int SpSymMatSetURValuesP(void*DS, double v[], int nn, int n){
  spdsmat*ds=(spdsmat*)DS;
  int i,j,k1,k2,*nnz=ds->nnz,*col=ds->col;
  double *an=ds->an;
  for (i=0;i<n;i++,nnz++){
    k1=*nnz; k2=*(nnz+1);
    for (j=k1;j<k2;j++,an++,col++){
      if ((*col)==i){ *an = v[*col]/2;}
      else { *an = v[*col]; }
    }
    v+=i+1;
  }
  return 0;
}

static int SpSymMatSetURValuesU(void*DS, double v[], int nn, int n){
  spdsmat*ds=(spdsmat*)DS;
  int i,j,k1,k2,*nnz=ds->nnz,*col=ds->col;
  double *an=ds->an;
  for (i=0;i<n;i++,nnz++){
    k1=*nnz; k2=*(nnz+1);
    for (j=k1;j<k2;j++,an++,col++){
      if ((*col)==i){ *an = v[*col]/2;}
      else { *an = v[*col]; }
    }
    v+=n;
  }
  return 0;
}

static int SpSymMatView(void *DS){
  spdsmat*ds=(spdsmat*)DS;
  int i,j,k1,k2,n=ds->n,*nnz=ds->nnz,*col=ds->col;
  double *an=ds->an;
  for (i=0;i<n;i++){
    k1=nnz[i]; k2=nnz[i+1];
    printf("Row %d: ",i);
    for (j=k1;j<k2;j++){
      if (col[j]==i){ printf("%d: %4.4f",col[j],2*an[j]); }
      else { printf("%d: %4.4f",col[j],an[j]);}
    }
    printf("\n");
  }
  return 0;
}
/*
static int SpSymMatShiftDiagonal(void *DS, double dd){
  spdsmat*ds=(spdsmat*)DS;
  int i,n=ds->n,*nnz=ds->nnz;
  double *an=ds->an;
  for (i=0;i<n;i++){
    an[nnz[i+1]-1] += dd/2;
  }
  return 0;
}
*/
static int SpSymMatDestroy(void *DS){
  spdsmat*ds=(spdsmat*)DS;
  int info;
  DSDPFREE(&ds->nnz,&info);if (info) return 1;
  DSDPFREE(&ds->col,&info);if (info) return 1;
  DSDPFREE(&ds->an,&info);if (info) return 1;
  DSDPFREE(&ds,&info);if (info) return 1;
  return 0;
}

static int SpSymMatGetSize(void *DS, int*n){
  spdsmat*ds=(spdsmat*)DS;
  *n=ds->n;
  return 0;
}

static int SpSymMatZero(void*DS){
  spdsmat*ds=(spdsmat*)DS;
  int nn=ds->nnz[ds->n];
  double *an=ds->an;
  memset((void*)an,0,nn*sizeof(double));
  return 0;
}

static int SpSymMatMult(void*DS, double x[], double y[], int n){
  spdsmat*ds=(spdsmat*)DS;
  int i,j,k1,k2,*nnz=ds->nnz,*col=ds->col;
  double *an=ds->an;
  memset((void*)y,0,n*sizeof(double));
  for (i=0;i<n;i++,nnz++){
    k1=*nnz; k2=*(nnz+1);
    for (j=k1;j<k2;j++,col++,an++){
      y[*col] += x[i] * (*an);
      y[i] += x[*col] * (*an);
    }
  }
  return 0;
}

static int SpSymMatVecVec(void*DS, double x[], int n, double *vAv){
  spdsmat*ds=(spdsmat*)DS;
  int i,j,k1,k2,*nnz=ds->nnz,*col=ds->col;
  double vv,*an=ds->an;
  *vAv=0;
  for (i=0;i<n;i++,nnz++){
    k1=*nnz; k2=*(nnz+1);
    vv=0;
    for (j=k1;j<k2;j++,col++,an++){
      vv+=x[*col]*(*an);
    }
    *vAv+=vv*x[i]*2;
  }
  return 0;
}
/*
static int SpSymMatAddRow(void *DS, int row, double dd, double v[], int n){
  spdsmat*ds=(spdsmat*)DS;
  int j,k1,k2,*nnz=ds->nnz,*col=ds->col;
  double *an=ds->an;
  k1=nnz[row]; k2=nnz[row+1];
  for (j=k1;j<k2;j++){
    if (row==col[j]){ an[j] += dd*v[col[j]]/2; } 
    else { an[j] += dd*v[col[j]]; }
  }
  return 0;
}
*/
static const char* dsmatname="SPARSE, SYMMETRIC MATRIX";
static int DSDPDSSparseInitializeOpsP(struct  DSDPDSMat_Ops* dsops){
  int info;
  if (!dsops) return 0;
  info=DSDPDSMatOpsInitialize(dsops); DSDPCHKERR(info);
  dsops->matseturmat=SpSymMatSetURValuesP;
  dsops->matview=SpSymMatView;
  dsops->matdestroy=SpSymMatDestroy;
  dsops->matgetsize=SpSymMatGetSize;
  dsops->matzeroentries=SpSymMatZero;
  dsops->matmult=SpSymMatMult;
  dsops->matvecvec=SpSymMatVecVec;
  dsops->id=6;
  dsops->matname=dsmatname;
  return 0;
}
static int DSDPDSSparseInitializeOpsU(struct  DSDPDSMat_Ops* dsops){
  int info;
  if (!dsops) return 0;
  info=DSDPDSMatOpsInitialize(dsops); DSDPCHKERR(info);
  dsops->matseturmat=SpSymMatSetURValuesU;
  dsops->matview=SpSymMatView;
  dsops->matdestroy=SpSymMatDestroy;
  dsops->matgetsize=SpSymMatGetSize;
  dsops->matzeroentries=SpSymMatZero;
  dsops->matmult=SpSymMatMult;
  dsops->matvecvec=SpSymMatVecVec;
  dsops->id=6;
  dsops->matname=dsmatname;
  return 0;
}

static struct  DSDPDSMat_Ops tdsdsopsp;
static struct  DSDPDSMat_Ops tdsdsopsu;
#undef __FUNCT__
#define __FUNCT__ "DSDPCreateSparseDSMat"
int DSDPSparseMatCreatePattern2P(int n, int rnnz[], int cols[], int tnnz,struct  DSDPDSMat_Ops* *dsmatops, void**dsmat){
  int i,info;
  spdsmat*ds;
  DSDPFunctionBegin;
  DSDPCALLOC1(&ds,spdsmat,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&ds->nnz,int,(n+1),&info);DSDPCHKERR(info);
  ds->nnz[0]=0;
  for (i=0;i<n;i++) ds->nnz[i+1]=ds->nnz[i]+rnnz[i];
  DSDPCALLOC2(&ds->col,int,tnnz,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&ds->an,double,tnnz,&info);DSDPCHKERR(info);
  for (i=0;i<tnnz;i++) ds->col[i]=cols[i];
  info=DSDPDSSparseInitializeOpsP(&tdsdsopsp); DSDPCHKERR(info);
  *dsmatops=&tdsdsopsp;
  *dsmat=(void*)ds;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCreateSparseDSMatU"
int DSDPSparseMatCreatePattern2U(int n, int rnnz[], int cols[], int tnnz,struct  DSDPDSMat_Ops* *dsmatops, void**dsmat){
  int i,info;
  spdsmat*ds;
  DSDPFunctionBegin;
  DSDPCALLOC1(&ds,spdsmat,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&ds->nnz,int,(n+1),&info);DSDPCHKERR(info);
  ds->nnz[0]=0;
  for (i=0;i<n;i++) ds->nnz[i+1]=ds->nnz[i]+rnnz[i];
  DSDPCALLOC2(&ds->col,int,tnnz,&info);DSDPCHKERR(info);
  DSDPCALLOC2(&ds->an,double,tnnz,&info);DSDPCHKERR(info);
  for (i=0;i<tnnz;i++) ds->col[i]=cols[i];
  info=DSDPDSSparseInitializeOpsU(&tdsdsopsu); DSDPCHKERR(info);
  *dsmatops=&tdsdsopsu;
  *dsmat=(void*)ds;
  DSDPFunctionReturn(0);
}
