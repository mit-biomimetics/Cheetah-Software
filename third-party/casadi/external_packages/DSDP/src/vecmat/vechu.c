#include "dsdpdatamat_impl.h"
#include "dsdpsys.h"
/*! \file vechu.c
\brief DSDPDataMat for sparse matrices in upper full symmetric format 
*/

typedef struct {
  int    neigs;
  double *eigval;
  double *an;
  int    *cols,*nnz;
} Eigen;

typedef struct {
  int    nnzeros;
  const int    *ind;
  const double *val;
  int ishift;
  double alpha;

  Eigen   *Eig;
  int factored;
  int owndata;
  int    n;
} vechmat;

#define GETI(a,b)    (int)((int)a/(int)b)
#define GETJ(a,b)    (int)((int)a%(int)b)

static void getij(int k, int n, int *i,int *j){
  *i=GETI(k,n);
  *j=GETJ(k,n);
  return;
}

#undef __FUNCT__  
#define __FUNCT__ "CreateVechMatWData"
static int CreateVechMatWdata(int n, int ishift, double alpha, const int *ind, const double *vals, int nnz, vechmat **A){
  int info;
  vechmat* V;
  DSDPCALLOC1(&V,vechmat,&info);DSDPCHKERR(info);
  V->n=n; V->ishift=ishift, V->ind=ind; V->val=vals;V->nnzeros=nnz;
  V->alpha=alpha;
  V->owndata=0;
  *A=V;
  return 0;
}

static int VechMatAddMultiple(void* AA, double scl, double r[], int nn, int n){
  vechmat* A=(vechmat*)AA;
  int k;
  const int *ind=A->ind,nnz=A->nnzeros;
  const double *val=A->val;
  double *rr=r-A->ishift;
  scl*=A->alpha;
  for (k=0; k<nnz; ++k){
    *(rr+(ind[k])) +=scl*(val[k]);
  }
  return 0;
}

static int VechMatDot(void* AA, double x[], int nn, int n, double *v){
  vechmat* A=(vechmat*)AA;
  int k,nnz=A->nnzeros; 
  const int *ind=A->ind;
  double vv=0,*xx=x-A->ishift;
  const double *val=A->val;
  for (k=0;k<nnz;++k,++ind,++val){
    vv+=(*val)*(*(xx+(*ind)));
  }
  *v=2*vv*A->alpha;
  return 0;
}

static int EigMatVecVec(Eigen*, double[], int, double*);
static int VechMatGetRank(void*,int*,int);

static int VechMatVecVec(void* AA, double x[], int n, double *v){
  vechmat* A=(vechmat*)AA;
  int info,rank=n,i=0,j,k,kk;
  const int *ind=A->ind,ishift=A->ishift,nnz=A->nnzeros;
  double vv=0,dd;
  const double *val=A->val;
  
  if (A->factored==3){
    info=VechMatGetRank(AA,&rank,n);
    if (nnz>3 && rank<nnz){
      info=EigMatVecVec(A->Eig,x,n,&vv);
      *v=vv*A->alpha;
      return 0;
    }
  }
  
  for (k=0; k<nnz; ++k,++ind,++val){
    kk=*ind-ishift;
    i=GETI(kk,n);
    j=GETJ(kk,n);
    dd=x[i]*x[j]*(*val);
    vv+=2*dd;
    if (i==j){ vv-=dd; }
  }
  *v=vv*A->alpha;

  return 0;
}


static int VechMatGetRowNnz(void* AA, int trow, int nz[], int *nnzz,int nn){
  vechmat* A=(vechmat*)AA;
  int i=0,j,k,t;
  const int *ind=A->ind, ishift=A->ishift,nnz=A->nnzeros;
  *nnzz=0;
  for (k=0; k<nnz; ++k /*,ind*/){ /* Commented out by @jaeandersson */
    t=ind[k]-ishift;
    i=GETI(t,nn);
    j=GETJ(t,nn);
    if (i==trow){
      nz[j]++;(*nnzz)++;
    } else if (j==trow){
      nz[i]++;(*nnzz)++;
    }
  }
  return 0;
}

static int VechMatFNorm2(void* AA, int n, double *fnorm2){
  vechmat* A=(vechmat*)AA;
  int i=0,j,k,t;
  const int *ind=A->ind,ishift=A->ishift,nnz=A->nnzeros;
  double fn2=0;
  const double *val=A->val;
  for (k=0; k<nnz; ++k){
    t=ind[k]-ishift;
    i=GETI(t,n);
    j=GETJ(t,n);
    if (i==j){
      fn2+= val[k]*val[k];
    } else {
      fn2+= 2*val[k]*val[k];
    }
  }
  *fnorm2=fn2*A->alpha*A->alpha;
  return 0;
}

static int VechMatAddRowMultiple(void* AA, int trow, double scl, double r[], int m){
  vechmat* A=(vechmat*)AA;
  int i=0,j,k,t,ishift=A->ishift,nnz=A->nnzeros;
  const int *ind=A->ind;
  const double *val=A->val;
  scl*=A->alpha;
  for (k=0; k<nnz; ++k){
    t=ind[k]-ishift;
    i=GETI(t,m);
    j=GETJ(t,m);
    if (i==trow){
      r[j]+=scl*val[k];
    } else if (j==trow){
      r[i]+=scl*val[k];
    }
  }
  return 0;
}

static int VechMatCountNonzeros(void* AA, int*nnz, int n){
  vechmat* A=(vechmat*)AA;
  *nnz=A->nnzeros;
  return 0;
}

#undef __FUNCT__  
#define __FUNCT__ "VechMatDestroy"
static int VechMatDestroy(void* AA){
  vechmat* A=(vechmat*)AA;
  int info;
  if (A->owndata){
    /*
    if (A->ind){ DSDPFREE(&A->ind,&info);DSDPCHKERR(info);}
    if (A->val){ DSDPFREE(&A->val,&info);DSDPCHKERR(info);}
    */
    return 1;
  }
  if (A->Eig){
    DSDPFREE(&A->Eig->eigval,&info);DSDPCHKERR(info);
    DSDPFREE(&A->Eig->an,&info);DSDPCHKERR(info);
    if (A->Eig->cols){DSDPFREE(&A->Eig->cols,&info);DSDPCHKERR(info);}
    if (A->Eig->nnz){DSDPFREE(&A->Eig->nnz,&info);DSDPCHKERR(info);}
    DSDPFREE(&A->Eig,&info);DSDPCHKERR(info);
  }
  DSDPFREE(&A,&info);DSDPCHKERR(info);
  return 0;
}



#undef __FUNCT__  
#define __FUNCT__ "DSDPCreateVechMatEigs"
static int CreateEigenLocker(Eigen **EE,int iptr[], int neigs, int n){
  int i,k,info;
  Eigen *E;

  for (k=0,i=0;i<neigs;i++) k+=iptr[i];
  if (k>n*neigs/4){

    DSDPCALLOC1(&E,Eigen,&info);DSDPCHKERR(info);
    DSDPCALLOC2(&E->eigval,double,neigs,&info);DSDPCHKERR(info);
    DSDPCALLOC2(&E->an,double,n*neigs,&info);DSDPCHKERR(info);
    E->neigs=neigs;
    E->cols=0;
    E->nnz=0;

  } else {

    DSDPCALLOC1(&E,Eigen,&info);DSDPCHKERR(info);
    DSDPCALLOC2(&E->eigval,double,neigs,&info);DSDPCHKERR(info);
    DSDPCALLOC2(&E->nnz,int,neigs,&info);DSDPCHKERR(info);
    DSDPCALLOC2(&E->an,double,k,&info);DSDPCHKERR(info);
    DSDPCALLOC2(&E->cols,int,k,&info);DSDPCHKERR(info);
    E->neigs=neigs;

    if (neigs>0) E->nnz[0]=iptr[0];
    for (i=1;i<neigs;i++){E->nnz[i]=E->nnz[i-1]+iptr[i];}
  }
  *EE=E;
  return 0;
}


static int EigMatSetEig(Eigen* A,int row, double eigv, int idxn[], double v[], int nsub,int n){
  int j,k,*cols=A->cols;
  double *an=A->an;
  A->eigval[row]=eigv;
  if (cols){
    k=0; if (row>0){ k=A->nnz[row-1];}
    cols+=k; an+=k;
    for (k=0,j=0; j<nsub; j++){
      if (v[j]==0.0) continue;
      cols[k]=idxn[j]; an[k]=v[j]; k++;
    }
  } else {
    an+=n*row;
    for (j=0; j<nsub; j++){
      if (v[j]==0.0) continue;
      an[idxn[j]]=v[j];
    }
  }
  return 0;
}


static int EigMatGetEig(Eigen* A,int row, double *eigenvalue, double eigenvector[], int n, int spind[], int *nind){
  int i,*cols=A->cols,bb,ee;
  double* an=A->an;
  *eigenvalue=A->eigval[row];
  *nind=0;
  if (cols){
    memset((void*)eigenvector,0,n*sizeof(double));
    if (row==0){ bb=0;} else {bb=A->nnz[row-1];} ee=A->nnz[row];
    for (i=bb;i<ee;i++){
      eigenvector[cols[i]]=an[i];
      spind[i-bb]=cols[i]; (*nind)++;
    }
  } else {
    memcpy((void*)eigenvector,(void*)(an+n*row),n*sizeof(double));
    for (i=0;i<n;i++)spind[i]=i;
    *nind=n;
  }
  return 0;
}

static int EigMatVecVec(Eigen* A, double v[], int n, double *vv){
  int i,rank,*cols=A->cols,neigs=A->neigs,*nnz=A->nnz,bb,ee;
  double* an=A->an,*eigval=A->eigval,dd,ddd=0;

  if (cols){
    for (rank=0;rank<neigs;rank++){
      if (rank==0){ bb=0;} else {bb=nnz[rank-1];} ee=nnz[rank];
      for (dd=0,i=bb;i<ee;i++){
	dd+=an[i]*v[cols[i]];
      }
      ddd+=dd*dd*eigval[rank];
    }
  } else {
    for (rank=0;rank<neigs;rank++){
      for (dd=0,i=0;i<n;i++){
	dd+=an[i]*v[i];
      }
      an+=n;
      ddd+=dd*dd*eigval[rank];
    }
  }
  *vv=ddd;
  return 0;
}


static int VechMatComputeEigs(vechmat*,double[],int,double[],int,double[],int,int[],int,double[],int,double[],int);

static int VechMatFactor(void*AA, double dmatp[], int nn0, double dwork[], int n, double ddwork[], int n1, int iptr[], int n2){

  vechmat*  A=(vechmat*)AA;
  int i,j,k,t,info,is_diag;
  const int *ind=A->ind,ishift=A->ishift,nonzeros=A->nnzeros;
  double *ss1=0,*ss2=0;
  int nn1=0,nn2=0;
  if (A->factored) return 0;

  memset((void*)iptr,0,3*n*sizeof(int));  
  /* Find number of nonzeros in each row */
  for (is_diag=1,k=0; k<nonzeros; k++){
    t=ind[k]-ishift;
    i=GETI(t,n);
    j=GETJ(t,n);
    iptr[i]++;
    if (i!=j) {is_diag=0;iptr[j]++;}
  }
  
  if (is_diag){ A->factored=1; return 0;}
  /* Find most nonzeros per row */
  for (j=0,i=0; i<n; i++){ if (iptr[i]>j) j=iptr[i]; }
  if (j<2){ A->factored=2; return 0; }
  
  info=VechMatComputeEigs(A,dmatp,nn0,dwork,n,ddwork,n1,iptr,n2,ss1,nn1,ss2,nn2);DSDPCHKERR(info);
  A->factored=3;
  return 0;
}

static int VechMatGetRank(void *AA,int *rank,int n){
  vechmat*  A=(vechmat*)AA;
  switch (A->factored){
  case 1:
    *rank=A->nnzeros;
    break;
  case 2:
    *rank=2*A->nnzeros;
    break;
  case 3:
    *rank=A->Eig->neigs;
    break;
  default:
    DSDPSETERR(1,"Vech Matrix not factored yet\n");
  }
  return 0;
}

static int VechMatGetEig(void* AA, int rank, double *eigenvalue, double vv[], int n, int indx[], int *nind){
  vechmat*  A=(vechmat*)AA;
  const double *val=A->val,tt=sqrt(0.5);
  int info,i,j,k,t;
  const int *ind=A->ind,ishift=A->ishift;

  *nind=0;
  switch (A->factored){
  case 1:
    memset(vv,0,n*sizeof(double));
    t=ind[rank]-ishift;
    i=GETI(t,n);
    j=GETJ(t,n);
    vv[i]=1.0;
    *eigenvalue=val[rank]*A->alpha;
    *nind=1;
    indx[0]=i;
    break;
  case 2:
    memset(vv,0,n*sizeof(double));
    k=rank/2;
    getij(ind[k]-ishift,n,&i,&j);
    if (i==j){
      if (k*2==rank){ 
	vv[i]=1.0; *eigenvalue=val[k]*A->alpha;
	*nind=1;
	indx[0]=i;
      } else {
	*eigenvalue=0;
      }
    } else {
      if (k*2==rank){
	vv[i]=tt;  vv[j]=tt; *eigenvalue=val[k]*A->alpha;
	*nind=2;
	indx[0]=i; indx[1]=j;
      } else {
	vv[i]=-tt; vv[j]=tt; *eigenvalue=-val[k]*A->alpha;
	*nind=2;
	indx[0]=i; indx[1]=j;
      }
    }
    break;
  case 3:
    info=EigMatGetEig(A->Eig,rank,eigenvalue,vv,n,indx,nind);DSDPCHKERR(info);
    *eigenvalue=*eigenvalue*A->alpha;
    break;
  default:
    DSDPSETERR(1,"Vech Matrix not factored yet\n");
  }

  return 0;  
}

static int VechMatView(void* AA){
  vechmat* A=(vechmat*)AA;
  int info,i=0,j,k,rank=0,ishift=A->ishift,n=A->n,nnz=A->nnzeros;
  const int *ind=A->ind;
  const double *val=A->val;
  for (k=0; k<nnz; k++){
    getij(ind[k]-ishift,n,&i,&j);
    printf("Row: %d, Column: %d, Value: %10.8f \n",i,j,A->alpha*val[k]);
  }
  if (A->factored>0){
    info=VechMatGetRank(AA,&rank,n);DSDPCHKERR(info);
    printf("Detected Rank: %d\n",rank);
  }
  return 0;
}


static struct  DSDPDataMat_Ops vechmatops;
static const char *datamatname="STANDARD VECH MATRIX";

static int VechMatOpsInitialize(struct  DSDPDataMat_Ops *sops){
  int info;
  if (sops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(sops); DSDPCHKERR(info);
  sops->matvecvec=VechMatVecVec;
  sops->matdot=VechMatDot;
  sops->matfnorm2=VechMatFNorm2;
  sops->mataddrowmultiple=VechMatAddRowMultiple;
  sops->mataddallmultiple=VechMatAddMultiple;
  sops->matview=VechMatView;
  sops->matdestroy=VechMatDestroy;
  sops->matfactor2=VechMatFactor;
  sops->matgetrank=VechMatGetRank;
  sops->matgeteig=VechMatGetEig;
  sops->matrownz=VechMatGetRowNnz;
  sops->matnnz=VechMatCountNonzeros;
  sops->id=3;
  sops->matname=datamatname;
  return 0;
}

/*!
\fn  int DSDPGetVecUMat(int n,int ishift,double alpha, const int ind[], const double val[],int nnz, struct  DSDPDataMat_Ops**sops, void**smat)
\brief Given data in full symmetric format, create a sparse matrix usuable by DSDP.
\param n number of rows and columns of the matrix
\param ishift the index of the first element in the matrix (usually 0)
\param alpha the multiple of these matrix.
\param ind array of matrix indices.
\param val array of matrix values.
\param nnz number of elements in array.
\param sops address of a pointer to a table of function pointers
\param smat address of a pointer to an opaque data type.
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetVecUMat"
int DSDPGetVecUMat(int n,int ishift,double alpha,const int ind[], const double val[],int nnz, struct  DSDPDataMat_Ops**sops, void**smat){ 
  int info,i,j,k,itmp,nn=n*n;
  double dtmp;
  vechmat* AA;
  DSDPFunctionBegin;
  for (k=0;k<nnz;++k){ 
    itmp=ind[k]-ishift; 
    if (itmp>=nn){
      getij(itmp,n,&i,&j);
      /*
      DSDPSETERR(2,"Illegal index value: Element %d in array has row %d (>0) or column %d (>0) is greater than %d. \n",k+1,i+1,j+1,n);
      */
      DSDPSETERR3(2,"Illegal index value: Element %d in array has index %d greater than or equal to %d. \n",k,itmp,nn);
    } else if (itmp<0){
      DSDPSETERR1(2,"Illegal index value: %d.  Must be >= 0\n",itmp);
    }
  }
  for (k=0;k<nnz;++k) dtmp=val[k];
  info=CreateVechMatWdata(n,ishift,alpha,ind,val,nnz,&AA); DSDPCHKERR(info);
  AA->factored=0;
  AA->Eig=0;
  info=VechMatOpsInitialize(&vechmatops); DSDPCHKERR(info);
  if (sops){*sops=&vechmatops;}
  if (smat){*smat=(void*)AA;}
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "VechMatComputeEigs"
static int VechMatComputeEigs(vechmat* AA,double DD[], int nn0, double W[], int n, double WORK[], int n1, int iiptr[], int n2, double ss1[],int nn1, double ss2[], int nn2){
  
  int i,j,k,nsub,neigs,info,*iptr,*perm,*invp;
  long int *i2darray=(long int*)DD;
  int ownarray1=0,ownarray2=0,ownarray3=0;
  int ishift=AA->ishift,nonzeros=AA->nnzeros;
  const int *ind=AA->ind;
  const double *val=AA->val;
  double *dmatarray=ss1,*dworkarray=ss2,maxeig,eps=1.0e-12,eps2=1.0e-12;
  
  iptr=iiptr;  perm=iptr+n;  invp=perm+n;
  /* These operations were done before calling this routine * /
  / * Integer arrays corresponding to rows with nonzeros and inverse map * /
  memset((void*)iiptr,0,3*n*sizeof(int));

  / * Find number of nonzeros in each row * /
  for (i=0,k=0; k<nonzeros; k++){
    getij(ind[k],i,n,&i,&j);
    iptr[i]++; iptr[j]++;
  } 
  */
  /* Number of rows with a nonzero.  Order the rows with nonzeros. */
  for (nsub=0,i=0; i<n; i++){
    if (iptr[i]>0){ invp[nsub]=i; perm[i]=nsub; nsub++;}
  }
  
  /* create a dense array in which to put numbers */
  if (nsub*nsub>nn1){
    DSDPCALLOC2(&dmatarray,double,(nsub*nsub),&info); DSDPCHKERR(info);
    ownarray1=1;
  }
  memset((void*)dmatarray,0,nsub*nsub*sizeof(double));
  if (nsub*nsub>nn2){
    DSDPCALLOC2(&dworkarray,double,(nsub*nsub),&info); DSDPCHKERR(info);
    ownarray2=1;
  }

  if (nsub*nsub*sizeof(long int)>nn0*sizeof(double)){
    DSDPCALLOC2(&i2darray,long int,(nsub*nsub),&info); DSDPCHKERR(info);
    ownarray3=1;
  }
  
  
  for (i=0,k=0; k<nonzeros; k++){
    getij(ind[k]-ishift,n,&i,&j);
    dmatarray[perm[i]*nsub+perm[j]] += val[k];
    if (i!=j){
      dmatarray[perm[j]*nsub+perm[i]] += val[k];
    }
  }
  /* Call LAPACK to compute the eigenvalues */
  memset((void*)W,0,n*sizeof(double));

  info=DSDPGetEigs(dmatarray,nsub,dworkarray,nsub*nsub,i2darray,nsub*nsub,
		   W,nsub,WORK,n1,iiptr+3*n,n2-3*n);
  if (info){
    memset((void*)dmatarray,0,nsub*nsub*sizeof(double));
    for (i=0,k=0; k<nonzeros; k++){
      getij(ind[k]-ishift,n,&i,&j);
      dmatarray[perm[i]*nsub+perm[j]] += val[k];
      if (i!=j){
	dmatarray[perm[j]*nsub+perm[i]] += val[k];
      }
    }
    info=DSDPGetEigs2(dmatarray,nsub,dworkarray,nsub*nsub,i2darray,nsub*nsub,
		      W,nsub,WORK,n1,iiptr+3*n,n2-3*n); DSDPCHKERR(info);
  }
  /*  dsyev_("V","L",&N,dmatarray,&LDA,W,WORK,&LWORK,&INFO); */
  
  for (maxeig=0,i=0;i<nsub;i++){
    if (fabs(W[i])>maxeig){ maxeig=fabs(W[i]); }
  }
  memset((void*)iptr,0,nsub*sizeof(int));
  /* Compute sparsity pattern for  eigenvalue and eigenvector structures */
  /* Count the nonzero eigenvalues */
  for (neigs=0,k=0; k<nsub; k++){
    if (fabs(W[k]) /* /maxeig */ > eps){
      for (j=0;j<nsub;j++){
	if (fabs(dmatarray[nsub*k+j]) >= eps2){iptr[neigs]++;
	} else { dmatarray[nsub*k+j]=0.0;}
      }
      neigs++;
      /*
    } else if (fabs(W[k])>1.0e-100){
      printf("SKIPPING EIGENVALUE: %4.4e, max is : %4.4e\n",W[k],maxeig);
      */
    } 
  }
  
  info=CreateEigenLocker(&AA->Eig,iptr,neigs,n);DSDPCHKERR(info);
  DSDPLogInfo(0,49," Data Mat has %d eigenvectors.",neigs);
  /* Copy into structure */
  for (neigs=0,i=0; i<nsub; i++){
    if (fabs(W[i]) > eps){
      info=EigMatSetEig(AA->Eig,neigs,W[i],invp,dmatarray+nsub*i,nsub,n);DSDPCHKERR(info);
      neigs++;
    }
  }
  
  if (ownarray1){ DSDPFREE(&dmatarray,&info);DSDPCHKERR(info);}
  if (ownarray2){ DSDPFREE(&dworkarray,&info);DSDPCHKERR(info);}
  if (ownarray3){ DSDPFREE(&i2darray,&info);DSDPCHKERR(info);}
  return 0;
}
 
