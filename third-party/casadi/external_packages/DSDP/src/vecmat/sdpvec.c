#include "dsdpsys.h"
#include "dsdpvec.h"
#include "dsdplapack.h"
/*!
\file sdpvec.c
\brief DSDPVec operations
*/
#if !defined (min)
#define min(a,b) ((a <= b)? (a) : (b))
#endif
#if !defined (max)
#define max(a,b) ((a >= b)? (a) : (b))
#endif

#define DSPPVecCheck(a,b)    {if (a.dim != b.dim) return 1; if (a.dim>0 && (a.val==NULL || b.val==NULL) ) return 2;}

static int nvecs=0;
#undef __FUNCT__
#define __FUNCT__ "DSDPVecCreateSeq"
int DSDPVecCreate(DSDPVec *V){
  int info;
  info = DSDPVecCreateSeq(0,V);DSDPCHKERR(info);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecCreateSeq"
int DSDPVecCreateSeq(int n ,DSDPVec *V){
  int info;
  V->dim=n;
  if (n>0){
    nvecs++;
    DSDPCALLOC2(&(V->val),double,n,&info);DSDPCHKERR(info);
    if (V->val==NULL) return 1;
  } else {
    V->val=NULL;
  }
  return 0;
}
/*
#undef __FUNCT__
#define __FUNCT__ "DSDPVecCreateWArray"
int DSDPVecCreateWArray(DSDPVec *V, double* vv, int n){
  V->dim=n;
  if (n>0){
    V->val=vv;
  } else {
    V->val=NULL;
  }
  return 0;
}
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPVecDestroy"
int DSDPVecDestroy(DSDPVec *V){
  int info;
  if ((*V).val){ 
    DSDPFREE(&(*V).val,&info);DSDPCHKERR(info);
    nvecs--;
  }

  (*V).dim=0;
  (*V).val=0;
  return 0;
}

/*
int DSDPVecGetSize(DSDPVec V, int *n){
  *n=V.dim;
  return 0;
}

int DSDPVecGetArray(DSDPVec V, double **dptr){
  *dptr=V.val;
  return 0;
}

int DSDPVecRestoreArray(DSDPVec V, double **dptr){
  *dptr=0;
  return 0;
}
*/


#undef __FUNCT__
#define __FUNCT__ "DSDPVecView"
int DSDPVecView(DSDPVec vec){
  int i;
  for (i=0; i<vec.dim; i++){
    printf("%3.3e ",vec.val[i]);
  }
  printf("\n");
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecISet"
int DSDPVecISet(int* ival,DSDPVec V){
  int i;
  for (i=0;i<V.dim;i++){
    V.val[i]=ival[i];
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecSetValue"
int DSDPVecSetValue(DSDPVec V,int row,double value){
  V.val[row]=value;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecZero"
int DSDPVecZero(DSDPVec V){
  int n=V.dim;
  double *v=V.val;
  memset((void*)v,0,n*sizeof(double));
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVecNormalize"
int DSDPVecNormalize(DSDPVec V){
  int info;
  double vnorm;
  info = DSDPVecNorm2(V,&vnorm);DSDPCHKERR(info);
  if (vnorm==0){ return 1;}
  vnorm=1.0/(vnorm);
  info = DSDPVecScale(vnorm,V);DSDPCHKERR(info);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecSetBasis"
int DSDPVecSetBasis(DSDPVec V,int row){
  int info;
  info=DSDPVecZero(V);
  V.val[row]=1.0;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecCopy"
int DSDPVecCopy( DSDPVec v1,  DSDPVec v2){

  int n=v1.dim;
  double *val1=v1.val,*val2=v2.val;
  DSPPVecCheck(v1,v2);
  if (val1!=val2){
    memcpy(val2,val1,n*sizeof(double));
  }
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVecSum"
int DSDPVecSum( DSDPVec v, double *vnorm){
  int i,n;
  n = v.dim;
  *vnorm = 0.0;
  for (i=0; i<n; i++){
    *vnorm += v.val[i];
  }
  if (*vnorm!=*vnorm) return 1;
  return 0;
}
#undef __FUNCT__
#define __FUNCT__ "DSDPVecNorm1"
int DSDPVecNorm1( DSDPVec v, double *vnorm){
  ffinteger N=v.dim,INCX=1;
  *vnorm=0;
  *vnorm=dasum(&N,v.val,&INCX);
  if (*vnorm!=*vnorm) return 1;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVecDot"
int DSDPVecDot(DSDPVec V1, DSDPVec V2, double *ans){
  ffinteger ione=1, nn=V1.dim;
  double *v1=V1.val,*v2=V2.val;
  *ans=ddot(&nn,v1,&ione,v2,&ione);
  if (*ans!=*ans) return 1;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVecNorm22"
int DSDPVecNorm22( DSDPVec VV, double *vnorm){
  ffinteger ione=1,nn=VV.dim;
  double dd,*v=VV.val;
  dd=dnrm2(&nn,v,&ione);
  *vnorm = dd*dd;
  if (*vnorm!=*vnorm) return 1;
  return 0;
}
#undef __FUNCT__
#define __FUNCT__ "DSDPVecNorm2"
int DSDPVecNorm2( DSDPVec VV, double *vnorm){
  ffinteger ione=1,nn=VV.dim;
  double dd,*v=VV.val;
  dd=dnrm2(&nn,v,&ione);
  *vnorm = dd;
  if (*vnorm!=*vnorm) return 1;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecScale"
int DSDPVecScale(double alpha, DSDPVec VV){
  ffinteger ione=1,nn=VV.dim;
  double *v=VV.val;
  dscal(&nn,&alpha,v,&ione);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecAXPY"
int DSDPVecAXPY(double alpha,  DSDPVec x,  DSDPVec y){
  ffinteger ione=1,nn=x.dim;
  double *yy=y.val,*xx=x.val;
  if (alpha==0) return 0;
  daxpy(&nn,&alpha,xx,&ione,yy,&ione);
  return 0;
}

/*
#undef __FUNCT__
#define __FUNCT__ "DSDPVecNorm22"
int DSDPVecNorm22( DSDPVec v, double *vnorm){
  int i,n=v.dim;
  double *val=v.val;

  *vnorm = 0.0;
  for (i=0; i<n; i++){
    *vnorm += val[i]*val[i];
  }
  return 0;
}
#undef __FUNCT__
#define __FUNCT__ "DSDPVecNorm2"
int DSDPVecNorm2( DSDPVec v, double *vnorm){
  int info;
  info=DSDPVecNorm22(v,vnorm); if (info) return 1;
  if (*vnorm!=*vnorm) return 1;
  *vnorm = sqrt(*vnorm);
  return 0;
}
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPVecNormInfinity"
int DSDPVecNormInfinity( DSDPVec v, double *vnorm){

  int i,n=v.dim;
  double *val=v.val;

  *vnorm = 0.0;

  for (i=0; i<n; i++){
    *vnorm = max(*vnorm,fabs(val[i]));
  }
  if (*vnorm!=*vnorm) return 1;

  return 0;
}
/*
#undef __FUNCT__
#define __FUNCT__ "DSDPVecScale"
int DSDPVecScale(double alpha, DSDPVec x){
  int i,ii,n;
  double *xx=x.val;
  n=x.dim;

  for (ii=0; ii<n/4; ++ii){
    i=ii*4;
    xx[i]*= alpha;
    xx[i+1]*= alpha;
    xx[i+2]*= alpha;
    xx[i+3]*= alpha;
  }
  for (i=4*(n/4); i<n; ++i){
    xx[i]*= alpha;
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecAXPY"
int DSDPVecAXPY(double alpha,  DSDPVec x,  DSDPVec y){

  int i,ii,n=x.dim;
  double *yy=y.val,*xx=x.val;

  DSPPVecCheck(x,y);

  for (ii=0; ii<n/4; ++ii){
    i=ii*4;
    yy[i] += (alpha)*xx[i];
    yy[i+1] += (alpha)*xx[i+1];
    yy[i+2] += (alpha)*xx[i+2];
    yy[i+3] += (alpha)*xx[i+3];
  }
  for (i=4*(n/4); i<n; ++i){
    yy[i] += (alpha)*xx[i];
  }

  return 0;
}
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPVecWAXPBY"
int DSDPVecWAXPBY(DSDPVec w, double alpha,  DSDPVec x,  double beta, DSDPVec y){

  int i,ii,n=x.dim;
  double *yy=y.val,*xx=x.val,*ww=w.val;
  DSPPVecCheck(x,y);
  DSPPVecCheck(x,w);

  for (ii=0; ii<n/4; ++ii){
    i=ii*4;
    ww[i] = (alpha)*xx[i] + (beta)*yy[i];
    ww[i+1] = (alpha)*xx[i+1] + (beta)*yy[i+1];
    ww[i+2] = (alpha)*xx[i+2] + (beta)*yy[i+2];
    ww[i+3] = (alpha)*xx[i+3] + (beta)*yy[i+3];
  }
  for (i=4*(n/4); i<n; ++i){
    ww[i] = (alpha)*xx[i] + (beta)*yy[i];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecWAXPY"
int DSDPVecWAXPY(DSDPVec w,double alpha, DSDPVec x, DSDPVec y){  
  int info;
  info=DSDPVecCopy(y,w);
  info=DSDPVecAXPY(alpha,x,w);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVecAYPX"
int DSDPVecAYPX(double alpha,  DSDPVec x,  DSDPVec y){

  int i,ii,n=x.dim;
  double *yy=y.val,*xx=x.val;

  DSPPVecCheck(x,y);
  for (ii=0; ii<n/4; ++ii){
    i=ii*4;
    yy[i]   = xx[i]+(alpha)*yy[i];
    yy[i+1] = xx[i+1]+(alpha)*yy[i+1];
    yy[i+2] = xx[i+2]+(alpha)*yy[i+2];
    yy[i+3] = xx[i+3]+(alpha)*yy[i+3];
  }
  for (i=4*(n/4); i<n; ++i){
    yy[i]   = xx[i]+(alpha)*yy[i];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecAYPX"
int DSDPVecScaleCopy(DSDPVec x,  double alpha, DSDPVec y){

  int i,ii,n=x.dim;
  double *yy=y.val,*xx=x.val;

  DSPPVecCheck(x,y);
  for (ii=0; ii<n/4; ++ii){
    i=ii*4;
    yy[i]   = (alpha)*xx[i];
    yy[i+1] = (alpha)*xx[i+1];
    yy[i+2] = (alpha)*xx[i+2];
    yy[i+3] = (alpha)*xx[i+3];
  }
  for (i=4*(n/4); i<n; ++i){
    yy[i]   = (alpha)*xx[i];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecDuplicate"
int DSDPVecDuplicate(DSDPVec V1,DSDPVec *V2){
  int info,n=V1.dim;
  info = DSDPVecCreateSeq(n ,V2);DSDPCHKERR(info);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DSDPVecSet"
int DSDPVecSet(double alpha, DSDPVec V){

  int i,ii,n=V.dim;
  double *val=V.val;

  if (alpha==0.0){
    memset((void*)val,0,n*sizeof(double));
    return 0;
  }
  for (ii=0; ii<n/4; ++ii){
    i=ii*4;
    val[i] = val[i+1] = val[i+2] = val[i+3] = alpha; 
  }
  for (i=4*(n/4); i<n; ++i){
    val[i]= alpha;
  }
  return 0;
}

/*
#undef __FUNCT__
#define __FUNCT__ "DSDPVecDot"
int DSDPVecDot(DSDPVec V1, DSDPVec V2, double *ans){

  int i,ii,m=V1.dim;
  double *v1=V1.val,*v2=V2.val;

  DSPPVecCheck(V1,V2);
  *ans=0.0;
  for (ii=0; ii<m/4; ++ii){
    i=ii*4;
    *ans += v1[i]*v2[i] + v1[i+1]*v2[i+1] + v1[i+2]*v2[i+2] + v1[i+3]*v2[i+3] ; 
  }
  for (i=4*(m/4); i<m; ++i){
    *ans += v1[i]*v2[i];
  }
  if (*ans!=*ans) return 1;

  return 0;
}
*/

#undef __FUNCT__
#define __FUNCT__ "DSDPVecPointwiseMin"
int DSDPVecPointwiseMin( DSDPVec V1, DSDPVec V2, DSDPVec V3){
  
  int i,n=V1.dim;
  double *v1=V1.val,*v2=V2.val,*v3=V3.val;

  DSPPVecCheck(V1,V3);
  DSPPVecCheck(V2,V3);
  for (i=0; i<n; ++i){
    v3[i]=DSDPMin(v2[i],v1[i]);
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecPointwiseMax"
int DSDPVecPointwiseMax( DSDPVec V1, DSDPVec V2, DSDPVec V3){
  
  int i,n=V1.dim;
  double *v1=V1.val,*v2=V2.val,*v3=V3.val;

  DSPPVecCheck(V1,V3);
  DSPPVecCheck(V2,V3);
  for (i=0; i<n; ++i){
    v3[i]=DSDPMax(v2[i],v1[i]);
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecPointwiseMult"
int DSDPVecPointwiseMult( DSDPVec V1, DSDPVec V2, DSDPVec V3){
  
  int ii,i,n=V1.dim;
  double *v1=V1.val,*v2=V2.val,*v3=V3.val;

  DSPPVecCheck(V1,V3);
  DSPPVecCheck(V2,V3);
  for (ii=0; ii<n/4; ++ii){
    i=ii*4;
    v3[i]=v1[i]*v2[i]; 
    v3[i+1]=v1[i+1]*v2[i+1]; 
    v3[i+2]=v1[i+2]*v2[i+2]; 
    v3[i+3]=v1[i+3]*v2[i+3];
  }
  for (i=4*(n/4); i<n; i++){
    v3[i]=v1[i]*v2[i];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecPointwiseDivide"
int DSDPVecPointwiseDivide( DSDPVec V1, DSDPVec V2, DSDPVec V3){
  
  int ii,i,n=V1.dim;
  double *v1=V1.val,*v2=V2.val,*v3=V3.val;

  DSPPVecCheck(V1,V3);
  DSPPVecCheck(V2,V3);
  for (ii=0; ii<n/4; ++ii){
    i=ii*4;
    v3[i]=v1[i]/v2[i]; v3[i+1]=v1[i+1]/v2[i+1]; v3[i+2]=v1[i+2]/v2[i+2]; v3[i+3]=v1[i+3]/v2[i+3];
  }
  for (i=4*(n/4); i<n; i++){
    v3[i]=v1[i]/v2[i];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecShift"
int DSDPVecShift(double alpha, DSDPVec V){
  int i,n=V.dim;
  double *v=V.val;
  for (i=0; i<n; i++){
    v[i]+= alpha;
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecSemiNorm"
int DSDPVecSemiNorm(DSDPVec V, double *ans){

  int i;
  double dtmp=0.0;

  for (i=0; i<V.dim; i++){
    dtmp=min(V.val[i],dtmp);
  }
  *ans = fabs(dtmp);
  if (*ans!=*ans) return 1;
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecReciprocal"
int DSDPVecReciprocal(DSDPVec V){

  int i,n=V.dim;
  double *val=V.val;

  for (i=0; i<n; i++){
    val[i]= 1.0/val[i];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecReciprocalSqrt"
int DSDPVecReciprocalSqrt(DSDPVec V){

  int i,n=V.dim;
  double *val=V.val;

  for (i=0; i<n; i++){
    val[i]= sqrt(1.0/val[i]);
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPVecAbsoluteValue"
int DSDPVecAbsoluteValue(DSDPVec V){

  int i,n=V.dim;
  double *val=V.val;
  for (i=0; i<n; i++){
    val[i]=fabs(val[i]);
  }
  return 0;
}


