#include "dsdpsys.h"
#include "sdpconevec.h"
#include "dsdplapack.h"
/*!
\file sdpconevec.c
\brief Implementation of the SDPCone vector operations.
*/

#define SDPConeVecCheck(a,b)    {if (a.dim != b.dim) return 1; if (a.dim>0 && (a.val==NULL || b.val==NULL) ) return 2;}
static int nvecs=0;

#undef __FUNCT__
#define __FUNCT__ "SDPConeVecCreate"
int SDPConeVecCreate(int n ,SDPConeVec *V){
  int info;
  V->dim=n;
  if (n>0){
    nvecs++;
    DSDPCALLOC2(&(V->val),double,n,&info);DSDPCHKERR(info);
  } else {
    V->val=NULL;
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeVecDestroy"
int SDPConeVecDestroy(SDPConeVec *V){
  int info;
  if ((*V).val){ 
    DSDPFREE(&(*V).val,&info);DSDPCHKERR(info);
    nvecs--;
  }

  (*V).dim=0;
  (*V).val=0;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeVecView"
/*!
\fn int SDPConeVecView(SDPConeVec V)

\brief Print the elements of the vector
\param V the vector

*/
int SDPConeVecView(SDPConeVec V){
  int i;
  for (i=0; i<V.dim; i++){
    printf("%3.3e ",V.val[i]);
  }
  printf("\n");
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeVecZero"
/*!
\fn int SDPConeVecZero(SDPConeVec V)

\brief Zero the elements of the vector
\param V the vector

*/
int SDPConeVecZero(SDPConeVec V){
  int n=V.dim;
  double *v=V.val;
  memset((void*)v,0,n*sizeof(double));
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeVecNormalize"
/*!
\fn int SDPConeVecNormalize(SDPConeVec V)

\brief Scale the vector to norm of 1.
\param V the vector

*/
int SDPConeVecNormalize(SDPConeVec V){
  int info;
  double vnorm;
  info = SDPConeVecNorm2(V,&vnorm);DSDPCHKERR(info);
  if (vnorm==0){ return 1;}
  vnorm=1.0/(vnorm);
  info = SDPConeVecScale(vnorm,V);DSDPCHKERR(info);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeVecCopy"
/*!
\fn int SDPConeVecCopy(SDPConeVec v1, SDPConeVec v2)

\brief Copy v1 to v2
\param v1 source
\param v2 destination
*/
int SDPConeVecCopy( SDPConeVec v1,  SDPConeVec v2){

  int n=v1.dim;
  double *val1=v1.val,*val2=v2.val;
  SDPConeVecCheck(v1,v2);
  if (val1!=val2){
    memcpy(val2,val1,n*sizeof(double));
  }
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeVecDot"
/*!
\fn int SDPConeVecDot(SDPConeVec V1, SDPConeVec V2, double *ans)

\brief Inner product of two vectors.
\param V1 a vec
\param V2 a vec
\param ans the inner product
*/
int SDPConeVecDot(SDPConeVec V1, SDPConeVec V2, double *ans){
  ffinteger ione=1, nn=V1.dim;
  double *v1=V1.val,*v2=V2.val;
  *ans=ddot(&nn,v1,&ione,v2,&ione);
  if (*ans!=*ans) return 1;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeVecNorm2"
/*!
\fn int SDPConeVecNorm2(SDPConeVec VV, double *vnorm)

\brief Compute the Euclidean norm.
\param VV a vec
\param vnorm its norm
*/
int SDPConeVecNorm2( SDPConeVec VV, double *vnorm){
  ffinteger ione=1,nn=VV.dim;
  double dd,*v=VV.val;
  dd=dnrm2(&nn,v,&ione);
  *vnorm = dd;
  if (*vnorm!=*vnorm) return 1;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeVecScale"
/*!
\fn int SDPConeVecScale(double alpha, SDPConeVec VV)

\brief Compute the Euclidean norm.
\param alpha scalar.
\param VV a vec
*/
int SDPConeVecScale(double alpha, SDPConeVec VV){
  ffinteger ione=1,nn=VV.dim;
  double *v=VV.val;
  dscal(&nn,&alpha,v,&ione);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeVecAXPY"
/*!
\fn int SDPConeVecAXPY(double alpha,  SDPConeVec x,  SDPConeVec y);

\brief Add a multiple of X to Y.
\param alpha scalar 
\param x a vec
\param y a vec
*/
int SDPConeVecAXPY(double alpha,  SDPConeVec x,  SDPConeVec y){
  ffinteger ione=1,nn=x.dim;
  double *yy=y.val,*xx=x.val;
  if (alpha==0) return 0;
  daxpy(&nn,&alpha,xx,&ione,yy,&ione);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeVecDuplicate"
/*!
\fn int SDPConeVecDuplicate(SDPConeVec V1,SDPConeVec *V2);

\brief Allocate another vector with the same structure as the first
\param V1 source vector
\param V2 new vector
*/
int SDPConeVecDuplicate(SDPConeVec V1,SDPConeVec *V2){
  int info,n=V1.dim;
  info = SDPConeVecCreate(n ,V2);DSDPCHKERR(info);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "SDPConeVecSet"
/*!
\fn int SDPConeVecSet(double alpha, SDPConeVec V)

\brief Set each element of vector to this number.
\param alpha scalar.
\param V a vec
*/
int SDPConeVecSet(double alpha, SDPConeVec V){

  int i,n=V.dim;
  double *val=V.val;

  if (alpha==0.0){
    memset((void*)val,0,n*sizeof(double));
    return 0;
  }
  for (i=0; i<n; ++i){
    val[i]= alpha;
  }
  return 0;
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPIndexInitialize"
/*!
\fn int DSDPIndexInitialize(DSDPIndex *IS);

\brief Set structure pointers to 0.
\param IS indices
*/
int DSDPIndexInitialize(DSDPIndex *IS){
  DSDPFunctionBegin;
  IS->indx=0;
  DSDPFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "DSDPIndexCreate"
/*!
\fn int DSDPIndexCreate(int n, DSDPIndex *IS);

\brief Allocate array for indices
\param n dimension of block or vector associated with it.
\param IS indices
*/
int DSDPIndexCreate(int n,DSDPIndex *IS){
  int info,*is;
  DSDPFunctionBegin;
  DSDPCALLOC2(&is,int,n+1,&info);
  IS->indx=is;is[0]=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPIndexDestroy"
/*!
\fn int DSDPIndexDestroy(DSDPIndex *IS);

\brief Deallocate memory 
\param IS indices
*/
int DSDPIndexDestroy(DSDPIndex *IS){
  int info;
  DSDPFunctionBegin;
  DSDPFREE(&IS->indx,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPIndexView"
/*!
\fn int DSDPIndexView(DSDPIndex IS);

\brief Print indices
\param IS indices
*/
int DSDPIndexView(DSDPIndex IS){
  int i;
  DSDPFunctionBegin;
  printf("Index Set with %d indices.\n",IS.indx[0]);
  for (i=0;i<IS.indx[0];i++){
    printf(" %d",IS.indx[i+1]);
  }
  printf(" \n");
  DSDPFunctionReturn(0);
}
