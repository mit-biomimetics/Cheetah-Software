#if !defined(__DSDP_VECTORS_H) 
#define __DSDP_VECTORS_H

#include <math.h>
/*! 
\file dsdpvec.h
\brief Vector operations used by the solver
*/
/* Define DSDP Vector Structure */

/*!
\typedef struct DSDPVec_C DSDPVec;

\brief This object hold m+2 variables: a scaling of C, the y
variables, and r.  

Unscaled, this vector is \f$ [ -1 \ y_1 \ \ldots \ y_m \ r ] \f$.

*/
struct  DSDPVec_C{
  int    dim;
  double *val;
};

typedef struct DSDPVec_C DSDPVec;

#define DSDPVecGetArray(a,b)              0;{ *(b)=((a).val); }
#define DSDPVecRestoreArray(a,b)          0;{ *(b)=0; }
#define DSDPVecGetSize(a,b)               0;{ *(b)=((a).dim); }
#define DSDPVecAddElement(a,b,c)          0;{ if (c){((a).val[b])+=(c);}  }
#define DSDPVecSetElement(a,b,c)          0;{  {((a).val[b])=(c); } }
#define DSDPVecGetElement(a,b,c)          0;{ *(c)=((a).val[b]); }
#define DSDPVecSetR(a,b)                  0;{  {((a).val[(a).dim-1])=(b); } }
#define DSDPVecAddR(a,b)                  0;{  if(b){((a).val[(a).dim-1])+=(b); } }
#define DSDPVecGetR(a,b)                  0;{ *(b)=((a).val[(a).dim-1]); }
#define DSDPVecSetC(a,b)                  0;{  {((a).val[0])=(b); } }
#define DSDPVecAddC(a,b)                  0;{  if(b){((a).val[0])+=(b); } }
#define DSDPVecGetC(a,b)                  0;{ *(b)=((a).val[0]); }
#define DSDPVecCreateWArray(a,b,c)        0;{ (*(a)).val=(b); (*(a)).dim=(c);}
/*
extern int DSDPVecGetArray(DSDPVec, double **);
extern int DSDPVecRestoreArray(DSDPVec, double **);
extern int DSDPVecGetSize(DSDPVec, int *);
extern int DSDPVecAddElement(DSDPVec, int, double);
extern int DSDPVecSetElement(DSDPVec, int, double);
extern int DSDPVecGetElement(DSDPVec, int, double*);
extern int DSDPVecCreateWArray(DSDPVec*, double*, int);
*/
#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPVecCreateSeq(int,DSDPVec *);
extern int DSDPVecDuplicate(DSDPVec,DSDPVec *);
extern int DSDPVecSet(double, DSDPVec );
extern int DSDPVecISet(int*,DSDPVec);
extern int DSDPVecZero(DSDPVec );
extern int DSDPVecNormalize(DSDPVec );
extern int DSDPVecSetValue(DSDPVec,int,double);
extern int DSDPVecSetBasis(DSDPVec,int);
extern int DSDPVecCopy( DSDPVec,  DSDPVec);
extern int DSDPVecScale(double, DSDPVec);
extern int DSDPVecScaleCopy(DSDPVec,  double, DSDPVec);
extern int DSDPVecAXPY(double,  DSDPVec,  DSDPVec);
extern int DSDPVecAYPX(double,  DSDPVec,  DSDPVec);
extern int DSDPVecWAXPY(DSDPVec,double,DSDPVec,DSDPVec);
extern int DSDPVecWAXPBY(DSDPVec,double,DSDPVec,double,DSDPVec);
extern int DSDPVecPointwiseMin( DSDPVec, DSDPVec, DSDPVec);
extern int DSDPVecPointwiseMax( DSDPVec, DSDPVec, DSDPVec);
extern int DSDPVecPointwiseMult( DSDPVec, DSDPVec, DSDPVec);
extern int DSDPVecPointwiseDivide( DSDPVec, DSDPVec, DSDPVec);
extern int DSDPVecReciprocalSqrt(DSDPVec);
extern int DSDPVecDot(DSDPVec, DSDPVec, double *);
extern int DSDPVecSum( DSDPVec, double *);
extern int DSDPVecNorm1( DSDPVec, double *);
extern int DSDPVecNorm2( DSDPVec, double *);
extern int DSDPVecNorm22( DSDPVec, double *);
extern int DSDPVecNormInfinity( DSDPVec, double *);
extern int DSDPVecAbsoluteValue( DSDPVec);
extern int DSDPVecShift(double, DSDPVec);
extern int DSDPVecView( DSDPVec);
extern int DSDPVecDestroy(DSDPVec*);

#ifdef __cplusplus
}
#endif

#endif

