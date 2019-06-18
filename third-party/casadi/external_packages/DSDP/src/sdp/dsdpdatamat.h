#if !defined(__DSDP_DATAMATRIXOPS_H) 
#define __DSDP_DATAMATRIXOPS_H

#include "sdpconevec.h"
/*! \file dsdpdatamat.h
\brief The interface between the SDPCone and the data matrices
 */

/*!
struct DSDPDataMat_C { void* matdata; struct DSDPDataMat_Ops* dsdpops; };

\brief Symmetric data matrix for one block in the semidefinite cone.
\sa DSDPDataMat
*/
struct DSDPDataMat_C{
  void* matdata;
  struct DSDPDataMat_Ops* dsdpops;
};

/*!
\typedef struct DSDPDataMat   DSDPDataMat;
\brief Represents a single symmetric
data matrix for one block in this semidefinite cone.
*/
typedef struct DSDPDataMat_C   DSDPDataMat;

#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPDataMatSetData(DSDPDataMat*, struct DSDPDataMat_Ops*,  void*);
extern int DSDPDataMatInitialize(DSDPDataMat*);
extern int DSDPDataMatGetType(DSDPDataMat, int *);
extern int DSDPDataMatTest(DSDPDataMat);

extern int DSDPDataMatVecVec(DSDPDataMat,SDPConeVec,double*);
extern int DSDPDataMatDot(DSDPDataMat,double[], int,int,double*);
extern int DSDPDataMatGetRowNonzeros(DSDPDataMat, int, int, int*, int*);
extern int DSDPDataMatCountNonzeros(DSDPDataMat,int*,int);
extern int DSDPDataMatFNorm2(DSDPDataMat,int,double*);
extern int DSDPDataMatMultiply(DSDPDataMat,SDPConeVec,SDPConeVec);
extern int DSDPDataMatView(DSDPDataMat);
extern int DSDPDataMatDestroy(DSDPDataMat*);
extern int DSDPDataMatGetRank(DSDPDataMat, int*,int);
extern int DSDPDataMatGetEig(DSDPDataMat, int, SDPConeVec, DSDPIndex, double *);
extern int DSDPDataMatFactor(DSDPDataMat,SDPConeVec, double[],int,double[],int,int[],int);
extern int DSDPDataMatAddMultiple(DSDPDataMat, double, double[],int,int);
extern int DSDPDataMatAddRowMultipleToVector(DSDPDataMat, int, double, SDPConeVec);

#ifdef __cplusplus
}
#endif

#endif


