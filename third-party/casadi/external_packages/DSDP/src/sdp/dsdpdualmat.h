#if !defined(__DSDP_DUALMATRIX_H) 
#define __DSDP_DUALMATRIX_H

#include "sdpconevec.h"
#include "dsdpbasictypes.h"
#include "dsdpxmat.h"
/*!
\file dsdpdualmat.h
\brief The interface between the SDPCone and the matrix S.
 */

/*!
struct DSDPDualMat_C { void* matdata; struct DSDPDualMat_Ops* dsdpops; };

\brief Represents an S matrix for one block in the semidefinite cone.
\sa DSDPDualMat
*/
struct DSDPDualMat_C{
  void* matdata;
  struct DSDPDualMat_Ops* dsdpops;
};

/*!
\typedef struct DSDPDualMat_C   DSDPDualMat;
\brief Represents an S matrix for one block in the semidefinite cone.
*/
typedef struct DSDPDualMat_C   DSDPDualMat;

#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPDualMatInitialize(DSDPDualMat*);
extern int DSDPDualMatSetData(DSDPDualMat*,struct DSDPDualMat_Ops*,void*);
extern int DSDPDualMatGetType(DSDPDualMat, int *);

extern int DSDPDualMatGetSize(DSDPDualMat, int*);
extern int DSDPDualMatTest(DSDPDualMat);
extern int DSDPDualMatDestroy(DSDPDualMat *);
extern int DSDPDualMatView(DSDPDualMat);

extern int DSDPDualMatCholeskyFactor(DSDPDualMat,DSDPTruth *);
extern int DSDPDualMatInvert(DSDPDualMat);
extern int DSDPDualMatInverseAdd(DSDPDualMat,double,DSDPVMat);
extern int DSDPDualMatInverseMultiply(DSDPDualMat, DSDPIndex, SDPConeVec, SDPConeVec);
extern int DSDPDualMatCholeskySolveForward(DSDPDualMat, SDPConeVec, SDPConeVec);
extern int DSDPDualMatCholeskySolveBackward(DSDPDualMat, SDPConeVec, SDPConeVec);
extern int DSDPDualMatCholeskyForwardMultiply(DSDPDualMat, SDPConeVec, SDPConeVec);
extern int DSDPDualMatCholeskyBackwardMultiply(DSDPDualMat, SDPConeVec, SDPConeVec);
extern int DSDPDualMatLogDeterminant(DSDPDualMat, double*);
extern int DSDPDualMatIsFull(DSDPDualMat,DSDPTruth*);
extern int DSDPDualMatSetArray(DSDPDualMat,DSDPVMat);
extern int DSDPDualMatCheck(DSDPDualMat,SDPConeVec,SDPConeVec,DSDPIndex,DSDPVMat);
extern int DSDPDualMatGetArray(DSDPDualMat,double*[],int*);

#ifdef __cplusplus
}
#endif

#endif


