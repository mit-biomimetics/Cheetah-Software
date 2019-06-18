#if !defined(__DSDP_DSMATRIX_H) 
#define __DSDP_DSMATRIX_H

#include "sdpconevec.h"
#include "dsdpxmat.h"

/*!
\file dsdpdsmat.h
\brief The interface between the SDPCone and the Delta S matrix
*/

/* DSDPDSMat objects are not used for much:  DS, X, eigenvalue stuff     */
/* These objects are good basically for assembling a matrix, accessing
   the data, and applying the operator to a vector                     */

/* DSDP Matrix Structure */
/*!
struct DSDPDSMat_C { void* matdata; struct DSDPDSMat_Ops* dsdpops; };

\brief Symmetric Delta S matrix for one block in the semidefinite cone.
\sa DSDPDSMat
*/
struct DSDPDSMat_C{
  void  *matdata;
  struct DSDPDSMat_Ops* dsdpops;
};

/*!
typedef struct DSDPDSMat_C DSDPDSMat;

\brief A symmetric Delta S matrix for one block in the semidefinite cone.
*/
typedef struct DSDPDSMat_C DSDPDSMat;

#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPDSMatGetType(DSDPDSMat, int *);
extern int DSDPDSMatSetData(DSDPDSMat *, struct DSDPDSMat_Ops*,  void*);
extern int DSDPDSMatInitialize(DSDPDSMat*);

extern int DSDPDSMatZeroEntries(DSDPDSMat);
extern int DSDPDSMatSetArray(DSDPDSMat,DSDPVMat);
extern int DSDPDSMatMult(DSDPDSMat, SDPConeVec, SDPConeVec);
extern int DSDPDSMatVecVec(DSDPDSMat, SDPConeVec, double*);

extern int DSDPDSMatTest(DSDPDSMat);
extern int DSDPDSMatGetSize(DSDPDSMat,int*);
extern int DSDPDSMatView(DSDPDSMat);
extern int DSDPDSMatDestroy(DSDPDSMat*);

extern int DSDPDSMatCheck(DSDPDSMat,SDPConeVec,SDPConeVec,DSDPVMat);
#ifdef __cplusplus
}
#endif


#endif


