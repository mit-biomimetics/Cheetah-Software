#if !defined(__DSDP_VSYMMETRICMATRIX_H) 
#define __DSDP_VSYMMETRICMATRIX_H

/*! 
\file dsdpxmat.h
\brief The interface between the SDPCone and the dense matrix array.
 */
#include "sdpconevec.h"

/* DSDP V Matrix Structure */
/*!
struct DSDPVMat_C { void* matdata; struct DSDPVMat_Ops* dsdpops; };

\brief Dense symmetric matrix for one block in the semidefinite cone.
\sa DSDPVMat
*/
struct DSDPVMat_C{
  void  *matdata;
  struct DSDPVMat_Ops* dsdpops;
};

/*!
\typedef struct DSDPVMat_C   DSDPVMat;
\brief Represents a dense symmetric matrix for one block in the semidefinite cone.
*/
typedef struct DSDPVMat_C DSDPVMat;

#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPVMatGetType(DSDPVMat, int *);
extern int DSDPVMatSetData(DSDPVMat *, struct DSDPVMat_Ops*,  void*);
extern int DSDPVMatInitialize(DSDPVMat*);

extern int DSDPVMatTest(DSDPVMat);
extern int DSDPVMatGetSize(DSDPVMat,int*);
extern int DSDPVMatView(DSDPVMat);
extern int DSDPVMatDestroy(DSDPVMat*);

extern int DSDPVMatExist(DSDPVMat,int*);
extern int DSDPVMatZeroEntries(DSDPVMat);
extern int DSDPVMatAddOuterProduct(DSDPVMat, double, SDPConeVec);
extern int DSDPVMatMult(DSDPVMat,SDPConeVec,SDPConeVec);
extern int DSDPVMatScaleDiagonal(DSDPVMat,double);
extern int DSDPVMatShiftDiagonal(DSDPVMat,double);
extern int DSDPVMatNormF2(DSDPVMat, double*);
extern int DSDPVMatGetArray(DSDPVMat,double**,int*);
extern int DSDPVMatRestoreArray(DSDPVMat,double**,int*);
extern int DSDPVMatMinEigenvalue(DSDPVMat,SDPConeVec,SDPConeVec,double*);
extern int DSDPVMatCheck(DSDPVMat, SDPConeVec, SDPConeVec);

#ifdef __cplusplus
}
#endif

#endif


