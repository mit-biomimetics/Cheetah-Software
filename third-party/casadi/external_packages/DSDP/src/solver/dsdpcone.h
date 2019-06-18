#ifndef __DSDPCONEOPERATIONS_H
#define __DSDPCONEOPERATIONS_H

/*!
\file dsdpcone.h
\brief The public interface between the cones and the solver.
*/

#include "dsdpbasictypes.h"
#include "dsdpvec.h"
#include "dsdpschurmat.h"


/*!
\typedef struct DSDPCone DSDPCone;

\brief This object holds the data of a SDP, LP, or other cone.  Its
structure is opaque to the DSDP Solver, but it must implement 
the interface below and provide a structure of function pointers.

*/
struct DSDPCone_C{
  void* conedata;
  struct DSDPCone_Ops* dsdpops;
};

typedef struct DSDPCone_C   DSDPCone;

#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPConeSetUp(DSDPCone,DSDPVec);
extern int DSDPConeSetUp2(DSDPCone, DSDPVec, DSDPSchurMat);
extern int DSDPConeGetDimension(DSDPCone,double*);
extern int DSDPConeSparsityInSchurMat(DSDPCone, int,int[],int);
extern int DSDPConeComputeS(DSDPCone,DSDPVec,DSDPDualFactorMatrix,DSDPTruth *);
extern int DSDPConeInvertS(DSDPCone);
extern int DSDPConeComputeHessian(DSDPCone,double,DSDPSchurMat,DSDPVec,DSDPVec);
extern int DSDPConeMultiplyAdd(DSDPCone,double,DSDPVec,DSDPVec,DSDPVec);
extern int DSDPConeComputeRHS(DSDPCone,double,DSDPVec,DSDPVec,DSDPVec);
extern int DSDPConeComputeMaxStepLength(DSDPCone,DSDPVec,DSDPDualFactorMatrix,double*);
extern int DSDPConeComputeLogSDeterminant(DSDPCone,double*,double*);
extern int DSDPConeComputeX(DSDPCone,double,DSDPVec,DSDPVec,DSDPVec,double*);
extern int DSDPConeSetXMaker(DSDPCone,double,DSDPVec,DSDPVec);
extern int DSDPConeView(DSDPCone);
extern int DSDPConeMonitor(DSDPCone,int);
extern int DSDPConeDestroy(DSDPCone*);
extern int DSDPConeANorm2(DSDPCone,DSDPVec);
extern int DSDPGetConeName(DSDPCone,char*,int);

extern int DSDPConeSetData(DSDPCone*,struct DSDPCone_Ops*,void*);
extern int DSDPConeInitialize(DSDPCone*);

#ifdef __cplusplus
}
#endif

#endif





