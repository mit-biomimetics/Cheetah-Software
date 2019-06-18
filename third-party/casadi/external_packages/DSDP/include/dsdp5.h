#if !defined(__DSDP_INTERFACE_H) 
#define __DSDP_INTERFACE_H
/*!
  \file dsdp5.h
  \brief The API to DSDP for those applications using DSDP as a subroutine library
*/

#include "dsdpbasictypes.h"
#include "dsdpsys.h"

/*! 
\typedef struct SDPCone_C* SDPCone;
\brief The SDPCone object points to blocks of data that specify
semidefinite matrix inequalities.
*/
/*! 
\typedef struct LPCone_C* LPCone;
\brief The LPCone object points to blocks of data that specify
linear scalar inequality constraints.
*/
/*! 
\typedef struct BCone_C* BCone;
\brief The BCone object points to lower and upper bounds on 
the variable y in (D).
*/
typedef struct SDPCone_C* SDPCone;
typedef struct LPCone_C* LPCone;
typedef struct BCone_C* BCone;
extern FILE *dsdpoutputfile;

#ifdef __cplusplus
extern "C" {
#endif


extern int DSDPCreate(int, DSDP *);
extern int DSDPSetup(DSDP);
extern int DSDPSolve(DSDP);
extern int DSDPComputeX(DSDP);
extern int DSDPComputeAndFactorS(DSDP,DSDPTruth*);
extern int DSDPDestroy(DSDP);

extern int DSDPCreateBCone(DSDP, BCone*);
extern int BConeAllocateBounds(BCone,int);
extern int BConeSetLowerBound(BCone,int, double);
extern int BConeSetUpperBound(BCone,int, double);
extern int BConeSetUnboundedLower(BCone,int);
extern int BConeSetUnboundedUpper(BCone,int);
extern int BConeSetPSlackVariable(BCone,int);
extern int BConeSetPSurplusVariable(BCone,int);
extern int BConeScaleBarrier(BCone,double);
extern int BConeView(BCone);
extern int BConeSetXArray(BCone,double[], int);
extern int BConeCopyXSingle(BCone,double [], int m);
extern int BConeCopyX(BCone,double[],double[],int);

extern int DSDPBoundDualVariables(DSDP,double, double);
extern int DSDPSetYBounds(DSDP,double,double);
extern int DSDPGetYBounds(DSDP,double*,double*);

extern int DSDPCreateLPCone(DSDP,LPCone*);
extern int LPConeSetData(LPCone,int,const int[],const int[],const double[]);
extern int LPConeSetData2(LPCone,int,const int[],const int[],const double[]);
extern int LPConeSetDataC(LPCone lpcone,int n, const double vals[]);
extern int LPConeGetData(LPCone,int,double[],int);
extern int LPConeScaleBarrier(LPCone,double);
extern int LPConeGetXArray(LPCone,double*[], int*);
extern int LPConeGetSArray(LPCone,double*[], int*);
extern int LPConeGetDimension(LPCone,int*);
extern int LPConeView(LPCone lpcone);
extern int LPConeView2(LPCone lpcone);
extern int LPConeCopyS(LPCone,double[],int);

extern int DSDPCreateSDPCone(DSDP, int, SDPCone*);
extern int SDPConeSetBlockSize(SDPCone, int, int);
extern int SDPConeGetBlockSize(SDPCone, int, int*);
extern int SDPConeSetStorageFormat(SDPCone, int, char);
extern int SDPConeGetStorageFormat(SDPCone, int, char*);
extern int SDPConeCheckStorageFormat(SDPCone, int, char);
extern int SDPConeSetSparsity(SDPCone, int, int);
extern int SDPConeView(SDPCone);
extern int SDPConeView2(SDPCone);
extern int SDPConeView3(SDPCone);

extern int SDPConeSetASparseVecMat(SDPCone,int,int,int,double,int,const int[],const double[],int);
extern int SDPConeSetADenseVecMat(SDPCone,int,int,int,double,double[],int);
/* extern int SDPConeSetDenseMatWSparseData(SDPCone,int, int, int, double*, int*, int*); */
extern int SDPConeSetARankOneMat(SDPCone,int,int,int,double,int,const int[],const double[],int);
extern int SDPConeSetConstantMat(SDPCone,int,int,int,double); 
extern int SDPConeSetZeroMat(SDPCone,int,int,int);
extern int SDPConeSetIdentity(SDPCone,int,int,int,double);
extern int SDPConeViewDataMatrix(SDPCone,int,int);
extern int SDPConeMatrixView(SDPCone,int);

extern int SDPConeAddASparseVecMat(SDPCone,int,int,int,double,int,const int[],const double[],int);
extern int SDPConeAddADenseVecMat(SDPCone,int,int,int,double,double[],int);
extern int SDPConeAddConstantMat(SDPCone,int,int,int,double); 
extern int SDPConeAddIdentity(SDPCone,int,int,int,double);
extern int SDPConeAddARankOneMat(SDPCone,int,int,int,double,int,const int[],const double[],int);

/* For backward compatibility */
extern int SDPConeAddSparseVecMat(SDPCone,int,int,int,int,const int[],const double[],int);
extern int SDPConeAddDenseVecMat(SDPCone,int,int,int,double[],int);
extern int SDPConeSetSparseVecMat(SDPCone,int,int,int,int,const int[],const double[],int);
extern int SDPConeSetDenseVecMat(SDPCone,int,int,int,double[],int);

extern int SDPConeSetXMat(SDPCone,int,int);
extern int SDPConeSetXArray(SDPCone,int,int,double[], int);
extern int SDPConeGetXArray(SDPCone,int,double*[],int*);
extern int SDPConeRestoreXArray(SDPCone,int,double*[],int*);
extern int SDPConeCheckData(SDPCone);
extern int SDPConeRemoveDataMatrix(SDPCone,int,int);
extern int SDPConeGetNumberOfBlocks(SDPCone, int*);
extern int SDPConeComputeS(SDPCone, int, double,double[], int, double,int,double[],int);
extern int SDPConeComputeX(SDPCone,int,int,double[],int);
extern int SDPConeAddADotX(SDPCone,int,double,double[],int,double[],int);
extern int SDPConeViewX(SDPCone,int,int,double[],int);
extern int SDPConeSetLanczosIterations(SDPCone,int);
extern int SDPConeScaleBarrier(SDPCone,int,double); 
extern int SDPConeXVMultiply(SDPCone,int,double[],double[],int);
extern int SDPConeComputeXV(SDPCone,int,int*);
extern int SDPConeAddXVAV(SDPCone,int,double[],int,double[],int);
extern int SDPConeUseLAPACKForDualMatrix(SDPCone,int);

extern int DSDPSetDualObjective(DSDP,int,double);
extern int DSDPAddObjectiveConstant(DSDP,double);
extern int DSDPGetDObjective(DSDP,double*);
extern int DSDPGetDDObjective(DSDP,double*);
extern int DSDPGetPObjective(DSDP,double*);
extern int DSDPGetPPObjective(DSDP,double*);
/*
extern int DSDPGetDualObjective(DSDP,double*); 
extern int DSDPGetPrimalObjective(DSDP,double*); 
*/
#define DSDPGetDualObjective   DSDPGetDDObj
#define DSDPGetPrimalObjective DSDPGetPPObj
extern int DSDPGetDualityGap(DSDP,double*); 
extern int DSDPGetScale(DSDP,double*);
extern int DSDPSetScale(DSDP,double);
extern int DSDPGetPenaltyParameter(DSDP,double*);
extern int DSDPGetPenalty(DSDP,double*);
extern int DSDPCopyB(DSDP,double[], int);

extern int DSDPSetR0(DSDP,double);
extern int DSDPGetR(DSDP, double *);
extern int DSDPSetRTolerance(DSDP,double);
extern int DSDPGetRTolerance(DSDP,double*);

extern int DSDPSetY0(DSDP,int, double);
extern int DSDPGetY(DSDP, double[], int);
extern int DSDPGetYMakeX(DSDP, double[], int);
extern int DSDPGetDYMakeX(DSDP, double[], int);
extern int DSDPGetMuMakeX(DSDP, double*);

#define DSDPSetInitialBarrierParameter DSDPSetBarrierParameter
#define DSDPGetInitialBarrierParameter DSDPGetBarrierParameter
extern int DSDPGetBarrierParameter(DSDP, double *);
extern int DSDPSetBarrierParameter(DSDP, double);
extern int DSDPReuseMatrix(DSDP,int);
extern int DSDPGetReuseMatrix(DSDP,int*);
extern int DSDPGetDimension(DSDP, double*);

extern int DSDPSetMaxIts(DSDP,int); 
extern int DSDPGetMaxIts(DSDP,int*); 
extern int DSDPSetStepTolerance(DSDP,double);
extern int DSDPGetStepTolerance(DSDP,double*);
extern int DSDPSetGapTolerance(DSDP,double);
extern int DSDPGetGapTolerance(DSDP,double*);
extern int DSDPSetPNormTolerance(DSDP,double);
extern int DSDPGetPNormTolerance(DSDP,double*);
extern int DSDPSetDualBound(DSDP,double);
extern int DSDPGetDualBound(DSDP,double*);
extern int DSDPSetPTolerance(DSDP,double);
extern int DSDPGetPTolerance(DSDP,double*);
extern int DSDPGetPInfeasibility(DSDP,double*);
extern int DSDPSetMaxTrustRadius(DSDP,double);
extern int DSDPGetMaxTrustRadius(DSDP,double*);
extern int DSDPStopReason(DSDP,DSDPTerminationReason *); 
extern int DSDPGetSolutionType(DSDP,DSDPSolutionType*);
extern int DSDPSetPotentialParameter(DSDP, double);
extern int DSDPGetPotentialParameter(DSDP, double*);
extern int DSDPUseDynamicRho(DSDP, int);
extern int DSDPGetPotential(DSDP,double*);
extern int DSDPUseLAPACKForSchur(DSDP,int);
extern int DSDPGetNumberOfVariables(DSDP,int*);
extern int DSDPGetFinalErrors(DSDP,double[6]);
extern int DSDPGetGapHistory(DSDP, double[], int); 
extern int DSDPGetRHistory(DSDP, double[], int); 
extern int DSDPGetIts(DSDP,int *); 
extern int DSDPGetPnorm(DSDP, double *);
extern int DSDPGetStepLengths(DSDP, double*,double*);
extern int DSDPSetMonitor(DSDP, int (*)(DSDP,void*),void*);
extern int DSDPSetStandardMonitor(DSDP,int);
extern int DSDPSetFileMonitor(DSDP,int);
extern int DSDPSetPenaltyParameter(DSDP,double);
extern int DSDPUsePenalty(DSDP,int);
extern int DSDPPrintLogInfo(int);
extern int DSDPComputeMinimumXEigenvalue(DSDP, double*);
extern int DSDPGetTraceX(DSDP dsdp, double*);
extern int DSDPSetZBar(DSDP,double); 
extern int DSDPSetDualLowerBound(DSDP, double);
extern int DSDPGetDataNorms(DSDP, double[3]);
extern int DSDPGetYMaxNorm(DSDP, double*);
extern int SDPConeUseFullSymmetricFormat(SDPCone, int);
extern int SDPConeUsePackedFormat(SDPCone, int);
extern int DSDPSetFixedVariable(DSDP,int,double);
extern int DSDPSetFixedVariables(DSDP,double[],double[],double[],int);
extern int DSDPGetFixedYX(DSDP,int,double*);
extern int DSDPView(DSDP);
extern int DSDPPrintOptions();
extern int DSDPPrintData(DSDP,SDPCone,LPCone);
extern int DSDPPrintSolution(FILE*,DSDP,SDPCone, LPCone);
extern int DSDPSetOptions(DSDP,char*[], int);
extern int DSDPReadOptions(DSDP, char[]);
extern int DSDPSetDestroyRoutine(DSDP, int (*)(void*), void*);

#ifdef __cplusplus
}
#endif

#endif
