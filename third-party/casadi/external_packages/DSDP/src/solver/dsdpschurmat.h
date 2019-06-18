#if !defined(__DSDP_SCHURMATRIXOPERATIONS_H) 
#define __DSDP_SCHURMATRIXOPERATIONS_H

/*!
\file dsdpschurmat.h
\brief Methods of a Schur Matrix.
*/

#include "dsdpvec.h"
#include "dsdpbasictypes.h"
#include "dsdpcg.h"


typedef struct {
  int *var;
  int nvars;
  int nmaxvars;
  double *fval;
  double *fdual;
  double *xout;
} FixedVariables;

typedef struct {
  FixedVariables fv;
  DSDPVec rhs3,dy3;
  double dd,r;
  int m;
} DSDPSchurInfo;

/*!
struct DSDPSchurMat_C
\brief Schur complement matrix whose solution is the Newton direction.
\sa DSDPSchurMat
*/
struct DSDPSchurMat_C{
  void* data;
  struct DSDPSchurMat_Ops *dsdpops;
  DSDPSchurInfo* schur;
};

/*!
\typedef struct DSDPSchurMat_C DSDPSchurMat;

\brief This object represents the Schur Matrix. Its structure is
opaque to the DSDP solver, but it must implement the interface below
and provide a structure of function pointers.

*/
typedef struct DSDPSchurMat_C   DSDPSchurMat;


#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPSchurMatSetData(DSDPSchurMat*,struct DSDPSchurMat_Ops*, void*);

extern int DSDPSchurMatInitialize(DSDPSchurMat*);
extern int DSDPSchurMatSetup(DSDPSchurMat, DSDPVec);
extern int DSDPSchurMatZeroEntries(DSDPSchurMat);
extern int DSDPSchurMatInParallel(DSDPSchurMat, DSDPTruth*);
extern int DSDPSchurMatShiftDiagonal(DSDPSchurMat, double);
extern int DSDPSchurMatAssemble(DSDPSchurMat);
extern int DSDPSchurMatMultiply(DSDPSchurMat, DSDPVec, DSDPVec);
extern int DSDPSchurMatMultR(DSDPSchurMat, DSDPVec, DSDPVec);
extern int DSDPSchurMatReducePVec(DSDPSchurMat, DSDPVec);
extern int DSDPSchurMatFactor(DSDPSchurMat,DSDPTruth*);
extern int DSDPSchurMatSolve(DSDPSchurMat, DSDPVec, DSDPVec);
extern int DSDPSchurMatDestroy(DSDPSchurMat*);
extern int DSDPSchurMatView(DSDPSchurMat);
extern int DSDPSchurMatSetR(DSDPSchurMat, double);

extern int DSDPSchurMatRowColumnScaling(DSDPSchurMat,int, DSDPVec,int*);
extern int DSDPSchurMatAddRow(DSDPSchurMat, int, double, DSDPVec);

extern int DSDPSchurMatVariableCompute(DSDPSchurMat, int, double*);
extern int DSDPSchurMatVariableComputeC(DSDPSchurMat, double*);
extern int DSDPSchurMatVariableComputeR(DSDPSchurMat, double*);
extern int DSDPSchurMatAddDiagonalElement(DSDPSchurMat, int, double);
extern int DSDPSchurMatAddC(DSDPSchurMat,int,double);
extern int DSDPSchurMatAddR(DSDPSchurMat,int,double);

extern int DSDPSchurMatDiagonalScaling(DSDPSchurMat, DSDPVec);
extern int DSDPSchurMatAddDiagonal(DSDPSchurMat, DSDPVec);

extern int DSDPSchurMatRowScaling(DSDPSchurMat, DSDPVec);

extern int DSDPZeroFixedVariables( DSDPSchurMat, DSDPVec);
extern int DSDPApplyFixedVariables( DSDPSchurMat, DSDPVec);
extern int DSDPIsFixed( DSDPSchurMat, int, DSDPTruth*);
extern int DSDPInitializeFixedVariable( FixedVariables *);
extern int DSDPAddFixedVariable( DSDPSchurMat, int, double);

#ifdef __cplusplus
}
#endif

#endif


