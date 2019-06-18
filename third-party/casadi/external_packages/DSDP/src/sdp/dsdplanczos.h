#ifndef __TAO_DSDPSTEP_H
#define __TAO_DSDPSTEP_H
#include "sdpconevec.h"
/*! 
\file dsdplanczos.h
\brief Lanczos procedure determines the maximum step length
 */

/*!
typedef struct DSDPLanczosStepLength
\brief Apply Lanczos prodedure to find distance to boundary.
*/
typedef struct {
  int lanczosm;
  int maxlanczosm;
  double *darray;  /* For SLOW BUT ROBUST METHOD ONLY */
  SDPConeVec *Q; /* Size 2 for FAST, Size lanczosm for SLOW */
  SDPConeVec Tv; /* For SLOW BUT ROBUST METHOD ONLY */

  double *dwork4n;
  int *iwork10n;
  int lwork,liwork;
  int n;
  int type;
} DSDPLanczosStepLength;

#ifdef __cplusplus
extern "C" {
#endif
extern int DSDPLanczosInitialize(DSDPLanczosStepLength*);
extern int DSDPSetMaximumLanczosIterations( DSDPLanczosStepLength *LZ, int);
extern int DSDPFastLanczosSetup(DSDPLanczosStepLength*,SDPConeVec);
extern int DSDPRobustLanczosSetup(DSDPLanczosStepLength*,SDPConeVec);
extern int DSDPLanczosStepSize( DSDPLanczosStepLength*, SDPConeVec, SDPConeVec, DSDPDualMat, DSDPDSMat, double *);
extern int DSDPLanczosDestroy( DSDPLanczosStepLength*);
extern int DSDPLanczosMinXEig( DSDPLanczosStepLength*, DSDPVMat, SDPConeVec, SDPConeVec, double *);

#ifdef __cplusplus
}
#endif

#endif
