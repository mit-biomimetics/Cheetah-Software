#include "dsdpvec.h"
#if !defined(__DSDP_CG_H) 
#define __DSDP_CG_H
/*!
  \file dsdpcg.h
  \brief Internal data structure for CG method.
*/
typedef struct{

  int    setup2;
  int    m;

  DSDPVec    Diag;
  DSDPVec    RHS2;
  DSDPVec    R;
  DSDPVec    BR;
  DSDPVec    P;

  DSDPVec    BP;
  DSDPVec    TTT;

} DSDPCG;


#ifdef __cplusplus
extern "C" {
#endif
int DSDPCGSetup(DSDPCG*, DSDPVec);
int DSDPCGDestroy(DSDPCG**);
int DSDPCGInitialize(DSDPCG **);
#ifdef __cplusplus
}
#endif

#endif
