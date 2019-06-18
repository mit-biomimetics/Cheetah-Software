#ifndef __DSDPCONE_H
#define __DSDPCONE_H

/*!
\file dsdpcone_impl.h
\brief Implementations of a cone (SDP,LP,...) must provide a structure of function pointers.
*/

#include "dsdpbasictypes.h"
#include "dsdpvec.h"
#include "dsdpschurmat.h"

struct  DSDPCone_Ops{
  int id;
  int (*conesize)(void*,double*);
  int (*conesetup)(void*,DSDPVec);
  int (*conesetup2)(void*,DSDPVec,DSDPSchurMat);
  int (*conecomputes)(void*,DSDPVec,DSDPDualFactorMatrix,DSDPTruth*);
  int (*coneinverts)(void*);
  int (*conelogpotential)(void*,double*,double*);
  int (*conesetxmaker)(void*,double,DSDPVec,DSDPVec);
  int (*conecomputex)(void*,double,DSDPVec,DSDPVec,DSDPVec,double*);
  int (*conehessian)(void*,double,DSDPSchurMat,DSDPVec,DSDPVec);
  int (*conehmultiplyadd)(void*,double,DSDPVec,DSDPVec,DSDPVec);
  int (*conerhs)(void*,double,DSDPVec,DSDPVec,DSDPVec);
  int (*conemaxsteplength)(void*,DSDPVec,DSDPDualFactorMatrix,double*);
  int (*coneanorm2)(void*,DSDPVec);
  int (*conesparsity)(void*,int,int*,int[],int);
  int (*conemonitor)(void*,int);
  int (*conedestroy)(void*);
  int (*coneview)(void*);
  const char *name;
};

extern int DSDPAddCone(DSDP,struct  DSDPCone_Ops*, void*);
extern int DSDPConeOpsInitialize(struct  DSDPCone_Ops*);

#endif
