#include "dsdpbasictypes.h"
#if !defined(__DSDP_SCHURMATRICES_H) 
#define __DSDP_SCHURMATRICES_H

/*!
\file dsdpschurmat_impl.h
\brief Function pointers that a Schur complement matrix (dense, sparse, parallel dense) must provide.
*/

struct  DSDPSchurMat_Ops{
  int id;
  int (*matzero)(void*);
  int (*matrownonzeros)(void*,int,double*,int*,int);
  int (*mataddrow)(void*,int,double,double[],int);
  int (*mataddelement)(void*,int,double);
  int (*matadddiagonal)(void*,double[],int);
  int (*matshiftdiagonal)(void*,double);
  int (*matassemble)(void*);
  int (*matscaledmultiply)(void*,double[],double[],int);
  int (*matmultr)(void*,double[],double[],int);
  int (*matfactor)(void*,int*);
  int (*matsolve)(void*,double[],double[],int);
  int (*matsetup)(void*,int);
  int (*pmatwhichdiag)(void*,double[],int);
  int (*pmatonprocessor)(void*,int,int*);
  int (*pmatlocalvariables)(void*,double[],int);
  int (*pmatreduction)(void*,double[],int);
  int (*pmatdistributed)(void*,int*);
  int (*matdestroy)(void*);
  int (*matview)(void*);
  const char *matname;
};

extern int DSDPSetSchurMatOps(DSDP,struct DSDPSchurMat_Ops*,void*);
extern int DSDPSchurMatOpsInitialize(struct DSDPSchurMat_Ops*);
extern int DSDPSparsityInSchurMat(DSDP,int,int[],int);
#endif


