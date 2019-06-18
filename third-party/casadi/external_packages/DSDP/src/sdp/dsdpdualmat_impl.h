#if !defined(__DSDP_DUALMATRIXOPS_H) 
#define __DSDP_DUALMATRIXOPS_H

/*!
\file dsdpdualmat_impl.h
\brief Structure of function pointers that each symmetric positive definite matrix type (dense, sparse) must implement.
*/

/*!
struct DSDPDualMat_Ops
\brief Table of function pointers that operate on the S matrix.
*/
struct  DSDPDualMat_Ops{
  int id;
  int (*matseturmat)(void*,double[],int,int); /* Set full array of values into matrix */
  int (*matgetarray)(void*,double*[],int*);   /* Set pointers to a dense array and its dimension */
  int (*matcholesky)(void*,int*);  /* Cholesky second argument is 0 for success, 1 otherwise */
  int (*matsolveforward)(void*,double[],double[],int); /* Called after Cholesky */
  int (*matsolvebackward)(void*,double[],double[],int); /* Called after Cholesky */
  int (*matinvert)(void*); /* Called after Cholesky factorization */
  int (*matinverseadd)(void*,double,double[],int,int); /* Add multiple of the inverse to array */
  int (*matinversemultiply)(void*,int[],int,double[],double[],int); /* Called after invert */
  int (*matforwardmultiply)(void*,double[],double[],int); /* Called after invert */
  int (*matbackwardmultiply)(void*,double[],double[],int); /* Called after invert */
  int (*matlogdet)(void*,double*); /* Called after Cholesky */
  int (*matfull)(void*,int*); /* Is fully dense or not? */
  int (*mattest)(void*);
  int (*matgetsize)(void*,int*);
  int (*matdestroy)(void*);
  int (*matview)(void*);
  const char *matname;
};

#ifdef __cplusplus
extern "C" {
#endif
extern int DSDPDualMatOpsInitialize(struct  DSDPDualMat_Ops*);
#ifdef __cplusplus
}
#endif

#endif


