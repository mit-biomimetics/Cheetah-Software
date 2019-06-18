#if !defined(__DSDP_DSMATRIXOPS_H) 
#define __DSDP_DSMATRIXOPS_H

/*!
\file dsdpdsmat_impl.h
\brief Structure of function pointers that each SDP Delta S matrix type 
(sparse, dense, diagonal, ...) must implement.
 */

/*!
struct DSDPDSMat_Ops

\brief Symmetric Delta S matrix for one block in the semidefinite cone.
*/
struct  DSDPDSMat_Ops{
  int id;
  int (*matzeroentries)(void*); 
  int (*matmult)(void*,double[],double[], int); /* Multiply by a vector */
  int (*matgetsize)(void*,int*);
  int (*matseturmat)(void*,double[],int,int); /* Set values from array */
  int (*matvecvec)(void*,double[],int,double*); /* v' * DS * v */
  int (*mattest)(void*);
  int (*matview)(void*);
  int (*matdestroy)(void*);
  const char *matname;
};

#ifdef __cplusplus
extern "C" {
#endif
extern int DSDPDSMatOpsInitialize(struct  DSDPDSMat_Ops*);
#ifdef __cplusplus
}
#endif

#endif


