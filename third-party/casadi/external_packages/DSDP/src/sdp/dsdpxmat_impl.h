#if !defined(__DSDP_VMATRIXOPS_H) 
#define __DSDP_VMATRIXOPS_H
/*!
\file dsdpxmat_impl.h
\brief Structure of function pointers that each dense matrix array type 
(upper full, packed symmetric, ...) must implement.
 */

/*!
struct DSDPVMat_Ops
\brief Table of function pointers that operate on the dense matrix.
*/
struct  DSDPVMat_Ops{
  int id;
  int (*matgetsize)(void*,int*);
  int (*mataddouterproduct)(void*,double,double[],int);
  int (*matmult)(void*,double[],double[],int);
  int (*matscalediagonal)(void*,double);
  int (*matshiftdiagonal)(void*,double);
  int (*matfnorm2)(void*,int,double*);
  int (*matzeroentries)(void*);
  int (*matgeturarray)(void*,double*[],int*);
  int (*matrestoreurarray)(void*,double*[],int*);
  int (*matmineig)(void*,double[],double[],int,double*);
  int (*mattest)(void*);
  int (*matdestroy)(void*);
  int (*matview)(void*);
  const char *matname;

};

#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPVMatOpsInitialize(struct  DSDPVMat_Ops*);

#ifdef __cplusplus
}
#endif

#endif


