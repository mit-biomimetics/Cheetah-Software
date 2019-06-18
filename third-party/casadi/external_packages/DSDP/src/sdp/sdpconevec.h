#if !defined(__SDPCONE_VECTORS_H) 
#define __SDPCONE_VECTORS_H
/*!
\file sdpconevec.h
\brief Each block of the SDPCone has two vectors of appropriate size.
*/
#include <math.h>
/*!
struct SDPConeVec_C
\brief Vector whose length corresponds to dimension of a block in a cone.
\sa SDPConeVec
*/
struct  SDPConeVec_C{
  int    dim;
  double *val;
};

typedef struct {
  int *indx;
} DSDPIndex;

/*!
\typedef struct SDPConeVec_C   SDPConeVec;
\brief SDPConeVec is a vector with the dimension of the block in the SDP cone.
*/
typedef struct SDPConeVec_C SDPConeVec;

#ifdef __cplusplus
extern "C" {
#endif

extern int SDPConeVecCreate(int,SDPConeVec *);
extern int SDPConeVecDuplicate(SDPConeVec,SDPConeVec *);
extern int SDPConeVecDot(SDPConeVec, SDPConeVec, double *);
extern int SDPConeVecView( SDPConeVec);
extern int SDPConeVecDestroy(SDPConeVec*);
#define SDPConeVecCreateWArray(a,b,c)        0;{ (*(a)).val=(b); (*(a)).dim=(c);}

extern int SDPConeVecSet(double, SDPConeVec );
extern int SDPConeVecZero(SDPConeVec );
extern int SDPConeVecNormalize(SDPConeVec );
extern int SDPConeVecAXPY(double,  SDPConeVec,  SDPConeVec);
extern int SDPConeVecNorm2( SDPConeVec, double *);
extern int SDPConeVecCopy( SDPConeVec,  SDPConeVec);
extern int SDPConeVecScale(double, SDPConeVec);

#define SDPConeVecGetArray(a,b)              0;{ *(b)=((a).val); }
#define SDPConeVecRestoreArray(a,b)          0;{ *(b)=0;}
#define SDPConeVecGetSize(a,b)               0;{ *(b)=((a).dim); }

/*!
typedef struct { int *indx; } DSDPIndex;

\brief Identifies sparsity in SDPConeVec.
*/

extern int DSDPIndexInitialize(DSDPIndex*);
extern int DSDPIndexDestroy(DSDPIndex*);
extern int DSDPIndexSetBasis(DSDPIndex, int);
extern int DSDPIndexCreate(int,DSDPIndex*);
extern int DSDPIndexView(DSDPIndex);
#ifdef __cplusplus
}
#endif


#endif
