#ifndef __TAO_DSDPSDP_H
#define __TAO_DSDPSDP_H
/*! 
\file dsdpsdp.h
\brief Internal SDPCone data structures and routines
 */

#include "dsdpschurmat.h"
#include "dsdpbasictypes.h"
#include "dsdpvec.h"

#include "sdpconevec.h"
#include "dsdpdatamat.h"
#include "dsdpdualmat.h"
#include "dsdpxmat.h"
#include "dsdpdsmat.h"
#include "dsdplanczos.h"

typedef enum { SDPCONEEXIST=1, SDPCONESETUP1=2 } SDPConeStatus;

/*!
struct DSDPDataTranspose
\brief Internal structure for transpose of data.
*/
typedef struct{
  int m;
  int *nnzblocks;
  int **nzblocks;
  int *ttnzmat;
  int **nnzmats;
  int **idA;
  int *idAP;
} DSDPDataTranspose;

/*!
struct SDPBlockData
\brief Internal structure for data in one block of semidefintie.
*/
typedef struct{
  int maxnnzmats;
  int nnzmats;
  int *nzmat;
  DSDPDataMat *A;
  double r;        /* Multiple of Identity added to S to make it psd */
  double scl;
} DSDPBlockData;

/*!
struct SDPBlk
\brief Internal structure for block of semidefinite cone.
*/
typedef struct{

  DSDPBlockData ADATA;
  DSDPLanczosStepLength Lanczos;    /* For Lanczos steplength routine */

  int n;                                                /* Dimensions */
  double gammamu;          /* Scale Barrier, used only by user option */
  double bmu;                                    /* For LMI, not used */
  char format;       /* Packed Symmetric, Full Symmetric, Lower,Upper */
  int   nnz;
  SDPConeStatus status;

  SDPConeVec W;
  SDPConeVec W2;
  DSDPIndex IS;

  DSDPDualMat S;      /* Dual variable matrices */
  DSDPDualMat SS;     /* Compute primal variable matrices */
  DSDPDSMat DS;       /* Dual variable step matrices */
  DSDPVMat T;         /* Work Array and Primal variable matrice X   */

} SDPblk;

/*!
struct SDPCone_C 
\brief Internal structure for semidefinite cone.
\sa SDPCone
*/
struct SDPCone_C {
  int keyid;

  /* Dimensions */
  int m, nn;

  /* Data in block format */
  int nblocks;
  SDPblk *blk;

  /* Transpose of Data */
  DSDPDataTranspose ATR;

  /* Work space */
  DSDPVec Work, Work2;

  /* Current Solution */
  DSDPVec YY,YX,DYX;
  double xmakermu;

  int optype;
  DSDP dsdp;
};

#define SDPCONEKEY  5438
#define SDPConeValid(a) {if (!(a)||((a)->keyid!=SDPCONEKEY)){ DSDPSETERR(101,"DSDPERROR: Invalid SDPCone object\n");}}

#define DSDPCHKBLOCKERR(a,b);  { if (b){ DSDPSETERR1(b,"Block Number: %d,\n",a);} }
#define DSDPCHKVARERR(a,b);  { if (b){ DSDPSETERR1(b,"Variable Number: %d,\n",a);} }

extern int DSDPSetDataMatZero(DSDPDataMat*);

#include "dsdp5.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Operations on the Data */
extern int DSDPBlockDataInitialize(DSDPBlockData*);
extern int DSDPBlockDataAllocate(DSDPBlockData*, int);
extern int DSDPBlockAddDataMatrix(DSDPBlockData*,int, struct DSDPDataMat_Ops*, void*);
extern int DSDPBlockSetDataMatrix(DSDPBlockData*,int, struct DSDPDataMat_Ops*, void*);
extern int DSDPBlockRemoveDataMatrix(DSDPBlockData*,int);
extern int DSDPBlockDataMarkNonzeroMatrices(DSDPBlockData*,int*);
extern int DSDPBlockDataRowSparsity(DSDPBlockData*,int,int[],int[],int);
extern int DSDPBlockASum(DSDPBlockData*,double,DSDPVec,DSDPVMat);
extern int DSDPBlockADot(DSDPBlockData*,double,DSDPVec,DSDPVMat,DSDPVec);
extern int DSDPBlockvAv(DSDPBlockData*,double,DSDPVec,SDPConeVec,DSDPVec);
extern int DSDPBlockDataDestroy(DSDPBlockData*);
extern int DSDPBlockCheck(DSDPBlockData*,SDPConeVec,DSDPVMat);
extern int DSDPBlockANorm2(DSDPBlockData*, DSDPVec, int);
extern int DSDPBlockView(DSDPBlockData*);
extern int DSDPBlockView2(DSDPBlockData*);
extern int DSDPBlockCountNonzeroMatrices(DSDPBlockData *,int*);
extern int DSDPBlockGetMatrix(DSDPBlockData*,int,int*,double*,DSDPDataMat*);
extern int DSDPBlockFactorData(DSDPBlockData*,DSDPVMat,SDPConeVec);
extern int DSDPBlockTakeDownData(DSDPBlockData*);
extern int DSDPBlockDataRank(DSDPBlockData*,int*,int);

extern int DSDPBlockTakeDown(SDPblk*);
extern int DSDPBlockInitialize(SDPblk*);

extern int DSDPBlockEventInitialize(void);
extern int DSDPBlockEventZero(void);

extern int DSDPDataMatCheck(DSDPDataMat,SDPConeVec,DSDPIndex,DSDPVMat);

/* Operations on the Transpose of the Data */
extern int DSDPDataTransposeInitialize(DSDPDataTranspose*);
extern int DSDPDataTransposeTakeDown(DSDPDataTranspose*);
extern int DSDPDataTransposeSetup(DSDPDataTranspose*,SDPblk*,int,int);

extern int DSDPUseDefaultDualMatrix(SDPCone);

extern int SDPConeSetup(SDPCone,DSDPVec);
extern int SDPConeSetup2(SDPCone,DSDPVec,DSDPSchurMat);
extern int SDPConeComputeHessian(SDPCone,double,DSDPSchurMat,DSDPVec,DSDPVec);
extern int SDPConeMultiply(SDPCone,int,double,DSDPVec,DSDPVec,DSDPVec);
extern int SDPConeComputeRHS(SDPCone, int, double, DSDPVec, DSDPVec,DSDPVec);
extern int SDPConeComputeXX(SDPCone,int,DSDPVec,double,DSDPDualMat, DSDPVMat);
extern int SDPConeDestroy(SDPCone);

extern int SDPConeCheckJ(SDPCone,int);
extern int SDPConeCheckN(SDPCone,int, int);
extern int SDPConeCheckM(SDPCone,int);
extern int SDPConeCheckStorageFormat(SDPCone,int, char);


extern int SDPConeComputeSS(SDPCone, int, DSDPVec, DSDPVMat);
extern int SDPConeComputeXDot(SDPCone,int,DSDPVec,DSDPVMat,DSDPVec,double*,double*, double *);
extern int SDPConeComputeX3(SDPCone,int,double,DSDPVec,DSDPVec,DSDPVMat);

/* extern int DSDPPrintSDPA(TAO_DSDP *); */
extern int DSDPMakeVMat(char, int, DSDPVMat*);
extern int DSDPMakeVMatWithArray(char, double[], int, int, DSDPVMat*);

extern int DSDPSetDualMatrix(SDPCone sdpcone,int (*createdualmatrix)(DSDPBlockData*,DSDPVec,DSDPVMat,DSDPVec,DSDPVec,DSDPDualMat*,DSDPDualMat*,DSDPDSMat*,void*),void*);

extern int DSDPAddSDP(DSDP,SDPCone);
#ifdef __cplusplus
}
#endif

extern int SDPConeSetRIdentity(SDPCone,int,int,double);
extern int DSDPDualMatEventInitialize(void);
extern int DSDPVMatEventInitialize(void);
extern int DSDPDualMatEventZero(void);
extern int DSDPVMatEventZero(void);

#endif
