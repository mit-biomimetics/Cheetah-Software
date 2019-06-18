#ifndef __DSDP_H
#define __DSDP_H

#include "dsdpbasictypes.h"
#include "dsdpvec.h"
#include "dsdpschurmat.h"
#include "dsdpcone.h"
#include "dsdpconverge.h"
/*!
  \file dsdp.h
  \brief Internal data structure for the DSDP solver
*/

/*!
\typedef struct LUBounds_C* YBoundCone;
\brief Cone with bounds on variables y.
*/
typedef struct LUBounds_C* YBoundCone;

/*!
\typedef typedef struct RDCone* RRCone;
\brief Cone with nonnegativity on variable r.
*/
typedef struct RDCone* RRCone;


#define MAX_DSDP_MONITORS 5
#define MAX_XMAKERS 4

typedef struct {
  DSDPVec y;
  DSDPVec dy;
  double  mu;
  double  pstep;
  DSDPVec rhs;
} XMaker;

typedef struct { /* This information is needed to compute the step Direction */
  DSDPVec y;
  double zbar;
  double mutarget;
  double logdet;
} DSDPState;

typedef struct {
  int (*f)(void*);
  void * ptr;
}DRoutine;

typedef struct {
  int (*monitor)(struct DSDP_C *, void*);
  void *monitorctx;
} DMonitor;

typedef struct {
  DSDPCone cone;
  int coneid;
} DCone;

/*!
\struct DSDP_C dsdp.h "src/solver/dsdp.h"
\brief Internal structures for the DSDP solver.
\sa DSDP
*/
struct DSDP_C{

  DSDPCG *sles;
  int slestype;

  double schurmu;
  DSDPSchurMat M;
  double Mshift;
  double maxschurshift;

  int ncones,maxcones;
  DCone* K;

  int keyid;

  int solvetime,cgtime,ptime,dtime,ctime;
  int reuseM;
  DSDPTruth goty0;
  DSDPTruth setupcalled;

  int    m;       /*     number of constraints                        */ 
  double np;      /*     Dimension of full variable matrix            */

  int    itnow;   /* current iterate                                  */
  int    maxiter; /* Maximum number of iterates                       */
  double pobj;    /* current primal objective value - use duality gap */
  double ppobj;   /* current primal objetive value - evaluate P       */
  double dobj,ddobj;  /* the current dual objective value             */
  double pstep,dstep;  /* current primal and dual step lengths        */
  double dualitygap;
  double mutarget;
  double mu,muold,mu0;      /* The current  mu */
  double rho,potential,logdet,rhon;
  double pnorm;   /* the current value of ||P||                       */
  double maxtrustradius;
  double cnorm,anorm,bnorm; 
  double tracex,tracexs;
  double rgap;
  double pstepold;

  DSDPVec y;      /* dual variables                                       */
  DSDPVec y0;     /* initial dual variables                               */
  DSDPVec ytemp;  /* temporary dual variables                             */
  DSDPVec dy1;    /* search direction 1 affine direction                  */
  DSDPVec dy2;    /* search direction 2 centering direction               */
  DSDPVec dy;     /* total search direction = constant*dy1-dy2            */
  DSDPVec rhs1;   /* objective vector b to determine step direction       */
  DSDPVec rhs2;   /* barrier vector A(S^{-1}) to determine step direction */
  DSDPVec rhs;    /* right-hand side of linear system                     */
  DSDPVec rhstemp;/* temporary rhs vector                                 */
  DSDPVec b;      /* dual objective vector                                */

  /* Multiple of identity matrix added to dual        */
  double r;
  int rflag;
  DSDPPenalty UsePenalty;
  RRCone rcone;

  DSDPTruth usefixedrho; /*   True if fixed rho used. */
  
  XMaker xmaker[MAX_XMAKERS]; /* step direction used to create X */
  DSDPVec xmakerrhs;

  YBoundCone ybcone;
  double pinfeas; /* Infeasible in P indirectly -- neglect numerical errors */
  double perror;  /* Infeasible in P computed directly   */

  DSDPSolutionType pdfeasible;
  double dinfeastol; /* Parameter: Classify (D) as feasible */
  double pinfeastol; /* Parameter: Classify (P) as feasible */

  ConvergenceMonitor conv;
  DSDPTerminationReason reason;

  DMonitor dmonitor[MAX_DSDP_MONITORS];
  int nmonitors;

  DRoutine droutine[10];
  int ndroutines;
};

typedef struct DSDP_C PD_DSDP;

#define DSDPKEY  5432

#define DSDPValid(a) {if (!(a)||((a)->keyid!=DSDPKEY)){ DSDPSETERR(101,"DSDPERROR: Invalid DSDP object\n");}}

#include "dsdpbasictypes.h"


extern int DSDPCreateLUBoundsCone(DSDP, YBoundCone*);
extern int BoundYConeSetBounds(YBoundCone, double, double);
extern int BoundYConeGetBounds(YBoundCone, double*, double*);
extern int BoundYConeAddX(YBoundCone,double,DSDPVec,DSDPVec,DSDPVec,double*);
extern int BoundYConeAddS(YBoundCone,DSDPVec,DSDPVec);

#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPComputeObjective(DSDP, DSDPVec, double *);
extern int DSDPComputeDY(DSDP, double, DSDPVec, double*);
extern int DSDPComputeNewY(DSDP, double, DSDPVec);
extern int DSDPComputeRHS(DSDP, double, DSDPVec);
extern int DSDPComputePDY1(DSDP,double,DSDPVec);
extern int DSDPComputePDY(DSDP, double, DSDPVec, double*);
extern int DSDPComputePY(DSDP,double,DSDPVec);
extern int DSDPComputeG(DSDP, DSDPVec,DSDPVec,DSDPVec);
extern int DSDPComputePNorm(DSDP, double,DSDPVec,double*);
extern int DSDPComputeDualityGap(DSDP, double, double*);
extern int DSDPSaveYForX(DSDP, double, double);
extern int DSDPSetY(DSDP,double,double,DSDPVec);
extern int DSDPComputePotential(DSDP, DSDPVec, double,double*);
extern int DSDPComputePotential2(DSDP, DSDPVec, double, double, double*);
extern int DSDPSetRR(DSDP,double);
extern int DSDPGetRR(DSDP,double*);
extern int DSDPGetMaxYElement(DSDP,double*);

extern int DSDPSolveDynamicRho(DSDP);
extern int DSDPRefineStepDirection(DSDP,DSDPVec, DSDPVec);

/* Cone operations */
extern int DSDPSetUpCones(DSDP);
extern int DSDPSetUpCones2(DSDP, DSDPVec, DSDPSchurMat);
extern int DSDPSchurSparsity(DSDP, int, int[], int);
extern int DSDPComputeSS(DSDP, DSDPVec, DSDPDualFactorMatrix, DSDPTruth*);
extern int DSDPInvertS(DSDP);
extern int DSDPComputeHessian( DSDP, DSDPSchurMat,  DSDPVec, DSDPVec);
extern int DSDPHessianMultiplyAdd( DSDP, DSDPVec, DSDPVec);
extern int DSDPComputeMaxStepLength(DSDP, DSDPVec, DSDPDualFactorMatrix, double*);
extern int DSDPComputeLogSDeterminant(DSDP, double*);
extern int DSDPComputeANorm2(DSDP,DSDPVec);
extern int DSDPViewCones(DSDP);
extern int DSDPMonitorCones(DSDP,int);
extern int DSDPDestroyCones(DSDP);
extern int DSDPPassXVectors(DSDP,double,DSDPVec,DSDPVec);
extern int DSDPComputeXVariables(DSDP,double,DSDPVec,DSDPVec,DSDPVec,double*);
extern int DSDPGetConicDimension(DSDP,double*);
extern int DSDPSchurSparsity(DSDP, int, int[], int);

extern int DSDPCGSolve(DSDP,DSDPSchurMat,DSDPVec,DSDPVec,double,DSDPTruth*);
extern int DSDPComputeDualStepDirection(DSDP, double, DSDPVec);
extern int DSDPComputeDualStepDirections(DSDP);

extern int DSDPSetCone(DSDP,DSDPCone);

extern int DSDPScaleData(DSDP);
extern int DSDPComputeDataNorms(DSDP);

extern int DSDPTakeDown(DSDP);
extern int DSDPSetDefaultStatistics(DSDP);
extern int DSDPSetDefaultParameters(DSDP);
extern int DSDPSetDefaultMonitors(DSDP);

/* DSDP subroutines */
extern int DSDPInitializeVariables(DSDP);
extern int DSDPObjectiveGH( DSDP, DSDPSchurMat, DSDPVec);
extern int DSDPCheckForUnboundedObjective(DSDP, DSDPTruth*);

extern int DSDPCheckConvergence(DSDP,DSDPTerminationReason *);
extern int DSDPMonitor(DSDP);
extern int DSDPPrintStats(DSDP, void*); 
extern int DSDPPrintStatsFile(DSDP, void*);

extern int DSDPGetSchurMatrix(DSDP,DSDPSchurMat*);


#ifdef __cplusplus
}
#endif

extern int DSDPAddRCone(DSDP,RRCone*);
extern int RConeSetType(RRCone, DSDPPenalty);
extern int RConeGetRX(RRCone, double*);

extern int DSDPGetConvergenceMonitor(DSDP, ConvergenceMonitor**);
extern int DSDPDefaultConvergence(DSDP,void *);

extern int DSDPSetDefaultSchurMatrixStructure(DSDP);
extern int DSDPFixedVariablesNorm(DSDPSchurMat,DSDPVec);
extern int DSDPComputeFixedYX(DSDPSchurMat,DSDPVec);

extern int DSDPAddBCone(DSDP,DSDPVec,double);

#endif
