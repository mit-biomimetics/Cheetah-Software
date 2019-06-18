#include "dsdp.h"
#include "dsdpsys.h"
#include "dsdp5.h"
/*!
  \file dsdpsetup.c
  \brief Create DSDP solver and its data strucutures.
*/

/*!
\fn int DSDPCreate(int m, DSDP *dsdpnew)
\brief Create a DSDP solver. FIRST DSDP routine!

\param m the number of variables y
\param *dsdpnew will be set to a new solver object
\sa DSDPSetDualObjective()
\sa DSDPSetup()
\sa DSDPSolve()
\sa DSDPDestroy()
\ingroup DSDPBasic

For example, to create a DSDP solver for a problem with 10 y variables,
\code
int m=10;
DSDP dsdp;
DSDPCreate(m,&dsdp);
\endcode
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPCreate"
int DSDPCreate(int m,DSDP* dsdpnew){

  DSDP dsdp;
  int info;

  DSDPFunctionBegin;

  DSDPCALLOC1(&dsdp,PD_DSDP,&info);DSDPCHKERR(info);
  *dsdpnew=dsdp;
  dsdp->keyid=DSDPKEY;

  /* Initialize some parameters */
  DSDPEventLogInitialize();
  dsdp->m=m;
  dsdp->maxcones=0;
  dsdp->ncones=0;
  dsdp->K=0;
  dsdp->setupcalled=DSDP_FALSE;
  dsdp->ybcone=0;
  dsdp->ndroutines=0;
  /*  info = DSDPSetStandardMonitor(dsdp);DSDPCHKERR(info);  */
  info = DSDPVecCreateSeq(m+2,&dsdp->b);DSDPCHKERR(info);
  info = DSDPVecZero(dsdp->b);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->b,&dsdp->y);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->b,&dsdp->ytemp);DSDPCHKERR(info);
  info = DSDPVecZero(dsdp->y);DSDPCHKERR(info);
  info = DSDPVecSetC(dsdp->y,-1.0);DSDPCHKERR(info);

  info = DSDPAddRCone(dsdp,&dsdp->rcone);DSDPCHKERR(info);
  info = DSDPCreateLUBoundsCone(dsdp,&dsdp->ybcone);DSDPCHKERR(info);

  info=DSDPSetDefaultStatistics(dsdp);DSDPCHKERR(info);
  info=DSDPSetDefaultParameters(dsdp);DSDPCHKERR(info);
  info=DSDPSetDefaultMonitors(dsdp);DSDPCHKERR(info);

  /*  info = DSDPMatInitialize(m,m,&dsdp->Q);DSDPCHKERR(info); */
  info = DSDPSchurMatInitialize(&dsdp->M);DSDPCHKERR(info);
  info = DSDPSetDefaultSchurMatrixStructure(dsdp); DSDPCHKERR(info);
  info = DSDPCGInitialize(&dsdp->sles); DSDPCHKERR(info);

  /* Set the one global variable
     sdat=dsdp;
  */
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDefaultStatistics"
/*! 
\fn int DSDPSetDefaultStatistics(DSDP dsdp);
\brief Set default statistics.
\param dsdp the solver
*/
int DSDPSetDefaultStatistics(DSDP dsdp){
  
  int i;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  dsdp->reason=CONTINUE_ITERATING;
  dsdp->pdfeasible=DSDP_PDUNKNOWN;
  dsdp->itnow=0;
  dsdp->pobj=  1.0e10;
  dsdp->ppobj=  1.0e10;
  dsdp->dobj= -1.0e+9;
  dsdp->ddobj= -1.0e+9;
  dsdp->dualitygap=dsdp->ppobj-dsdp->ddobj;
  dsdp->pstep=1.0;
  dsdp->dstep=0.0;
  for (i=0;i<MAX_XMAKERS;i++){
    dsdp->xmaker[i].mu=1.0e200;
    dsdp->xmaker[i].pstep=0.0;
  }
  dsdp->pnorm=0.001;
  dsdp->mu=1000.0;
  dsdp->np=0;
  dsdp->anorm=0;
  dsdp->bnorm=0;
  dsdp->cnorm=0;
  dsdp->tracex=0;
  dsdp->tracexs=0;
  dsdp->Mshift=0;
  dsdp->goty0=DSDP_FALSE;
  DSDPFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDefaultParameters"
/*! 
\fn int DSDPSetDefaultParameters(DSDP dsdp);
\brief Set default parameters.
\param dsdp the solver
*/
int DSDPSetDefaultParameters(DSDP dsdp){
  
  int info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);

  /* Stopping parameters */
  info=DSDPSetMaxIts(dsdp,500);DSDPCHKERR(info);
  info=DSDPSetGapTolerance(dsdp,1.0e-6);DSDPCHKERR(info);
  info=DSDPSetPNormTolerance(dsdp,1.0e30);DSDPCHKERR(info);
  if (dsdp->m<100){info=DSDPSetGapTolerance(dsdp,1.0e-7);DSDPCHKERR(info);}
  if (dsdp->m>3000){info=DSDPSetGapTolerance(dsdp,5.0e-6);DSDPCHKERR(info);}
  info=RConeSetType(dsdp->rcone,DSDPInfeasible);DSDPCHKERR(info);
  info=DSDPSetDualBound(dsdp,1.0e20);DSDPCHKERR(info);
  info=DSDPSetStepTolerance(dsdp,5.0e-2);DSDPCHKERR(info);
  info=DSDPSetRTolerance(dsdp,1.0e-6);DSDPCHKERR(info);
  info=DSDPSetPTolerance(dsdp,1.0e-4);DSDPCHKERR(info);
  /* Solver options */
  info=DSDPSetMaxTrustRadius(dsdp,1.0e10);DSDPCHKERR(info);
  info=DSDPUsePenalty(dsdp,0);DSDPCHKERR(info);
  info=DSDPSetInitialBarrierParameter(dsdp,-1.0);DSDPCHKERR(info);
  info=DSDPSetPotentialParameter(dsdp,3.0);DSDPCHKERR(info);
  info=DSDPUseDynamicRho(dsdp,1);DSDPCHKERR(info);
  info=DSDPSetR0(dsdp,-1.0);DSDPCHKERR(info);
  info=DSDPSetPenaltyParameter(dsdp,1.0e8);DSDPCHKERR(info);
  info=DSDPReuseMatrix(dsdp,4);DSDPCHKERR(info);
  if (dsdp->m>100){info=DSDPReuseMatrix(dsdp,7);DSDPCHKERR(info);}
  if (dsdp->m>1000){info=DSDPReuseMatrix(dsdp,10);DSDPCHKERR(info);}
  if (dsdp->m<=100){info=DSDPSetPotentialParameter(dsdp,5.0);DSDPCHKERR(info);DSDPCHKERR(info);}
  dsdp->maxschurshift=1.0e-11;
  dsdp->mu0=-1.0;
  dsdp->slestype=2;
  info = DSDPSetYBounds(dsdp,-1e7,1e7);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDefaultMonitors"
/*! 
\fn int DSDPSetDefaultMonitors(DSDP dsdp);
\brief Set convergence monitor.
\param dsdp the solver
*/
int DSDPSetDefaultMonitors(DSDP dsdp){
  
  int info;

  DSDPFunctionBegin;
  DSDPValid(dsdp);
  dsdp->nmonitors=0;
  info=DSDPSetMonitor(dsdp,DSDPDefaultConvergence,(void*)&dsdp->conv); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetup(DSDP dsdp)
\brief Set up data structures in the solver and the cones associated with it.

\param dsdp the solver
\sa DSDPCreate()
\sa DSDPSolve()
\sa DSDPDestroy()
\ingroup DSDPBasic

This routine must be called before DSDPSolve().  Do not create
SDP, LP or other cones after calling this routines, and do not
set data into the cones after calling this routine.

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetUp"
int DSDPSetup(DSDP dsdp){
  
  int i,info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  
  /* Create the Work Vectors */
  info = DSDPVecDuplicate(dsdp->y,&dsdp->rhs1);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->rhs2);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->rhs);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->rhstemp);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->dy1);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->dy2);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->dy);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->y0);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->xmakerrhs);DSDPCHKERR(info);
  for (i=0;i<MAX_XMAKERS;i++){
    info = DSDPVecDuplicate(dsdp->y,&dsdp->xmaker[i].y);DSDPCHKERR(info);
    info = DSDPVecDuplicate(dsdp->y,&dsdp->xmaker[i].dy);DSDPCHKERR(info);
    info = DSDPVecDuplicate(dsdp->y,&dsdp->xmaker[i].rhs);DSDPCHKERR(info);
  }

  /* Create M */
  info = DSDPSetUpCones(dsdp);DSDPCHKERR(info);
  info = DSDPSchurMatSetup(dsdp->M,dsdp->ytemp);DSDPCHKERR(info); 

  info = DSDPCGSetup(dsdp->sles,dsdp->ytemp); DSDPCHKERR(info);

  info = DSDPSetUpCones2(dsdp,dsdp->y,dsdp->M);DSDPCHKERR(info);
  info = DSDPGetConicDimension(dsdp,&dsdp->np);DSDPCHKERR(info);

  info=DSDPComputeDataNorms(dsdp);DSDPCHKERR(info);
  dsdp->pinfeas=dsdp->bnorm+1;
  dsdp->perror=dsdp->bnorm+1;
  info=DSDPScaleData(dsdp);DSDPCHKERR(info);

  info=DSDPGetConicDimension(dsdp,&dsdp->np);DSDPCHKERR(info);
  dsdp->solvetime=0;
  dsdp->cgtime=0;
  dsdp->ptime=0;
  dsdp->dtime=0;
  dsdp->ctime=0;
  info=DSDPEventLogRegister("Primal Step",&dsdp->ptime);
  info=DSDPEventLogRegister("Dual Step",&dsdp->dtime);
  info=DSDPEventLogRegister("Corrector Step",&dsdp->ctime);
  info=DSDPEventLogRegister("CG Solve",&dsdp->cgtime);
  info=DSDPEventLogRegister("DSDP Solve",&dsdp->solvetime);
  dsdp->setupcalled=DSDP_TRUE;
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DSDPGetSchurMatrix"
int DSDPGetSchurMatrix(DSDP dsdp, DSDPSchurMat *M){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *M=dsdp->M;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetConvergenceMonitor"
/*!
\fn int DSDPGetConvergenceMonitor(DSDP dsdp, ConvergenceMonitor**ctx);
\brief Get the structure containing convergence parameters.
\param dsdp the solver
\param *ctx will point to the structure.
\ingroup DSDPRoutines

\note 
This structure part of the DSDP structure.

 */
int DSDPGetConvergenceMonitor(DSDP dsdp, ConvergenceMonitor**ctx){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *ctx=&dsdp->conv;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeDataNorms"
/*! 
\fn int DSDPComputeDataNorms(DSDP dsdp);
\brief Compute norms of A,C, and b.
\param dsdp the solver
*/
int DSDPComputeDataNorms(DSDP dsdp){
  int info;
  DSDPVec ytemp=dsdp->ytemp;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info = DSDPComputeANorm2(dsdp,ytemp);DSDPCHKERR(info);
  info = DSDPFixedVariablesNorm(dsdp->M,ytemp);DSDPCHKERR(info);
  info = DSDPVecGetC(ytemp,&dsdp->cnorm);DSDPCHKERR(info);
  dsdp->cnorm=sqrt(dsdp->cnorm);
  info = DSDPVecSetR(ytemp,0);DSDPCHKERR(info);
  info = DSDPVecSetC(ytemp,0);DSDPCHKERR(info);
  info = DSDPVecNorm1(ytemp,&dsdp->anorm);DSDPCHKERR(info);
  dsdp->anorm=sqrt(dsdp->anorm);
  DSDPLogInfo(0,2,"Norm of data: %4.2e\n",dsdp->anorm);
  info=DSDPVecCopy(dsdp->b,ytemp);DSDPCHKERR(info);
  info = DSDPVecSetR(ytemp,0);DSDPCHKERR(info);
  info = DSDPVecSetC(ytemp,0);DSDPCHKERR(info);
  info = DSDPVecNorm2(ytemp,&dsdp->bnorm);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPScaleData"
/*! 
\fn int DSDPScaleData(DSDP dsdp);
\brief Scale the matrix C.
\param dsdp the solver
*/
int DSDPScaleData(DSDP dsdp){
  int info;
  double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  scale=1.0*dsdp->anorm;
  if (dsdp->bnorm){ scale/=dsdp->bnorm;}
  if (dsdp->cnorm){ scale/=dsdp->cnorm;}
  scale=DSDPMin(scale,1.0);
  scale=DSDPMax(scale,1.0e-6);
  if (dsdp->cnorm==0){  scale=1;}
  info=DSDPSetScale(dsdp,scale);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSolve(DSDP dsdp)
\brief Apply DSDP to the problem.

Call this routine after DSDPCreate() and DSDPSetup(),
and after setting the data.

\param dsdp is the solver
\sa DSDPCreate()
\sa DSDPGetSolutionType()
\sa DSDPGetDObjective()
\sa DSDPGetY()
\sa DSDPStopReason()
\ingroup DSDPBasic
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSolve"
int DSDPSolve(DSDP dsdp){
  int info;
  DSDPFunctionBegin;
  info=DSDPEventLogBegin(dsdp->solvetime);
  dsdp->pdfeasible=DSDP_PDUNKNOWN;
  info=DSDPSetConvergenceFlag(dsdp,CONTINUE_ITERATING);DSDPCHKERR(info); 
  info=DSDPInitializeVariables(dsdp);DSDPCHKERR(info);
  info=DSDPSolveDynamicRho(dsdp);DSDPCHKERR(info);
  if (dsdp->pstep==1){info=DSDPRefineStepDirection(dsdp,dsdp->xmakerrhs,dsdp->xmaker[0].dy);DSDPCHKERR(info);}
  if (dsdp->pdfeasible==DSDP_PDUNKNOWN) dsdp->pdfeasible=DSDP_PDFEASIBLE;
  info=DSDPEventLogEnd(dsdp->solvetime);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPCallMonitors"
/*!
\fn int DSDPCallMonitors(DSDP dsdp,DMonitor dmonitor[], int ndmonitors);
\param dsdp solver
\param dmonitor array of monitors
\param ndmonitors number of monitors.
\brief Call the monitor routines.
 */
int DSDPCallMonitors(DSDP dsdp,DMonitor dmonitor[], int ndmonitors){
  int i,info;
  DSDPFunctionBegin;
  for (i=0; i<ndmonitors;i++){
    info=(dmonitor[i].monitor)(dsdp,dmonitor[i].monitorctx);  DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}
/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "DSDPCheckConvergence"
/*!
\fn int DSDPCheckConvergence(DSDP dsdp,DSDPTerminationReason *reason);
\param dsdp solver
\param reason termination reason
\brief Check for convergence and monitor solution
 */
int DSDPCheckConvergence(DSDP dsdp,DSDPTerminationReason *reason){
  int info;
  DSDPTruth unbounded;

  DSDPFunctionBegin;
  info = DSDPGetConicDimension(dsdp,&dsdp->np);DSDPCHKERR(info);
  dsdp->rgap=(dsdp->ppobj-dsdp->ddobj)/(1.0+fabs(dsdp->ppobj)+fabs(dsdp->ddobj));
  dsdp->pstepold=dsdp->pstep;
  if (dsdp->reason==CONTINUE_ITERATING){
    if (dsdp->itnow>0){
      info=DSDPCheckForUnboundedObjective(dsdp,&unbounded);DSDPCHKERR(info);
      if (unbounded==DSDP_TRUE){
	dsdp->pdfeasible=DSDP_UNBOUNDED;
	info=DSDPSetConvergenceFlag(dsdp,DSDP_CONVERGED); DSDPCHKERR(info); 
      }
    }
    if (dsdp->reason==CONTINUE_ITERATING){
      if (dsdp->muold<dsdp->mutarget && dsdp->pstep==1 && dsdp->dstep==1 && dsdp->rgap<1e-5){
	info=DSDPSetConvergenceFlag(dsdp,DSDP_NUMERICAL_ERROR); DSDPCHKERR(info);
	DSDPLogInfo(0,2,"DSDP Finished: Numerical issues: Increase in Barrier function. \n");}
      if (dsdp->itnow >= dsdp->maxiter){
	info=DSDPSetConvergenceFlag(dsdp,DSDP_MAX_IT); DSDPCHKERR(info);} 
      if (dsdp->Mshift>dsdp->maxschurshift){
	info = DSDPSetConvergenceFlag(dsdp,DSDP_INDEFINITE_SCHUR_MATRIX); DSDPCHKERR(info);
      }
    } 
    info=DSDPCallMonitors(dsdp,dsdp->dmonitor,dsdp->nmonitors);DSDPCHKERR(info);
    info=DSDPMonitorCones(dsdp,0); DSDPCHKERR(info);
  }
  dsdp->muold=dsdp->mutarget;
  info = DSDPStopReason(dsdp,reason); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}



/* ---------------------------------------------------------- */ 
#undef __FUNCT__  
#define __FUNCT__ "DSDPTakeDown"
/*!
\fn int DSDPTakeDown(DSDP dsdp);
\param dsdp solver
\brief Destroy internal data structures.
 */
int DSDPTakeDown(DSDP dsdp){

  int i,info;

  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info = DSDPVecDestroy(&dsdp->rhs);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->rhs1);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->rhs2);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->rhstemp);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->y);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->ytemp);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->dy1);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->dy2);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->dy);DSDPCHKERR(info);
  for (i=0;i<MAX_XMAKERS;i++){
    info = DSDPVecDestroy(&dsdp->xmaker[i].y);DSDPCHKERR(info);
    info = DSDPVecDestroy(&dsdp->xmaker[i].dy);DSDPCHKERR(info);
    info = DSDPVecDestroy(&dsdp->xmaker[i].rhs);DSDPCHKERR(info);
  }
  info = DSDPVecDestroy(&dsdp->xmakerrhs);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->y0);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->b);DSDPCHKERR(info);

  info = DSDPCGDestroy(&dsdp->sles);DSDPCHKERR(info);
  info = DSDPDestroyCones(dsdp);DSDPCHKERR(info);
  info = DSDPSchurMatDestroy(&dsdp->M);DSDPCHKERR(info);
  info = DSDPGetConicDimension(dsdp,&dsdp->np);DSDPCHKERR(info);
  dsdp->setupcalled=DSDP_FALSE;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetDestroyRoutine(DSDP dsdp, int (*fd)(void*), void* ctx);

\brief Set a routine that will be called during DSDPDestroy().
\param dsdp the solver
\param fd function pointer
\param ctx pointer to structure.
\sa DSDPDestroy()
*/
int DSDPSetDestroyRoutine(DSDP dsdp, int (*fd)(void*), void* ctx){
  int nd=dsdp->ndroutines;
  if (nd<10){
    dsdp->droutine[nd].f=fd;
    dsdp->droutine[nd].ptr=ctx;
    dsdp->ndroutines++;
  } else {
    printf("TOO MANY Destroy routines\n");
    return 1;
  }
  return 0;
}


/*!
\fn int DSDPDestroy(DSDP dsdp)
\brief Free the internal data structures of the solver and
the cones associated with it.

\param dsdp the solver
\sa DSDPCreate()
\sa DSDPSolve()
\sa DSDPSetup()
\ingroup DSDPBasic
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPDestroy"
int DSDPDestroy(DSDP dsdp){
  int i,info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  for (i=0;i<dsdp->ndroutines;i++){
    info=(*dsdp->droutine[i].f)(dsdp->droutine[i].ptr);DSDPCHKERR(info);
  }
  info=DSDPTakeDown(dsdp);DSDPCHKERR(info);
  DSDPFREE(&dsdp,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
