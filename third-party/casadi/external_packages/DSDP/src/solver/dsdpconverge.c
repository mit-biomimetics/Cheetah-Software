#include "dsdp5.h"
#include "dsdpconverge.h"
/*!
  \file dsdpconverge.c
  \brief Monitor convergence.
*/


extern int DSDPGetConvergenceMonitor(DSDP, ConvergenceMonitor**);

/*!
\fn int DSDPDefaultConvergence(DSDP dsdp, void *ctx)
\brief Check for Convergence.

\param dsdp is the solver
\param ctx is a pointer to a structure containing convergence parameters
\sa DSDPSetGapTolerance()
\sa DSDPSetStepTolerance()
\ingroup DSDPConverge

\note DSDP calls this routine before each iteration.

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPCheckConvergence"
int DSDPDefaultConvergence(DSDP dsdp,void *ctx){

  ConvergenceMonitor *conv=(ConvergenceMonitor*)ctx;
  int info,i,iter;
  double mu,mu2;
  double rgap,rgap2,rgaptol=conv->rgaptol;
  double infeastol=0;
  double pnorm,dstep,pstep,steptol=conv->steptol,pnormtol=conv->pnormtol;
  double ppobj,ddobj, gap, dualbound=conv->dualbound;
  double res,np;
  DSDPTerminationReason reason;
  
  DSDPFunctionBegin;
  info = DSDPGetStepLengths(dsdp,&pstep,&dstep); DSDPCHKERR(info);
  info = DSDPGetPnorm(dsdp,&pnorm); DSDPCHKERR(info);
  info = DSDPGetIts(dsdp,&iter); DSDPCHKERR(info);
  info = DSDPGetDDObjective(dsdp,&ddobj);  DSDPCHKERR(info);
  info = DSDPGetPPObjective(dsdp,&ppobj); DSDPCHKERR(info);
  info = DSDPGetR(dsdp,&res); DSDPCHKERR(info);
  info = DSDPGetBarrierParameter(dsdp,&mu); DSDPCHKERR(info);
  info = DSDPGetDimension(dsdp,&np); DSDPCHKERR(info);
  info = DSDPStopReason(dsdp,&reason); DSDPCHKERR(info);
  info = DSDPGetRTolerance(dsdp,&infeastol); DSDPCHKERR(info);
  info = DSDPGetDualityGap(dsdp,&gap); DSDPCHKERR(info);
  rgap=(gap)/(1.0+fabs(ddobj)/2+fabs(ppobj)/2);
  rgap2=(mu*np)/(1.0+fabs(ddobj)/2+fabs(ppobj)/2);
  if (iter==0){
    conv->history = DSDPHistory;
    for (i=0; i<DSDPHistory; i++){
      conv->alpha[i] = 0.0;
      conv->gaphist[i] = 0.0;
      conv->infhist[i] = 0.0;
    }
  }
  if (iter<conv->history && iter>0){
    conv->gaphist[iter-1]=(ppobj-ddobj);
    conv->infhist[iter-1]=res;
  }

  if ( 0==1 ){
    
  } else if ( ddobj!=ddobj || pnorm < 0){
    reason = DSDP_NUMERICAL_ERROR;
    DSDPLogInfo(0,2,"Stop due to Numerical Error\n");
  } else if ( rgap <=rgaptol/1.01 && res<=infeastol ){
    if (pnorm>pnormtol){
      mu2=gap/np;
      info = DSDPSetBarrierParameter(dsdp,mu2); DSDPCHKERR(info);
    } else {
      reason = DSDP_CONVERGED;
      DSDPLogInfo(0,2,"DSDP Converged:  Relative Duality Gap %4.2e < %4.2e, Primal Feasible, Dual Infeasiblity %4.2e < %4.2e \n",rgap,rgaptol,res,infeastol);
    }
  } else if ( rgap2 <=rgaptol/100 && rgap<0.01){
    reason = DSDP_CONVERGED;
    DSDPLogInfo(0,2,"DSDP Converged:  Relative Duality Gap %4.2e < %4.2e. Check Feasiblity \n",rgap,rgaptol);
  } else if ( ddobj > dualbound && res<=infeastol){
    reason = DSDP_UPPERBOUND;
    DSDPLogInfo(0,2,"DSDP Converged: Dual Objective: %4.2e > upper bound %4.2e\n",pnorm,dualbound);
  } else if ( iter > 5 && dstep<steptol && dstep*pnorm< steptol && rgap <= 1.0e-3 ) {
    reason = DSDP_SMALL_STEPS;
    DSDPLogInfo(0,2,"DSDP Terminated:  Small relative gap and small steps detected (3)\n");
  }

  info=DSDPSetConvergenceFlag(dsdp,reason); DSDPCHKERR(info); 

  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetGapTolerance(DSDP dsdp, double gaptol)
\brief Terminate the solver when the relative duality gap is less than this tolerance.

If pp is the objective value of (PP) and dd is the objective value of (DD),
the relative duality gap is defined to be (pp-dd)/(1 + fabs(pp)/2 + fabs(dd)/2);

\param dsdp is the solver
\param gaptol is the tolerance
\sa DSDPGetDualityGap()
\sa DSDPGetGapTolerance()
\sa DSDP_CONVERGED
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetGapTolerance"
int DSDPSetGapTolerance(DSDP dsdp,double gaptol){ 
  int info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  if (gaptol > 0) conv->rgaptol  = gaptol;
  DSDPLogInfo(0,2,"Set Relative Gap Tolerance: %4.4e\n",gaptol);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetGapTolerance(DSDP dsdp, double *gaptol)
\brief Get the termination tolerance.

\param dsdp is the solver
\param *gaptol will be set to the termination tolerance
\sa DSDPGetDualityGap()
\sa DSDPSetGapTolerance()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetGapTolerance"
int DSDPGetGapTolerance(DSDP dsdp,double *gaptol){ 
  int info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  DSDPFunctionBegin;
  *gaptol=conv->rgaptol;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetPNormTolerance(DSDP dsdp, double ptol)
\brief Terminate the solver when the relative duality gap is suffiently
small and the PNorm is less than this quantity. Smaller values imply
the final solution will be nearer to the central path.

\param dsdp is the solver
\param ptol is the tolerance, (>0, default is very big)
\sa DSDPGetPnorm()
\sa DSDPGetPNormTolerance()
\sa DSDPSetGapTolerance()
\sa DSDP_CONVERGED
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetPNormTolerance"
int DSDPSetPNormTolerance(DSDP dsdp,double ptol){ 
  int info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  if (ptol > 0) conv->pnormtol  = ptol;
  DSDPLogInfo(0,2,"Set Relative PNorm Tolerance: %4.4e\n",ptol);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetPNormTolerance(DSDP dsdp, double *ptol)
\brief Get the termination tolerance.

\param dsdp is the solver
\param *ptol will be set to the termination tolerance
\sa DSDPGetPnorm()
\sa DSDPSetPNormTolerance()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPNormTolerance"
int DSDPGetPNormTolerance(DSDP dsdp,double *ptol){ 
  int info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  DSDPFunctionBegin;
  *ptol=conv->pnormtol;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetDualBound(DSDP dsdp, double dbound)
\brief Terminate the solver if the objective value in (DD) is greater than this tolerance.

This parameter is helpful in branch and bound applictions.

\param dsdp is the solver
\param dbound is the bound
\sa DSDPGetDualBound()
\sa DSDPGetDDObjective()
\sa DSDP_UPPERBOUND
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDualBound"
int DSDPSetDualBound(DSDP dsdp,double dbound){ 
  int info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  conv->dualbound=dbound;
  DSDPLogInfo(0,2,"Set DualBound of %4.4e \n",dbound);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetDualBound(DSDP dsdp, double *dbound)
\brief Get the termination parameter.

\param dsdp is the solver
\param *dbound is a bound on (DD)
\sa DSDPSetDualBound()
\sa DSDPGetDDObjective()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDualBound"
int DSDPGetDualBound(DSDP dsdp,double *dbound){ 
  int info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  *dbound=conv->dualbound;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetStepTolerance(DSDP dsdp, double steptol)
\brief Terminate the solver if the step length in (DD) is below this tolerance

This heuristic will only be applied when the objective values
in (PP) and (DD) match to three significant digits.

\param dsdp is the solver
\param steptol is the tolerance
\sa DSDPGetStepTolerance()
\sa DSDPGetStepLengths()
\sa DSDP_SMALL_STEPS
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetStepTolerance"
int DSDPSetStepTolerance(DSDP dsdp,double steptol){ 
  int info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  if (steptol > 0) conv->steptol   = steptol;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetStepTolerance(DSDP dsdp, double *steptol)
\brief Get the current tolerance.

\param dsdp is the solver
\param *steptol will be set to the current tolerance
\sa DSDPSetStepTolerance()
\sa DSDPGetStepLengths()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetStepTolerance"
int DSDPGetStepTolerance(DSDP dsdp,double *steptol){ 
  int info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  *steptol=conv->steptol;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetRHistory(DSDP dsdp, double hist[], int length)
\brief Copy a history of the infeasibility in (D) into an array

Elememt i of the array will be set the the infeasibility in (D) at
iteration i (unless i exceeds the length of the array)

\param dsdp is the solver
\param hist is an array
\param length is the length of the array
\sa DSDPGetR()
\sa DSDPGetGapHistory()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetRHistory"
int DSDPGetRHistory(DSDP dsdp, double hist[], int length){ 
  int i,info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  for (i=0;i<length;i++) hist[i]=0.0;
  for (i=0;i<DSDPMin(length,DSDPHistory);i++) hist[i]=conv->infhist[i];
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetGapHistory(DSDP dsdp, double hist[], int length)
\brief Copy a history of the duality gap into an array

\param dsdp is the solver
\param hist is an array
\param length is the length of the array
\sa DSDPGetDualityGap()
\sa DSDPGetRHistory()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetGapHistory"
int DSDPGetGapHistory(DSDP dsdp, double hist[], int length){ 
  int i,info;
  ConvergenceMonitor *conv;
  DSDPFunctionBegin;
  info=DSDPGetConvergenceMonitor(dsdp,&conv); DSDPCHKERR(info);
  for (i=0;i<length;i++) hist[i]=0.0;
  for (i=0;i<DSDPMin(length,DSDPHistory);i++) hist[i]=conv->gaphist[i];
  DSDPFunctionReturn(0);
}

