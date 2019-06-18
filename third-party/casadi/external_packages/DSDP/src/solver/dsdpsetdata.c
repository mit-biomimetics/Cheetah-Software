/*!
  \file dsdpsetdata.c
  \brief Set parameters for the solver and retrieve statistics about the current solution.
*/

#include "dsdp.h"
#include "dsdp5.h"
#include "dsdpsys.h"

/*!
\fn int DSDPSetDualObjective(DSDP dsdp, int i, double bi)
\brief Set the objective vector b in (D).
\param dsdp is the solver
\param i is the variable number from 1 through m
\param bi is the objective value associated with variable i
\sa DSDPSetY0()
\sa DSDPGetDObjective()
\ingroup DSDPBasic

The dual objective function is \f$ \mbox{maximize} \ \ {\displaystyle \sum_{i=1}^m b_i \ y_i } \f$.

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDualObjective"
int DSDPSetDualObjective(DSDP dsdp,int i, double bi){ 
  int info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (i>dsdp->m || i<=0){
    DSDPSETERR2(1,"Invalid variable number: Is 1 <= %d <= %d?\n",i,dsdp->m);}
  info=DSDPVecSetElement(dsdp->b,i,bi);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPCopyB(DSDP dsdp, double bb[], int m)
\brief Copies the variables b from solver into an array.
\param dsdp is the solver
\param bb is an array
\param m is the length of the array and the dimension of y
\sa DSDPSetDualObjective()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPCopyB"
int DSDPCopyB(DSDP dsdp,double bb[], int m){ 
  int i,info;
  double *b;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (dsdp->m < m) DSDPFunctionReturn(1);
  info=DSDPVecGetArray(dsdp->b,&b);DSDPCHKERR(info);
  for (i=0;i<m;i++) bb[i]=b[i+1];
  info=DSDPVecRestoreArray(dsdp->b,&b);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPSetY0(DSDP dsdp,int i, double yi0)
\brief Set the initial values of variables y in (D).

To improve performance consider setting the initial 
values of the variables y in (D).

\param dsdp is the solver
\param i is the variable number from 1 through m
\param yi0 is the initial value af that variable 
\sa DSDPGetY()
\sa DSDPSetR0()
\sa DSDPSetPotentialParameter()
\sa DSDPReuseMatrix()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetY0"
int DSDPSetY0(DSDP dsdp,int i, double yi0){ 
  int info;double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (i>dsdp->m || i<=0){
    DSDPSETERR2(1,"Invalid variable number: Is 1<= %d <= %d\n",i,dsdp->m);}
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  info=DSDPVecSetElement(dsdp->y,i,scale*yi0);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetY(DSDP dsdp, double y[], int m)
\brief Copies the variables y into an array.
\param dsdp is the solver
\param y is an array
\param m is the length of the array and the dimension of y
\sa DSDPSetY0()
\sa DSDPComputeX()
\ingroup DSDPBasic
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetY"
int DSDPGetY(DSDP dsdp,double y[], int m){ 
  int i,info;
  double scale,*yy;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (dsdp->m < m-1) DSDPFunctionReturn(1);
  if (dsdp->m > m) DSDPFunctionReturn(1);
  info=DSDPVecCopy(dsdp->xmaker[0].y,dsdp->ytemp); DSDPCHKERR(info);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  info=DSDPVecGetArray(dsdp->ytemp,&yy);DSDPCHKERR(info);
  for (i=0;i<m;i++) y[i]=yy[i+1]/scale;
  info=DSDPVecRestoreArray(dsdp->ytemp,&yy);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}



/*!
\fn int DSDPGetScale(DSDP dsdp, double *scale)
\brief Copy the internal scaling factor from the solver

\param dsdp is the solver
\param *scale will be set to the scaling factor used in the solver
\sa DSDPSetScale()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetScale"
int DSDPGetScale(DSDP dsdp,double *scale){ 
  int info;double sscale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPVecGetC(dsdp->y,&sscale);DSDPCHKERR(info);
  *scale=fabs(sscale);
  if (sscale==0) *scale=1.0;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetScale(DSDP dsdp, double scale)
\brief Set the internal scaling factor.

DSDP scales the data C and solves the scaled problem.
DSDP determines a default scaling from a combination
of the norms on the data.

\param dsdp is the solver
\param scale is the scaling factor used in the solver (>0)
\sa DSDPGetScale()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetScale"
int DSDPSetScale(DSDP dsdp,double scale){ 
  int info;double sscale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  scale=fabs(scale);
  if (scale==0) scale=1.0;
  info=DSDPGetScale(dsdp,&sscale);DSDPCHKERR(info);
  sscale=scale/sscale;
  info=DSDPVecScale(sscale,dsdp->y);
  dsdp->mutarget*=sscale;
  dsdp->pobj*=sscale;
  dsdp->dobj*=sscale;
  dsdp->ppobj*=sscale;
  dsdp->ddobj*=sscale;
  dsdp->mu*=sscale;
  DSDPLogInfo(0,2,"Set DSDP C Scaling: %4.4e\n",scale);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPAddObjectiveConstant(DSDP dsdp, double c)
\brief Add a constant to the objective.
\param dsdp is the solver
\param c is the constant.
\sa DSDPGetDObjective()
\note This parameter does not affect the performance of the
solver.  It can, however, make the standout output more
consistent with the underlying application.
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPAddObjectiveConstant"
int DSDPAddObjectiveConstant(DSDP dsdp,double c){ 
  int info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPVecSetC(dsdp->b,-c);DSDPCHKERR(info);
  DSDPLogInfo(0,2,"Add Objective Constant: %4.4e\n",c);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetMaxIts(DSDP dsdp, int maxits)
\brief Terminate the solver after this number of iterations.
\param dsdp is the solver
\param maxits is the maximum number of DSDP iterations (>0)
\sa DSDPGetMaxIts()
\sa DSDPGetIts()
\sa DSDPSetGapTolerance()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetMaxIts"
int DSDPSetMaxIts(DSDP dsdp,int its){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (its >= 0) dsdp->maxiter = its;
  DSDPLogInfo(0,2,"Set Maximum Iterates: %4d\n",its);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetMaxIts(DSDP dsdp, int *maxits)
\brief Copy the maximum number of iterations from the solver
\param dsdp is the solver
\param *maxits will be the maximum number of iterations in DSDP
\sa DSDPSetMaxIts()
\sa DSDPGetIts()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetMaxIts"
int DSDPGetMaxIts(DSDP dsdp,int *its){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *its=dsdp->maxiter;
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPSetMaxTrustRadius(DSDP dsdp, double rad)
\brief Set a maximum trust radius on the step direction.

\param dsdp is the solver
\param rad is radius of the trust region.(default: 1e30)
\sa DSDPGetMaxTrustRadius()
\ingroup DSDPSolver

\note By default this tolerance is very large and does not 
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetMaxTrustRadius"
int DSDPSetMaxTrustRadius(DSDP dsdp,double rad){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (rad > 0) dsdp->maxtrustradius   = rad;
  DSDPLogInfo(0,2,"Set Maximum Trust Radius: %4.4e\n",rad);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetMaxTrustRadius(DSDP dsdp, double *rad)
\brief Copy the current radius of the trust region.

\param dsdp is the solver
\param *rad will be set to radius of the trust region
\sa DSDPSetMaxTrustRadius()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetMaxTrustRadius"
int DSDPGetMaxTrustRadius(DSDP dsdp,double *rad){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *rad=dsdp->maxtrustradius;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetZBar(DSDP dsdp, double ppobj)
\brief Set an upper bound on the objective value at the solution

\param dsdp is the solver
\param ppobj is the initial objective v (default: 1e30)
\sa DSDPGetPPObjective()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetZBar"
int DSDPSetZBar(DSDP dsdp,double ppobj){ 
  int info;
  double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  dsdp->ppobj=ppobj*scale;
  DSDPLogInfo(0,2,"Set Primal Objective and Upper bound on solution: %4.4e. \n",ppobj);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetR0(DSDP dsdp, double r0)
\brief Set an initial value for the variable r in (DD)

A negative value asks DSDP to choose this parameter.
The default heuristic generally uses very large values.
Smaller values may significantly improve performance.

\param dsdp is the solver
\param r0 is the initial objective v (default: -1)
\sa DSDPSetPenaltyParameter()
\sa DSDPGetR()
\sa DSDPSetY0()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetR0"
int DSDPSetR0(DSDP dsdp,double res){ 
  int info;
  double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  info=DSDPSetRR(dsdp,scale*res); DSDPCHKERR(info);
  if (res>=0)dsdp->goty0=DSDP_TRUE;
  DSDPLogInfo(0,2,"Set Dual Initial Infeasibility to %4.4e times Identity Matrix. \n",res);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetBarrierParameter(DSDP dsdp, double mu)
\brief Set the current barrier parameter.

The barrier parameter is defined as the difference between the
objective functions (PP) and (DD) divided by the potential
parameter rho.

\param dsdp is the solver
\param mu is the initial objective v
\sa DSDPGetBarrierParameter()
\sa DSDPSetZBar()
\sa DSDPSetPotentialParameter()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetBarrierParameter"
int DSDPSetBarrierParameter(DSDP dsdp,double mu){ 
  int info;double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  dsdp->mutarget = mu*scale;
  DSDPLogInfo(0,2,"Set InitialBarrierParameter: %4.4e \n",mu);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetBarrierParameter(DSDP dsdp, double *mu)
\brief Copy the current barrier parameter.

\param dsdp is the solver
\param *mu barrier parameter
\sa DSDPSetBarrierParameter()
\sa DSDPGetPPObjective()
\sa DSDPGetDDObjective()
\sa DSDPGetPotentialParameter()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetBarrierParameter"
int DSDPGetBarrierParameter(DSDP dsdp, double *mu){ 
  int info;double scale;
  DSDPFunctionBegin;
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  *mu=dsdp->mutarget/scale;
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPUsePenalty(DSDP dsdp, int yesorno)
\brief Use penalty parameter to enforce feasibility.
\param dsdp is the solver
\param yesorno is the decision
\sa DSDPSetPenaltyParameter()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPUsePenalty"
int DSDPUsePenalty(DSDP dsdp,int yesorno){ 
  DSDPPenalty UsePenalty;
  int info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (yesorno>0){
    UsePenalty=DSDPAlways;
  } else if (yesorno<0){
    UsePenalty=DSDPNever;
  } else {
    UsePenalty=DSDPInfeasible;
  }
  dsdp->UsePenalty=UsePenalty;
  info=RConeSetType(dsdp->rcone,UsePenalty);DSDPCHKERR(info);
  DSDPLogInfo(0,2,"Set UsePenalty: %d \n",yesorno);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetPenaltyParameter(DSDP dsdp, double Gamma)
\brief Set the penalty parameter Gamma.

DSDP uses a penalty parameter to enforce feasibility in (D).
The default value is 1e8, but other values may improve performance.
This value must exceed the trace of the solution X.

\param dsdp is the solver
\param Gamma is the penalty parameter
\sa DSDPGetPenaltyParameter()
\sa DSDPGetR()
\sa DSDPGetTraceX()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetPenaltyParameter"
int DSDPSetPenaltyParameter(DSDP dsdp,double Gamma){ 
  int info;
  double scale,ppenalty;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  ppenalty=fabs(Gamma*scale);
  info=DSDPVecSetR(dsdp->b,ppenalty);DSDPCHKERR(info);
  DSDPLogInfo(0,2,"Set Penalty Parameter: %4.4e\n",Gamma);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetPenaltyParameter(DSDP dsdp, double *Gamma)
\brief Copy the penalty parameter Gamma.

\param dsdp is the solver
\param *Gamma wil be set to the penalty parameter
\sa DSDPSetPenaltyParameter()
\sa DSDPGetR()
\sa DSDPGetDDObjective()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPenaltyParameter"
int DSDPGetPenaltyParameter(DSDP dsdp,double *Gamma){ 
  int info;
  double ppenalty;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPVecGetR(dsdp->b,&ppenalty);DSDPCHKERR(info);
  *Gamma=fabs(ppenalty);
  DSDPFunctionReturn(0);
}

/* Not current; not documented
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPenalty"
int DSDPGetPenalty(DSDP dsdp,double *penalty){ 
  int info;double ppenalty;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPVecGetR(dsdp->b,&ppenalty);DSDPCHKERR(info);
  *penalty=fabs(ppenalty);
  DSDPFunctionReturn(0);
}



/*!
\fn int DSDPGetPPObjective(DSDP dsdp, double *ppobj)
\brief Copy the objective value (PP).
\param dsdp is the solver
\param *ppobj will be the objective value in (PP)
\sa DSDPGetDDObjective()
\sa DSDPGetPObjective()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPPObjective"
int DSDPGetPPObjective(DSDP dsdp,double *ppobj){ 
  int info;
  double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  *ppobj=(dsdp->ppobj)/scale;
  if (dsdp->cnorm==0) *ppobj=0;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetDObjective(DSDP dsdp, double *dobj)
\brief Copy the objective value (D).
\param dsdp is the solver
\param *dobj will be the objective value in (D)
\sa DSDPGetPObjective()
\sa DSDPGetDDObjective()
\sa DSDPGetY()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDObjective"
int DSDPGetDObjective(DSDP dsdp,double *dobj){ 
  int info; double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  *dobj = (dsdp->dobj)/scale;
  if (dsdp->cnorm==0) *dobj=-fabs(*dobj);
  DSDPFunctionReturn(0);
}
/*!
\fn int DSDPGetDDObjective(DSDP dsdp, double *ddobj)
\brief Copy the objective value (DD).
\param dsdp is the solver
\param *ddobj will be the objective value in (DD)
\sa DSDPGetPPObjective()
\sa DSDPGetDObjective()
\sa DSDPGetY()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDDObjective"
int DSDPGetDDObjective(DSDP dsdp,double *ddobj){ 
  int info; double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  *ddobj = (dsdp->ddobj)/scale;
  if (dsdp->cnorm==0) *ddobj=-fabs(*ddobj);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetDualityGap(DSDP dsdp, double *dgap)
\brief Copy the difference between the objective values.
\param dsdp is the solver
\param *dgap will be set to the difference between the objective values in (PP) and (DD)
\sa DSDPGetDDObjective()
\sa DSDPGetPPObjective()
\sa DSDPGetDimension()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDualityGap"
int DSDPGetDualityGap(DSDP dsdp,double *dgap){ 
  int info; double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  *dgap = (dsdp->dualitygap)/scale;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetIts(DSDP dsdp, int *its)
\brief Copy the current iteration number.
\param dsdp is the solver
\param *its will be set to the current iteration number.
\sa DSDPSetMaxIts()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetIts"
int DSDPGetIts(DSDP dsdp,int *its){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *its=dsdp->itnow;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPStopReason(DSDP dsdp, DSDPTerminationReason *reason)
\brief Copy the reason why the solver terminated.
\param dsdp is the solver
\param *reason will be set to the proper enumerated type.
\sa DSDPSetMaxIts()
\sa DSDPSetGapTolerance()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPStopReason"
int DSDPStopReason(DSDP dsdp,DSDPTerminationReason *reason){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *reason=dsdp->reason;
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPGetR(DSDP dsdp, double *res)
\brief Copy the infeasibility in (D), or the variable r in (DD).
\param dsdp is the solver
\param *res will be set to the value of r.
\sa DSDPSetMaxIts()
\sa DSDPSetGapTolerance()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetR"
int DSDPGetR(DSDP dsdp, double *res){ 
  int info;double rr,scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetRR(dsdp,&rr);DSDPCHKERR(info);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  *res=rr/scale;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetDataNorms(DSDP dsdp, double dnorm[3])
\brief Copy the norms of the data C, A, and b into an array.
\param dsdp is the solver
\param dnorm will be set the norms the data C, A, and b.
\sa DSDPSetDualObjective()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDataNorms"
int DSDPGetDataNorms(DSDP dsdp, double dnorm[3]){
  int info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (dsdp->setupcalled==DSDP_FALSE){
    info=DSDPComputeDataNorms(dsdp);DSDPCHKERR(info);
  }
  dnorm[0]=dsdp->cnorm;
  dnorm[1]=dsdp->anorm;
  dnorm[2]=dsdp->bnorm;
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPGetMaxYElement(DSDP dsdp, double *ymax)
\brief Copy the the infinity norm of the variables y.
\param dsdp is the solver
\param *ymax will be set to the magnitude of the largest variable y.
\sa DSDPSetYBounds()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetMaxYElement"
int DSDPGetMaxYElement(DSDP dsdp,double* ymax){
  int info;
  DSDPFunctionBegin;
  info=DSDPGetYMaxNorm(dsdp,ymax);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDimension"
/*! \fn int DSDPGetDimension(DSDP dsdp, double *n);
\brief Copy the dimension of the cones, or the number of constraints in (D).
\param dsdp the solver
\param *n will be set to the dimension (a whole number)
\ingroup DSDPSolution
\sa DSDPGetNumberOfVariables()
*/
int DSDPGetDimension(DSDP dsdp, double *n){
  int info;
  DSDPFunctionBegin;
  info=DSDPGetConicDimension(dsdp,n);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetYMaxNorm(DSDP dsdp, double *ynorm)
\brief Copy the the infinity norm of the variables y.
\param dsdp is the solver
\param *ynorm will be set to the magnitude of the largest variable y.
\sa DSDPSetYBounds()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetYMaxNorm"
int DSDPGetYMaxNorm(DSDP dsdp, double *ynorm){ 
  int info;
  double cc,rr,yy;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPVecGetC(dsdp->y,&cc);DSDPCHKERR(info);
  info=DSDPVecGetR(dsdp->y,&rr);DSDPCHKERR(info);
  info=DSDPVecSetC(dsdp->y,0);DSDPCHKERR(info);
  info=DSDPVecSetR(dsdp->y,0);DSDPCHKERR(info);
  info=DSDPVecNormInfinity(dsdp->y,&yy);DSDPCHKERR(info);
  info=DSDPVecSetC(dsdp->y,cc);DSDPCHKERR(info);
  info=DSDPVecSetR(dsdp->y,rr);DSDPCHKERR(info);
  if (cc) yy/=fabs(cc);
  if (ynorm) *ynorm=yy;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetNumberOfVariables(DSDP dsdp, int *m)
\brief Copy the number of variables y.

\param dsdp the solver
\param *m will be set the number of variables y
\sa DSDPCreate()
\sa DSDPGetDimension()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetNumberOfVariables"
int DSDPGetNumberOfVariables(DSDP dsdp, int *m){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *m=dsdp->m;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetPnorm(DSDP dsdp, double *pnorm)
\brief Copy the proximity of the solution to the central path.
\param dsdp is the solver
\param *pnorm will be set a norm of the gradient of the barrier function
\sa DSDPSetGapTolerance()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPnorm"
int DSDPGetPnorm(DSDP dsdp, double *pnorm){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *pnorm=dsdp->pnorm;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetStepLengths(DSDP dsdp, double *pstep, double *dstep)
\brief Copy the step sizes in the current iteration.
\param dsdp is the solver
\param *pstep will be set to the step size in (PP)
\param *dstep will be set to the step size in (DD)
\sa DSDPSetStepTolerance()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetStepLengths"
int DSDPGetStepLengths(DSDP dsdp, double *pstep, double *dstep){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *dstep=dsdp->dstep;
  *pstep=dsdp->pstep;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetPotentialParameter(DSDP dsdp, double rho)
\brief Set the potential parameter.

The parameter rho in the solver will be set to this multiple
of the total dimension of the cones.  The default value is
3, but values of 4 or more may significantly improve performance.

\param dsdp is the solver
\param rho the potential parameter.
\sa DSDPGetPotentialParameter()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetPotentialParameter"
int DSDPSetPotentialParameter(DSDP dsdp, double rho){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (rho>1) dsdp->rhon=rho;
  DSDPLogInfo(0,2,"Set Potential Parameter %4.4f\n",rho);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetPotentialParameter(DSDP dsdp, double *rho)
\brief Copy the potential parameter.

\param dsdp is the solver
\param *rho will be set to the potential parameter
\sa DSDPSetPotentialParameter()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPotentialParameter"
int DSDPGetPotentialParameter(DSDP dsdp, double *rho){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *rho=dsdp->rhon;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetPotential(DSDP dsdp, double *potential)
\brief Copy the potential of the current solution.

\param dsdp is the solver
\param *potential will be set to the value of the potential function
\sa DSDPSetPotentialParameter()
\sa DSDPGetDDObjective()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPotential"
int DSDPGetPotential(DSDP dsdp, double *potential){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *potential=dsdp->potential;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPUseDynamicRho(DSDP dsdp, int yesorno)
\brief Use a dynamic strategy to choose parameter rho.
\param dsdp is the solver
\param yesorno is nonzero to use the dynamic strategy (default 1)
\sa DSDPSetPotentialParameter()
\sa DSDPGetPotential()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPUseDynamicRho"
int DSDPUseDynamicRho(DSDP dsdp, int yesorno){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (yesorno)  dsdp->usefixedrho=DSDP_FALSE;
  else  dsdp->usefixedrho=DSDP_TRUE;
  DSDPLogInfo(0,2,"Set UseDynamicRho: %d \n",yesorno);
  DSDPFunctionReturn(0);
}

/* Not Current or documented
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPBoundDualVariables"
/* !
\fn int DSDPBoundDualVariables(DSDP dsdp, double lbound, double ubound)
\brief Bounds on the variables y.

\param dsdp is the solver
\param lbound will be the lower bound of the variables y
\param ubound will be the upper bound of the variables y
\sa DSDPSetYBounds()
\ingroup DSDPSolver
*/
int DSDPBoundDualVariables(DSDP dsdp,double lbound, double ubound){
  int info;
  double bbound;
  DSDPFunctionBegin;
  bbound=DSDPMax(fabs(lbound),fabs(ubound));
  DSDPLogInfo(0,2,"Bound Variables between %4.4e and %4.4e \n",-bbound,bbound);
  info = BoundYConeSetBounds(dsdp->ybcone,-bbound,bbound);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

/*!
\fn int DSDPGetYBounds(DSDP dsdp, double *lbound, double *ubound)
\brief Copy the bounds on the variables y.

\param dsdp is the solver
\param *lbound will be set to the lower bound of the variables y
\param *ubound will be set to the upper bound of the variables y
\sa DSDPSetYBounds()
\ingroup DSDPSolver
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPGetYBounds"
int DSDPGetYBounds(DSDP dsdp,double *lbound, double *ubound){
  int info;
  DSDPFunctionBegin;
  info=BoundYConeGetBounds(dsdp->ybcone,lbound,ubound);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

/*!
\fn int DSDPSetYBounds(DSDP dsdp, double lbound, double ubound)
\brief Bound the variables y.

\param dsdp is the solver
\param lbound is the lower bound for each variable y.
\param ubound is the upper bound for each variable y.
\sa DSDPSetYBounds()
\ingroup DSDPSolver
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSetYBounds"
int DSDPSetYBounds(DSDP dsdp,double lbound, double ubound){
  int info;
  DSDPFunctionBegin;
  info=BoundYConeSetBounds(dsdp->ybcone,lbound,ubound);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}



/*!
\fn int DSDPReuseMatrix(DSDP dsdp, int rm)
\brief Reuse the Hessian of the barrier function multiple times at each DSDP iteration.

\param dsdp is the solver
\param rm is the maximum number of times the matrix will be used in each DSDP iteration
\sa DSDPGetReuseMatrix()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPReuseMatrix"
int DSDPReuseMatrix(DSDP dsdp, int rm){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  dsdp->reuseM=rm;
  DSDPLogInfo(0,2,"Reuse the Schur Matrix: %d times\n",rm);
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPGetReuseMatrix(DSDP dsdp, int *rm)
\brief Copy this parameter.

\param dsdp is the solver
\param *rm will be set to the maximum number of times the matrix will be reused
\sa DSDPReuseMatrix()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetReuseMatrix"
int DSDPGetReuseMatrix(DSDP dsdp, int *rm){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *rm=dsdp->reuseM;
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPSetMonitor(DSDP dsdp, int (*monitor)(DSDP,void*), void* monitorctx )
\brief Monitor each iteration of the solver.

\param dsdp is the solver
\param monitor is a function that will be called at each iteration
\param monitorctx is a pointer that will be passed to the function
\sa DSDPSetStandardMonitor()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetMonitor"
int DSDPSetMonitor(DSDP dsdp, int (*monitor)(DSDP,void*), void* monitorctx){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (dsdp->nmonitors<MAX_DSDP_MONITORS){
    DSDPLogInfo(0,2,"Set Monitor\n");
    dsdp->dmonitor[dsdp->nmonitors].monitor=monitor;
    dsdp->dmonitor[dsdp->nmonitors].monitorctx=monitorctx;
    dsdp->nmonitors++;
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetConvergenceFlag(DSDP dsdp, DSDPTerminationReason reason ){
\brief Monitor each iteration of the solver.

\param dsdp is the solver
\param reason is the termination reason
\sa DSDPStopReason()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetConvergenceFlag"
int DSDPSetConvergenceFlag(DSDP dsdp, DSDPTerminationReason reason ){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  dsdp->reason=reason;
  if (reason==DSDP_INFEASIBLE_START){
    DSDPLogInfo(0,2,"Initial Point Infeasible, Check variable bounds? \n",0);
  }
  DSDPFunctionReturn(0);
}

