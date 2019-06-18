#include "dsdp.h"
#include "dsdpsys.h"
#include "dsdp5.h"
/*!
\file dsdpx.c
\brief X variables, tolerances, errors, and feasibility.
*/

#undef __FUNCT__  
#define __FUNCT__ "DSDPInspectXY"
int DSDPInspectXY(DSDP dsdp, double xmakermu, DSDPVec xmakery, DSDPVec xmakerdy, DSDPVec AX, double *tracexs2, double *pobj2, double *rpinfeas2){
  int info;
  DSDPFunctionBegin;

  info=BoundYConeAddX(dsdp->ybcone,xmakermu,xmakery,xmakerdy,AX,tracexs2); DSDPCHKERR(info);
  info=DSDPVecGetC(AX,pobj2);DSDPCHKERR(info);

  info=DSDPVecSetC(AX,0);DSDPCHKERR(info);
  info=DSDPVecSetR(AX,0);DSDPCHKERR(info);
  info=DSDPVecNorm1(AX,rpinfeas2); DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}

/*!
\fn int DSDPComputeX(DSDP dsdp)
\brief Compute the X variables.

This routine explicitly computes X and determines the 
feasibility and unboundedness of the solution.
This routine should be called after DSDPSolve().

The cost of the dual-scaling direction is less than
the cost of other interior-point directions because the X
matrix does not have to be computed explicitly at each iteration.

\param dsdp the solver
\sa DSDPSolve()
\sa DSDPGetSolutionType()
\sa DSDPGetPObjective()
\ingroup DSDPBasic

These four routines can usually be called together.
\code
DSDP dsdp;
DSDPSolutionType type;

DSDPSetup(dsdp);
DSDPSolve(dsdp);
DSDPComputeX(dsdp);
DSDPGetSolutionType(dsdp,&type);
\endcode
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeX"
int DSDPComputeX(DSDP dsdp){
  int i,info;
  double pobj=0,ppobj2=0,ddobj,tracexs=0,tracexs2=0,rpinfeas=0,rpinfeas2=0,rpobjerr=0;
  double err1,cc,rrr,bigM,ymax,pfeastol=dsdp->pinfeastol;
  DSDPTerminationReason reason;
  DSDPVec AX=dsdp->ytemp;

  DSDPFunctionBegin;
  info=DSDPStopReason(dsdp,&reason);DSDPCHKERR(info);
  info=DSDPGetDDObjective(dsdp,&ddobj);DSDPCHKERR(info);
  info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
  info=DSDPGetR(dsdp,&rrr); DSDPCHKERR(info);
  info=DSDPGetPenalty(dsdp,&bigM);DSDPCHKERR(info);
  info=DSDPGetScale(dsdp,&cc);DSDPCHKERR(info);

  dsdp->pdfeasible=DSDP_PDFEASIBLE;
  for (i=0;i<MAX_XMAKERS;i++){
    if (i>0 && dsdp->xmaker[i].pstep<1) continue;
    info=DSDPComputeXVariables(dsdp,dsdp->xmaker[i].mu,dsdp->xmaker[i].y,dsdp->xmaker[i].dy,AX,&tracexs);DSDPCHKERR(info);
    info=DSDPVecGetC(AX,&pobj); DSDPCHKERR(info);    
    info=DSDPVecGetR(AX,&dsdp->tracex); DSDPCHKERR(info);
    info=DSDPVecSetC(AX,0);DSDPCHKERR(info);
    info=DSDPVecSetR(AX,0);DSDPCHKERR(info);
    info=DSDPVecNormInfinity(AX,&rpinfeas);DSDPCHKERR(info);
    rpinfeas=rpinfeas/(dsdp->tracex+1);

    DSDPLogInfo(0,2,"POBJ: %4.4e, DOBJ:  %4.4e\n",pobj,ddobj/cc);
    
    info=DSDPVecNorm2(AX,&err1);DSDPCHKERR(info);
    dsdp->tracexs=tracexs;
    dsdp->perror=err1;
    dsdp->pobj=cc*pobj;
    
    info=DSDPInspectXY(dsdp,dsdp->xmaker[i].mu,dsdp->xmaker[i].y,dsdp->xmaker[i].dy,AX,&tracexs2,&ppobj2,&rpinfeas2);DSDPCHKERR(info);
    rpinfeas2=rpinfeas2/(dsdp->tracex+1);
    /* rpinfeas is infeasibility of (P) while rpinfeas2 is infeasibility of (PP) */

    DSDPLogInfo(0,2,"X P Infeas: %4.2e , PObj: %4.8e\n",rpinfeas,pobj*(cc));
    DSDPLogInfo(0,2,"TOTAL  P Infeas: %4.2e PObj: %4.8e\n",rpinfeas2,ppobj2*(cc));
    rpobjerr= fabs(pobj-dsdp->ppobj)/(1+fabs(dsdp->ppobj));
    
    if (rpinfeas2 < pfeastol){ /* (PP) must be close to feasible */
      
      if (dsdp->rgap<0.1){
	if (rpinfeas>pfeastol/100 && fabs(rrr)>dsdp->dinfeastol){
	  dsdp->pdfeasible=DSDP_PDUNKNOWN;
	  DSDPLogInfo(0,2,"Warning: Try Increasing penalty parameter\n");
	} else if (rpinfeas>pfeastol && ddobj>0 && ppobj2<0 && fabs(rrr)<dsdp->dinfeastol){
	  dsdp->pdfeasible=DSDP_UNBOUNDED;
	  DSDPLogInfo(0,2,"Warning: D probably unbounded\n");
	  
	} else if (/* fabs(bigM)-dsdp->tracex < fabs(rrr) && rpinfeas<pfeastol */ fabs(rrr)>dsdp->dinfeastol){
	  dsdp->pdfeasible=DSDP_INFEASIBLE;
	  DSDPLogInfo(0,2,"Warning: D probably infeasible \n");
	}
      }
      i=i+10;
      break;

    } else { 
      /* Step direction was not accurate enough to compute X from Schur complement */
      DSDPLogInfo(0,2,"Try backup X\n");
      info=DSDPSetConvergenceFlag(dsdp,DSDP_NUMERICAL_ERROR); DSDPCHKERR(info); 
    }
    
  }

  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DSDPSaveBackupYForX"
int DSDPSaveBackupYForX(DSDP dsdp, int count,double mu, double pstep){
  int info;
  double pnorm;
  DSDPFunctionBegin;
  info=DSDPVecCopy(dsdp->y,dsdp->xmaker[count].y);DSDPCHKERR(info);
  info=DSDPComputeDY(dsdp,mu,dsdp->xmaker[count].dy,&pnorm); DSDPCHKERR(info);
  dsdp->xmaker[count].pstep=pstep;
  dsdp->xmaker[count].mu = mu;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSaveYForX(DSDP dsdp, double mu, double pstep);
\brief Save the current solution for later computation of X

\param dsdp the solver
\param mu barrier function
\param pstep primal step length, hopefully equals 1
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSaveYForX"
int DSDPSaveYForX(DSDP dsdp, double mu, double pstep){
  int info,isgood=0;
  double pnorm,newgap,ymax,dd=0;
  DSDPFunctionBegin;
  dsdp->pstepold=dsdp->pstep;
  info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
  if (pstep==0){
    info=DSDPVecCopy(dsdp->y,dsdp->xmaker[0].y);DSDPCHKERR(info);
    dsdp->xmaker[0].pstep=pstep;
  } else if (dsdp->Mshift*ymax>dsdp->pinfeastol*10){
    info=DSDPComputeDualityGap(dsdp,mu,&newgap);DSDPCHKERR(info);
    if (pstep==1 && newgap>0){
      dsdp->ppobj = dsdp->ddobj + newgap; dsdp->mu=(newgap)/(dsdp->np);
      dsdp->dualitygap=newgap;
    }
    info=DSDPVecZero(dsdp->rhstemp); DSDPCHKERR(info);
    info=BoundYConeAddX(dsdp->ybcone,dsdp->xmaker[0].mu,dsdp->xmaker[0].y,dsdp->xmaker[0].dy,dsdp->rhstemp,&dd); DSDPCHKERR(info);
    info=DSDPVecSetC(dsdp->rhstemp,0);
    info=DSDPVecSetR(dsdp->rhstemp,0);
    info=DSDPVecNormInfinity(dsdp->rhstemp,&dsdp->pinfeas); DSDPCHKERR(info);
    dsdp->pinfeas+=dsdp->Mshift*ymax;
    if (0==1){info=DSDPVecView(dsdp->rhstemp);}
    /* Not good enough to save */
  } else {
    info=DSDPVecCopy(dsdp->y,dsdp->xmaker[0].y);DSDPCHKERR(info);
    dsdp->xmaker[0].pstep=pstep;
    info=DSDPComputeRHS(dsdp,mu,dsdp->xmakerrhs); DSDPCHKERR(info);
    info=DSDPComputeDY(dsdp,mu,dsdp->xmaker[0].dy,&pnorm); DSDPCHKERR(info);
    dsdp->xmaker[0].mu = mu;
    info=DSDPComputeDualityGap(dsdp,mu,&newgap);DSDPCHKERR(info);
    if (pstep==1 && newgap>0){
      dsdp->ppobj = dsdp->ddobj + newgap; dsdp->mu=(newgap)/(dsdp->np);
      dsdp->dualitygap=newgap;

      info=DSDPVecZero(dsdp->rhstemp); DSDPCHKERR(info);
      info=BoundYConeAddX(dsdp->ybcone,dsdp->xmaker[0].mu,dsdp->xmaker[0].y,dsdp->xmaker[0].dy,dsdp->rhstemp,&dd); DSDPCHKERR(info);
      info=DSDPVecSetC(dsdp->rhstemp,0);
      info=DSDPVecSetR(dsdp->rhstemp,0);
      info=DSDPVecNormInfinity(dsdp->rhstemp,&dsdp->pinfeas); DSDPCHKERR(info);
      dsdp->pinfeas+=dsdp->Mshift*ymax;
      if (0==1){info=DSDPVecView(dsdp->rhstemp);}

    }
    isgood=1;
  }

  if (isgood==1){
    double penalty,r,rx;
    info=DSDPPassXVectors(dsdp,dsdp->xmaker[0].mu,dsdp->xmaker[0].y,dsdp->xmaker[0].dy); DSDPCHKERR(info);
    info=DSDPGetRR(dsdp,&r);DSDPCHKERR(info);
    if (r&& dsdp->rgap<0.1){  /* Estimate error in X */
      info=RConeGetRX(dsdp->rcone,&rx);DSDPCHKERR(info);
      info=DSDPGetPenalty(dsdp,&penalty);DSDPCHKERR(info);
      dsdp->pinfeas = dsdp->pinfeas *(1+fabs(penalty-rx));
    }
  }

  if (pstep==1.0 && dsdp->rgap>5.0e-1) {
    info=DSDPSaveBackupYForX(dsdp,MAX_XMAKERS-1,mu,pstep);DSDPCHKERR(info);
  }
  if (pstep==1.0 && dsdp->rgap>1.0e-3) {
    info=DSDPSaveBackupYForX(dsdp,2,mu,pstep);DSDPCHKERR(info);
  }
  if (pstep==1.0 && dsdp->rgap>1.0e-5) {
    info=DSDPSaveBackupYForX(dsdp,1,mu,pstep);DSDPCHKERR(info);
  }
  
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetPObjective(DSDP dsdp, double *pobj)
\brief Copy the objective value (P).
\param dsdp is the solver
\param *pobj will be the objective value in (P)
\sa DSDPGetDObjective()
\sa DSDPGetPPObjective()
\sa DSDPComputeX()
\note This value is correct only after calling DSDPComputeX()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPObjective"
int DSDPGetPObjective(DSDP dsdp,double *pobj){ 
  int info;
  double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  *pobj=(dsdp->pobj)/scale;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetSolutionType(DSDP dsdp, DSDPSolutionType *pdfeasible)
\brief Solutions can be bounded, infeasible, or unbounded.
\param dsdp is the solver
\param *pdfeasible will be set to the proper enumerated type.
\sa DSDPSetPenaltyParameter()
\sa DSDPSetYBounds()
\sa DSDPStopReason()
\ingroup DSDPBasic
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetSolutionType"
int DSDPGetSolutionType(DSDP dsdp,DSDPSolutionType *pdfeasible){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *pdfeasible=dsdp->pdfeasible;
  DSDPFunctionReturn(0);
}
/*!
\fn int DSDPGetTraceX(DSDP dsdp, double *tracex)
\brief Copy the trace of the variables X in (P).

For SDP blocks, this number corresponds to the trace
of the blocks, and for LP, it corresponds the the sum
of the variables x.  If this number is near the penalty
paramter, the problem (P) may be unbounded or the
penalty parameter may have to be increased.

\param dsdp is the solver
\param *tracex will be set the trace of the variables in (P)
\sa DSDPSetPenaltyParameter()
\sa DSDPComputeX()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetTraceX"
int DSDPGetTraceX(DSDP dsdp, double *tracex){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *tracex=dsdp->tracex;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetFinalErrors(DSDP dsdp, double err[6])
\brief Copy six different error measurements into an array.
\param dsdp is the solver
\param err will be set to the six error measurements
\sa DSDPGetDualityGap()
\sa DSDPGetR()
\sa DSDPGetPInfeasibility()
\ingroup DSDPSolution
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetFinalErrors"
int DSDPGetFinalErrors(DSDP dsdp, double err[6]){
  int info;
  double scale,rr,bnorm,dobj=0,pobj=0;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  info=DSDPVecGetR(dsdp->y,&rr); DSDPCHKERR(info);
  info=DSDPGetPObjective(dsdp,&pobj);DSDPCHKERR(info);
  info=DSDPGetDObjective(dsdp,&dobj);DSDPCHKERR(info);
  err[0]=dsdp->perror;
  err[1]=0;
  err[2]=fabs(rr)/scale;
  err[3]=0;
  err[4]=pobj - dobj;
  err[5]=dsdp->tracexs/scale;

  err[2] /= (1.0+dsdp->cnorm);
  info=DSDPVecCopy(dsdp->b,dsdp->ytemp);DSDPCHKERR(info);
  info=DSDPVecSetC(dsdp->ytemp,0);DSDPCHKERR(info);
  info=DSDPVecSetR(dsdp->ytemp,0);DSDPCHKERR(info);
  info=DSDPVecNormInfinity(dsdp->ytemp,&bnorm);
  err[0]=dsdp->perror/(1.0+bnorm);

  err[4]=(err[4])/(1.0+fabs(pobj)+fabs(dobj));
  err[5]=(err[5])/(1.0+fabs(pobj)+fabs(dobj));
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetPInfeasibility(DSDP dsdp, double *pperror)
\brief Copy the infeasibility in (P).

This infeasibility is the maximum difference between the
values xl and xu correponding to the bounds on the variables y.
Due to roundoff error, this number is usually much less than
the true infeasiblity in (P).  The true infeasibility is
not available at each iteration due to its high computational cost.

\param dsdp is the solver
\param *pperror will be set to the infeasibility in (P)
\sa DSDPSetYBounds()
\sa DSDPSetPTolerance()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPInfeasibility"
int DSDPGetPInfeasibility(DSDP dsdp, double *pperror){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (pperror) *pperror=dsdp->pinfeas;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetPTolerance(DSDP dsdp, double inftol)
\brief Classify (P) as feasible only if the infeasibility is less than this tolerance

\param dsdp is the solver
\param inftol will be set to the infeasibility in (P) (>0)
\sa DSDPSetYBounds()
\sa DSDPGetPInfeasibility()
\sa DSDPGetPTolerance()
\sa DSDPGetRTolerance()
\sa DSDPSolutionType
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetPTolerance"
int DSDPSetPTolerance(DSDP dsdp,double inftol){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (inftol > 0) dsdp->pinfeastol  = inftol;
  DSDPLogInfo(0,2,"Set P Infeasibility Tolerance: %4.4e\n",inftol);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetPTolerance(DSDP dsdp, double *inftol)
\brief Copy the feasibility tolerance.

\param dsdp is the solver
\param *inftol will be set to the infeasibility in (P)
\sa DSDPSetYBounds()
\sa DSDPGetPInfeasibility()
\sa DSDPSetPTolerance()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetPTolerance"
int DSDPGetPTolerance(DSDP dsdp,double *inftol){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (inftol) *inftol=dsdp->pinfeastol;
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPSetRTolerance(DSDP dsdp, double inftol)
\brief Classify (D) as feasible only if the variable r is
less than this tolerance.

\param dsdp is the solver
\param inftol is the tolerance (>0)
\sa DSDPSetPenaltyParameter()
\sa DSDPGetR()
\sa DSDPGetRTolerance()
\sa DSDPSolutionType
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetRTolerance"
int DSDPSetRTolerance(DSDP dsdp,double inftol){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (inftol > 0) dsdp->dinfeastol  = inftol;
  DSDPLogInfo(0,2,"Set D Infeasibility Tolerance: %4.4e\n",inftol);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetRTolerance(DSDP dsdp, double *inftol)
\brief Copy the maximum infeasibility allowed (D).

DSDP will classify the solution to (D) as feasible only
if the variable r in (DD) is less than this tolerance.
Small values for r are enforced through a penalty parameter.

\param dsdp is the solver
\param *inftol will be set to the tolerance for r in (DD)
\sa DSDPSetPenaltyParameter()
\sa DSDPGetR()
\sa DSDPSetRTolerance()
\ingroup DSDPConverge
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetRTolerance"
int DSDPGetRTolerance(DSDP dsdp,double *inftol){ 
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *inftol=dsdp->dinfeastol;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetYMakeX(DSDP dsdp, double y[], int m)
\brief Copies the variables y used to construct X into an array.
\param dsdp is the solver
\param y is an array
\param m is the length of the array and the dimension of y
\sa DSDPGetY()
\sa DSDPComputeX()
\sa DSDPGetDYMakeX()
\sa DSDPGetMuMakeX()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetYMakeX"
int DSDPGetYMakeX(DSDP dsdp,double y[], int m){ 
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
\fn int DSDPGetDYMakeX(DSDP dsdp, double dy[], int m)
\brief Copies the variables dy used to construct X into an array.
\param dsdp is the solver
\param dy is an array
\param m is the length of the array and the dimension of dy
\sa DSDPComputeX()
\sa DSDPGetYMakeX()
\sa DSDPGetMuMakeX()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetDYMakeX"
int DSDPGetDYMakeX(DSDP dsdp,double dy[], int m){ 
  int i,info;
  double scale,*yy;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  if (dsdp->m < m-1) DSDPFunctionReturn(1);
  if (dsdp->m > m) DSDPFunctionReturn(1);
  info=DSDPVecCopy(dsdp->xmaker[0].dy,dsdp->ytemp); DSDPCHKERR(info);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  info=DSDPVecGetArray(dsdp->ytemp,&yy);DSDPCHKERR(info);
  for (i=0;i<m;i++) dy[i]=yy[i+1]/scale;
  info=DSDPVecRestoreArray(dsdp->ytemp,&yy);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPGetMuMakeX(DSDP dsdp, double *mu)
\brief Copies the value of mu used to construct X.
\param dsdp is the solver
\param mu is the barrier parameter
\sa DSDPComputeX()
\sa DSDPGetYMakeX()
\sa DSDPGetDYMakeX()
\sa DSDPGetBarrierParameter()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPGetMuMakeX"
int DSDPGetMuMakeX(DSDP dsdp,double *mu){ 
  int info;
  double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPGetScale(dsdp,&scale);DSDPCHKERR(info);
  *mu=dsdp->xmaker[0].mu/scale;
  DSDPFunctionReturn(0);
}

