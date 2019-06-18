#include "dsdp.h"
#include "dsdpsys.h"
/*!
\file dualimpl.c
\brief Dual-scaling operations needed in the solver.
*/

/*!
\fn int DSDPComputeObjective(DSDP dsdp, DSDPVec Y, double *ddobj);
\brief Compute the objective function (DD).

\param dsdp is the solver
\param Y Current variables
\param ddobj objective value

\sa DSDPComputeNewY()
\sa DSDPGetDDObjective()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeObjective"
int DSDPComputeObjective(DSDP dsdp, DSDPVec Y, double *ddobj){
  int info;
  DSDPFunctionBegin;
  info = DSDPVecDot(Y,dsdp->b,ddobj);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPComputeDY(DSDP dsdp, double mu, DSDPVec DY, double *pnorm);
\brief Compute the step direction.

\param dsdp the solver
\param mu barrier parameter
\param DY Step direction
\param pnorm distance to the target

Assuming the affine direction and centering direction have alread been
computed, combine them with the appropriate barrier parameter.

\sa DSDPComputeRHS()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeDY"
int DSDPComputeDY(DSDP dsdp, double mu, DSDPVec DY, double *pnorm){
  int info;
  double ppnorm,ddy1=fabs(1.0/mu*dsdp->schurmu),ddy2=-1.0;
  DSDPFunctionBegin;
  info=DSDPComputeRHS(dsdp,mu,dsdp->rhs); DSDPCHKERR(info);
  info=DSDPVecWAXPBY(DY,ddy1,dsdp->dy1,ddy2,dsdp->dy2);DSDPCHKERR(info);
  info=DSDPComputePNorm(dsdp,mu,DY,&ppnorm);DSDPCHKERR(info);
  if (ppnorm<0){ /* If pnorm < 0 there are SMW numerical issues */
    DSDPLogInfo(0,2,"Problem with PNORM: %4.4e < 0 \n",ppnorm);
    /*    ppnorm=1.0; */
  }
  *pnorm=ppnorm;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputePDY(DSDP dsdp, double mu, DSDPVec DY, double *pnorm);
\brief Compute the step direction.

\param dsdp the solver
\param mu barrier parameter
\param DY Step direction
\param pnorm distance to the target

Assuming the affine direction and centering direction have alread been
computed, combine them with the appropriate barrier parameter.

\sa DSDPComputeRHS()
\sa DSDPComputeDY()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputePDY"
int DSDPComputePDY(DSDP dsdp, double mu, DSDPVec DY, double *pnorm){
  int info;
  double ppnorm,ddy1=-fabs(1.0/mu*dsdp->schurmu),ddy2=1.0;
  DSDPFunctionBegin;
  info=DSDPComputeRHS(dsdp,mu,dsdp->rhs); DSDPCHKERR(info);
  info=DSDPVecWAXPBY(DY,ddy1,dsdp->dy1,ddy2,dsdp->dy2);DSDPCHKERR(info);
  info=DSDPComputePNorm(dsdp,mu,DY,&ppnorm);DSDPCHKERR(info);
  if (ppnorm<0){ /* If pnorm < 0 there are SMW numerical issues */
    DSDPLogInfo(0,2,"Problem with PNORM: %4.4e < 0 \n",ppnorm);
    /*    ppnorm=1.0; */
  }
  *pnorm=ppnorm;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputePDY1(DSDP dsdp, double mur, DSDPVec DY1);
\brief Compute an affine step direction dy1.

\param dsdp the solver
\param mur reciprocal of barrier parameter
\param DY1 Step direction

Assuming the affine direction has alread been computed, scale it.
\sa DSDPComputeDY()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputePDY1"
int DSDPComputePDY1(DSDP dsdp, double mur, DSDPVec DY1){
  int info;
  double ddy1=-fabs(mur*dsdp->schurmu);
  DSDPFunctionBegin;
  info=DSDPVecScaleCopy(dsdp->dy1,ddy1,DY1); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputeNewY(DSDP dsdp, double beta, DSDPVec Y);
\brief Update the Y variables.

\param dsdp the solver
\param beta step length
\param Y the new solution

Add a multiple of the step direction to the current solution.
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeNewY"
int DSDPComputeNewY(DSDP dsdp, double beta, DSDPVec Y){
  int info;
  double rtemp;
  DSDPFunctionBegin;
  info=DSDPVecWAXPY(Y,beta,dsdp->dy,dsdp->y);DSDPCHKERR(info);
  info=DSDPVecGetR(Y,&rtemp);DSDPCHKERR(info);
  rtemp=DSDPMin(0,rtemp);
  info=DSDPSchurMatSetR(dsdp->M,rtemp);DSDPCHKERR(info);
  info=DSDPVecSetR(Y,rtemp);DSDPCHKERR(info);
  info=DSDPApplyFixedVariables(dsdp->M,Y);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputePY(DSDP dsdp, double beta, DSDPVec PY);
\brief Compute PY = Y - beta DY for use in computing X.

\param dsdp the solver
\param beta step length
\param PY the new value

\sa DSDPComputeNewY()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputePY"
int DSDPComputePY(DSDP dsdp, double beta, DSDPVec PY){
  int info;
  DSDPFunctionBegin;
  info=DSDPVecWAXPY(PY,beta,dsdp->dy,dsdp->y);DSDPCHKERR(info);
  info=DSDPApplyFixedVariables(dsdp->M,PY);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputeRHS(DSDP dsdp, double mu, DSDPVec RHS);
\brief Compute the right-hand side of the linear system that determines
the step direction.

\param dsdp the solver
\param mu barrier parameter
\param RHS right-hand side direction

Assuming that the gradient of the objective and the gradient of the
barrier have already been
computed, combine them with the appropriate barrier parameter.

This vector is basically \f$ b - mu * A(S^{-1}) \f$

\sa DSDPComputeDY()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeRHS"
int DSDPComputeRHS(DSDP dsdp, double mu, DSDPVec RHS){
  int info;
  double ddrhs1=1.0/mu*dsdp->schurmu,ddrhs2=-( mu/fabs(mu) );
  DSDPFunctionBegin;
  info=DSDPVecWAXPBY(RHS,ddrhs1,dsdp->rhs1,ddrhs2,dsdp->rhs2);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPComputePNorm"
/*!
\fn int DSDPComputePNorm(DSDP dsdp, double mu, DSDPVec DY, double *pnorm);

\brief Compute proximity to a point on the central path.

\param dsdp the solver
\param mu barrier parameter
\param DY Newton step direction
\param pnorm the norm 

\sa DSDPComputeDY()

*/
int DSDPComputePNorm(DSDP dsdp, double mu, DSDPVec DY, double *pnorm){
  int info;
  double ppnorm=0;
  DSDPFunctionBegin;
  info=DSDPComputeRHS(dsdp,mu,dsdp->rhs); DSDPCHKERR(info);
  info = DSDPVecDot(dsdp->rhs,DY,&ppnorm);DSDPCHKERR(info);
  ppnorm/=dsdp->schurmu;
  if (ppnorm >= 0){  /* Theoretically True */
    *pnorm=sqrt(ppnorm);
  } else {
    DSDPLogInfo(0,2,"Problem with PNORM: %4.4e is not positive.\n",ppnorm);
    *pnorm=ppnorm;
  }
  if (*pnorm!=*pnorm){DSDPSETERR1(1,"Problem with PNORM: %4.4e is not positive.\n",ppnorm);}
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputeDualityGap(DSDP dsdp, double mu, double *gap);
\brief Compute the current duality gap.

\param dsdp the solver
\param mu barrier parameter
\param gap the duality gap

\sa DSDPGetDualityGap()
\sa DSDPGetPPObjective()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeDualityGap"
int DSDPComputeDualityGap(DSDP dsdp, double mu, double *gap){
  int info;
  double newgap=0,pnorm;
  double smu=1.0/dsdp->schurmu;
  DSDPFunctionBegin;
  info=DSDPComputeDY(dsdp,mu,dsdp->dy,&pnorm); DSDPCHKERR(info);
  info=DSDPVecDot(dsdp->dy,dsdp->rhs2,&newgap);DSDPCHKERR(info);
  newgap = (newgap*smu+dsdp->np)*mu;
  if (newgap<=0){
    DSDPLogInfo(0,2,"GAP :%4.4e<0: Problem\n",newgap);
  } else {
    DSDPLogInfo(0,2,"Duality Gap: %12.8e, Update primal objective: %12.8e\n",newgap,dsdp->ddobj+newgap);
  }
  newgap=DSDPMax(0,newgap);
  *gap=newgap;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputePotential(DSDP dsdp, DSDPVec y, double logdet, double *potential);
\brief Compute the potential of the given point.

\param dsdp the solver
\param y variables
\param logdet logarithmic barrier function of the given point
\param potential return the potential of the given point.

\sa DSDPSetPotentialParameter()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputePotential"
int DSDPComputePotential(DSDP dsdp, DSDPVec y, double  logdet, double *potential){
  int info;
  double dpotential,gap,ddobj;
  DSDPFunctionBegin;
  info=DSDPComputeObjective(dsdp,y,&ddobj);DSDPCHKERR(info);
  gap=dsdp->ppobj-ddobj;
  if (gap>0) dpotential=dsdp->rho*log(gap)-logdet;
  else {dpotential=dsdp->potential+1;}
  *potential=dpotential;
  DSDPLogInfo(0,9,"Gap: %4.4e, Log Determinant: %4.4e, Log Gap: %4.4e\n",gap,logdet,log(gap));
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputePotential2(DSDP dsdp, DSDPVec y, double mu, double logdet, double *potential);
\brief Compute the objective function plus the barrier function.

\param dsdp the solver
\param y variables
\param mu barrier function
\param logdet logarithmic barrier function of the given point
\param potential return the potential of the given point.\
\sa DSDPGetBarrierParameter()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputePotential2"
int DSDPComputePotential2(DSDP dsdp, DSDPVec y, double mu, double logdet, double *potential){
  int info;
  double ddobj;
  DSDPFunctionBegin;
  info=DSDPComputeObjective(dsdp,y,&ddobj);DSDPCHKERR(info);
  *potential=-(ddobj + mu*logdet)*dsdp->schurmu;
  *potential=-(ddobj/mu + logdet)*dsdp->schurmu;
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPSetY(DSDP dsdp, double beta, double logdet, DSDPVec ynew);
\brief Update the solver with these y variables

\param dsdp the solver
\param beta most recent step length
\param logdet logarithmic barrier function of the given point
\param ynew current solution.
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSetY"
int DSDPSetY(DSDP dsdp, double beta, double logdet, DSDPVec ynew){
  int info;
  double r1,r2,rr,pp;
  DSDPFunctionBegin;
  info=DSDPVecGetR(dsdp->y,&r1);DSDPCHKERR(info);
  info=DSDPVecGetR(ynew,&r2);DSDPCHKERR(info);
  if (r2==0&&r1!=0){dsdp->rflag=1;} else {dsdp->rflag=0;};
  info=DSDPVecCopy(ynew,dsdp->y);DSDPCHKERR(info);
  info = DSDPComputeObjective(dsdp,dsdp->y,&dsdp->ddobj);DSDPCHKERR(info);
  /* Correct ppobj if ppobj < ddobj , which can happen when dual infeasibility is present */
  if (dsdp->ppobj<=dsdp->ddobj){
    dsdp->ppobj=dsdp->ddobj+2*dsdp->mu * dsdp->np;
    DSDPLogInfo(0,2,"Primal Objective Not Right.  Assigned: %8.8e\n",dsdp->ppobj);
  }
  info=DSDPVecGetR(ynew,&rr);DSDPCHKERR(info);
  info=DSDPVecGetR(dsdp->b,&pp);DSDPCHKERR(info);
  dsdp->dobj=dsdp->ddobj-rr*pp;
  DSDPLogInfo(0,2,"Duality Gap: %4.4e, Potential: %4.4e \n",dsdp->dualitygap,dsdp->potential);
  dsdp->dualitygap=dsdp->ppobj-dsdp->ddobj;
  dsdp->mu=(dsdp->dualitygap)/(dsdp->np);
  dsdp->dstep=beta;
  dsdp->logdet=logdet;
  info=DSDPComputePotential(dsdp,dsdp->y,dsdp->logdet,&dsdp->potential);DSDPCHKERR(info);
  DSDPLogInfo(0,2,"Duality Gap: %4.4e, Potential: %4.4e \n",dsdp->dualitygap,dsdp->potential);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPSetRR"
/*!
\fn int DSDPSetRR(DSDP dsdp,double res);
\param dsdp solver
\param res variable r
\brief Set variable r.
 */
int DSDPSetRR(DSDP dsdp,double res){ 
  int info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPVecSetR(dsdp->y,-res);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetRR"
/*!
\fn int DSDPGetRR(DSDP dsdp,double *res);
\param dsdp solver
\param *res set variable r
\brief Get variable r.
 */
int DSDPGetRR(DSDP dsdp,double *res){ 
  int info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info=DSDPVecGetR(dsdp->y,res);DSDPCHKERR(info);
  *res=-*res;
  if (*res==0) *res=0;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPObjectiveGH"
/*!
\fn int DSDPObjectiveGH( DSDP dsdp , DSDPSchurMat M, DSDPVec vrhs1);
\param dsdp solver
\param M Schur matrix.
\param vrhs1 gradient vector
\brief Compute gradient of dual objective.
*/
int DSDPObjectiveGH( DSDP dsdp , DSDPSchurMat M, DSDPVec vrhs1){
  int i,info,m;
  double rtemp,dd;

  DSDPFunctionBegin;
  info=DSDPVecGetSize(vrhs1,&m); DSDPCHKERR(info);
  for (i=0;i<m;i++){
    info=DSDPSchurMatVariableCompute(M,i,&dd); DSDPCHKERR(info);
    if (dd){
      info=DSDPVecGetElement(dsdp->b,i,&rtemp);DSDPCHKERR(info);
      info=DSDPVecAddElement(vrhs1,i,rtemp);DSDPCHKERR(info);
    }
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPCheckForUnboundedObjective"
int DSDPCheckForUnboundedObjective(DSDP dsdp, DSDPTruth *unbounded){
  int info;
  double dtemp;
  DSDPTruth psdefinite;
  DSDPFunctionBegin;
  *unbounded=DSDP_FALSE;
  info = DSDPVecDot(dsdp->b,dsdp->dy2,&dtemp);DSDPCHKERR(info);
  if ( dtemp < 0 /* && dsdp->r==0 && dsdp->ddobj > 0 */) {  
    info = DSDPVecScaleCopy(dsdp->dy2,-1.0,dsdp->ytemp); DSDPCHKERR(info);
    info = DSDPComputeSS(dsdp,dsdp->ytemp,PRIMAL_FACTOR,&psdefinite);DSDPCHKERR(info);
    if (psdefinite == DSDP_TRUE){
      psdefinite=DSDP_FALSE;
      while(psdefinite==DSDP_FALSE){ /* Dual point should be a certificate of dual unboundedness, and be dual feasible */
	info=DSDPComputeSS(dsdp,dsdp->ytemp,PRIMAL_FACTOR,&psdefinite);DSDPCHKERR(info);
	if (psdefinite == DSDP_TRUE) break;
	info=DSDPVecScale(2.0,dsdp->ytemp); DSDPCHKERR(info);
      }
      info = DSDPVecCopy(dsdp->ytemp,dsdp->y); DSDPCHKERR(info);
      info = DSDPSaveYForX(dsdp,1.0e-12,1.0);DSDPCHKERR(info);
      info = DSDPComputeObjective(dsdp,dsdp->y,&dsdp->ddobj);DSDPCHKERR(info);
      info = DSDPVecNormalize(dsdp->y); DSDPCHKERR(info);
      *unbounded=DSDP_TRUE;
    }
  }
  DSDPFunctionReturn(0);
}

