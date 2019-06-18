#include "dsdp.h"
#include "dsdpsys.h"
/*!
  \file dualalg.c
  \brief Implements the dual-scaling algorithm.
*/

int DSDPChooseBarrierParameter(DSDP,double,double*,double*);
int DSDPYStepLineSearch(DSDP,double,double,DSDPVec);
int DSDPYStepLineSearch2(DSDP,double,double,DSDPVec);
int DSDPResetY0(DSDP);

#undef __FUNCT__
#define __FUNCT__ "DSDPYStepLineSearch"
/*!
\fn int DSDPYStepLineSearch(DSDP dsdp, double mutarget, double dstep0, DSDPVec dy);
\param dsdp the solver
\param mutarget barrier parameter
\param dstep0 initial step length
\param dy step direction
\brief Used for Newton step, the merit function of this line search is
the dual potential function.
*/
int DSDPYStepLineSearch(DSDP dsdp, double mutarget, double dstep0, DSDPVec dy){
  /* The merit function is the dual potential function */
  int info,attempt,maxattempts=30;
  double dstep,newpotential,logdet;
  double better=0.05, steptol=1e-8,maxmaxstep=0;
  DSDPTruth psdefinite;
  DSDPFunctionBegin;
  info=DSDPComputeMaxStepLength(dsdp,dy,DUAL_FACTOR,&maxmaxstep);DSDPCHKERR(info);
  info=DSDPComputePotential(dsdp,dsdp->y,dsdp->logdet,&dsdp->potential);DSDPCHKERR(info);
  if (dsdp->pnorm<0.5) better=0.0;
  dstep=DSDPMin(dstep0,0.95*maxmaxstep);
  if (dstep * dsdp->pnorm > dsdp->maxtrustradius) dstep=dsdp->maxtrustradius/dsdp->pnorm;
  DSDPLogInfo(0,8,"Full Dual StepLength %4.4e, %4.4e\n",maxmaxstep,dstep);
  psdefinite=DSDP_FALSE;
  for (psdefinite=DSDP_FALSE,attempt=0; attempt<maxattempts && psdefinite==DSDP_FALSE; attempt++){
    info=DSDPComputeNewY(dsdp,dstep,dsdp->ytemp);DSDPCHKERR(info);
    info=DSDPComputeSS(dsdp,dsdp->ytemp,DUAL_FACTOR,&psdefinite);DSDPCHKERR(info);
    if (psdefinite==DSDP_TRUE){      
      info=DSDPComputeLogSDeterminant(dsdp,&logdet);DSDPCHKERR(info);
      info=DSDPComputePotential(dsdp,dsdp->ytemp,logdet,&newpotential);DSDPCHKERR(info);
      if (newpotential>dsdp->potential-better && dstep > 0.001/dsdp->pnorm ){
	DSDPLogInfo(0,2,"Not sufficient reduction. Reduce stepsize.  Trust Radius: %4.4e\n",dstep*dsdp->pnorm);
	psdefinite=DSDP_FALSE; dstep=0.3*dstep;
      } 
    } else {
      dstep=dstep/3.0;
      DSDPLogInfo(0,2,"Dual Matrix not Positive Definite: Reduce step %4.4e",dstep);
    }
    if (dstep*dsdp->pnorm < steptol && dstep < steptol) break;
  } /* Hopefully, the initial step size works and only go through loop once */
  if (psdefinite==DSDP_TRUE){
    info=DSDPSetY(dsdp,dstep,logdet,dsdp->ytemp);DSDPCHKERR(info);
  } else {
    info=DSDPSetY(dsdp,0,dsdp->logdet,dsdp->y);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPYStepLineSearch2"
/*!
\fn int DSDPYStepLineSearch2(DSDP dsdp, double mutarget, double dstep0, DSDPVec dy);
\param dsdp the solver
\param mutarget barrier parameter
\param dstep0 initial step length
\param dy step direction
\brief Used for centering steps, the merit function of this line search is
the objective function plus the barrier term.
*/
int DSDPYStepLineSearch2(DSDP dsdp, double mutarget, double dstep0, DSDPVec dy){
  /* The merit function is the objective (DD) plus the barrier function */
  /* This line search is used in the corrector steps */
  int info, attempt, maxattempts=10;
  double dstep,newpotential,bdotdy,oldpotential,logdet;
  double maxmaxstep=0,steptol=1e-6;
  double a,b;
  DSDPTruth psdefinite;
  DSDPFunctionBegin;
  info=DSDPComputeMaxStepLength(dsdp,dy,DUAL_FACTOR,&maxmaxstep);DSDPCHKERR(info);
  info=DSDPComputePotential2(dsdp,dsdp->y,mutarget, dsdp->logdet,&oldpotential);DSDPCHKERR(info);
  info=DSDPVecDot(dsdp->rhs,dy,&bdotdy);DSDPCHKERR(info);
  dstep=DSDPMin(dstep0,0.95*maxmaxstep);
  if (dstep * dsdp->pnorm > dsdp->maxtrustradius) dstep=dsdp->maxtrustradius/dsdp->pnorm;
  DSDPLogInfo(0,8,"Full Dual StepLength %4.4e, %4.4e\n",maxmaxstep,dstep);
  for (psdefinite=DSDP_FALSE,attempt=0; attempt<maxattempts && psdefinite==DSDP_FALSE; attempt++){
    if (dstep < steptol) break;
    info=DSDPComputeNewY(dsdp,dstep,dsdp->ytemp);DSDPCHKERR(info);
    info=DSDPComputeSS(dsdp,dsdp->ytemp,DUAL_FACTOR,&psdefinite);DSDPCHKERR(info);
    if (psdefinite==DSDP_TRUE){
      info=DSDPComputeLogSDeterminant(dsdp,&logdet);DSDPCHKERR(info);
      info=DSDPComputePotential2(dsdp,dsdp->ytemp,mutarget,logdet,&newpotential);DSDPCHKERR(info);
      b=bdotdy; a=2*(newpotential-oldpotential+bdotdy*dstep)/(dstep*dstep);
      if (newpotential>oldpotential-0.1*dstep*bdotdy ){
	DSDPLogInfo(0,2,"Not sufficient reduction. Reduce stepsize.  Step:: %4.4e\n",dstep);
	psdefinite=DSDP_FALSE;
	if (b/a<dstep && b/a>0){ dstep=b/a;} else { dstep=dstep/2; } 
      } 
    } else {
      dstep=dstep/2.0;
      DSDPLogInfo(0,2,"Dual Matrix not Positive Definite: Reduce step %4.4e",dstep);
    }
  } /* Hopefully, the initial step size works and only go through loop once */
  if (psdefinite==DSDP_TRUE && dstep>=steptol){
    info=DSDPSetY(dsdp,dstep,logdet,dsdp->ytemp);DSDPCHKERR(info);
  } else {
    info=DSDPSetY(dsdp,0,dsdp->logdet,dsdp->y);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPSolveDynmaicRho"
/*!
\fn int DSDPSolveDynamicRho(DSDP dsdp);
\brief Apply dual-scaling algorithm.
\param dsdp the solver
*/
int DSDPSolveDynamicRho(DSDP dsdp){

  int info,attempt,maxattempts;
  double dd1,dd2,mutarget,ppnorm;
  DSDPTerminationReason reason;
  DSDPTruth cg1;
  DSDPTruth psdefinite;

  DSDPFunctionBegin;

  info=DSDPVecCopy(dsdp->y,dsdp->y0);DSDPCHKERR(info);
  for (dsdp->itnow=0; dsdp->itnow <= dsdp->maxiter+1 ; dsdp->itnow++){

    /* Check Convergence, and print information if desired */
    info=DSDPCheckConvergence(dsdp,&reason);DSDPCHKERR(info);
    if (reason != CONTINUE_ITERATING){break;}
    if (dsdp->mu0>0){dsdp->mutarget=DSDPMin(dsdp->mutarget,dsdp->mu0);}
 
    /* Compute the Gram matrix M and rhs */
    info=DSDPComputeDualStepDirections(dsdp); DSDPCHKERR(info);
    if (dsdp->reason==DSDP_INDEFINITE_SCHUR_MATRIX){continue;}

    info=DSDPComputePDY(dsdp,dsdp->mutarget,dsdp->dy,&dsdp->pnorm); DSDPCHKERR(info);

    DSDPEventLogBegin(dsdp->ptime);
    info=DSDPComputePY(dsdp,1.0,dsdp->ytemp);DSDPCHKERR(info);
    info=DSDPComputeSS(dsdp,dsdp->ytemp,PRIMAL_FACTOR,&psdefinite);DSDPCHKERR(info);
    if (psdefinite==DSDP_TRUE){
      dsdp->pstep=1.0;
      info=DSDPSaveYForX(dsdp,dsdp->mutarget,dsdp->pstep);DSDPCHKERR(info);
    } else {
      dsdp->pstep=0.0;
    }

    if (dsdp->usefixedrho==DSDP_TRUE){
      dsdp->rho=dsdp->rhon*dsdp->np;
      mutarget=(dsdp->ppobj-dsdp->ddobj)/(dsdp->rho);
      dsdp->pstep=0.5;
    } else {
      info = DSDPChooseBarrierParameter(dsdp,dsdp->mutarget,&dsdp->pstep,&mutarget);DSDPCHKERR(info);
      dsdp->rho=(dsdp->ppobj-dsdp->ddobj)/(mutarget);
    }
    DSDPEventLogEnd(dsdp->ptime);
    
    DSDPLogInfo(0,6,"Current mu=%4.8e, Target X with mu=%4.8e, Rho: %8.4e, Primal Step Length: %4.8f, pnorm: %4.8e\n",dsdp->mu,mutarget,dsdp->rho/dsdp->np,dsdp->pstep, dsdp->pnorm);


    /* Take Dual Step */
    /* Compute distance from chosen point on central path Pnorm */
    DSDPEventLogBegin(dsdp->dtime);
    info=DSDPComputeDY(dsdp,mutarget,dsdp->dy,&dsdp->pnorm); DSDPCHKERR(info);
    if (dsdp->pnorm<0.1){ mutarget/=10;  info=DSDPComputeDY(dsdp,mutarget,dsdp->dy,&dsdp->pnorm); DSDPCHKERR(info);}

    info=DSDPYStepLineSearch(dsdp, mutarget, 1.0, dsdp->dy);DSDPCHKERR(info);
    DSDPEventLogEnd(dsdp->dtime);

    maxattempts=dsdp->reuseM;
    if (dsdp->dstep<1 && dsdp->rgap<1e-5) maxattempts=0;
    if (dsdp->dstep<1e-13) maxattempts=0;
    if (dsdp->rgap<1e-6) maxattempts=0;
    if (maxattempts>dsdp->reuseM) maxattempts=dsdp->reuseM;
    for (attempt=0;attempt<maxattempts;attempt++){
      double cgtol=1e-6;
      if (attempt>0 && dsdp->pnorm < 0.1) break;
      if (attempt > 0 && dsdp->dstep<1e-4) break;
      if (dsdp->rflag) break;
      DSDPEventLogBegin(dsdp->ctime);
      DSDPLogInfo(0,2,"Reuse Matrix %d: Ddobj: %12.8e, Pnorm: %4.2f, Step: %4.2f\n",attempt,dsdp->ddobj,dsdp->pnorm,dsdp->dstep);
      info=DSDPInvertS(dsdp);DSDPCHKERR(info);
      info=DSDPComputeG(dsdp,dsdp->rhstemp,dsdp->rhs1,dsdp->rhs2);DSDPCHKERR(info);
      if (dsdp->slestype==2 || dsdp->slestype==3){
	if (dsdp->rflag){info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs1,dsdp->dy1,cgtol,&cg1);DSDPCHKERR(info);}
	info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs2,dsdp->dy2,cgtol,&cg1);DSDPCHKERR(info);
      }
      info=DSDPVecDot(dsdp->b,dsdp->dy1,&dd1);DSDPCHKERR(info);
      info=DSDPVecDot(dsdp->b,dsdp->dy2,&dd2);DSDPCHKERR(info);
      if (dd1>0 && dd2>0){
	mutarget=DSDPMin(mutarget,dd1/dd2*dsdp->schurmu);
      }
      mutarget=mutarget*(dsdp->np/(dsdp->np+sqrt(dsdp->np)));	  
      info=DSDPComputeDY(dsdp,mutarget,dsdp->dy, &ppnorm);DSDPCHKERR(info);
      if (ppnorm<=0){ DSDPEventLogEnd(dsdp->ctime);  break; }
      dsdp->pnorm=ppnorm;
      info=DSDPYStepLineSearch2(dsdp, mutarget, dsdp->dstep, dsdp->dy);DSDPCHKERR(info);
      DSDPEventLogEnd(dsdp->ctime);
    }
    if (attempt>0)dsdp->dstep=1.0;
    
    dsdp->mutarget=DSDPMin(dsdp->mu,mutarget);

    info=DSDPGetRR(dsdp,&dd1);DSDPCHKERR(info);
    if (dsdp->itnow==0 && dsdp->xmaker[0].pstep<1.0 && dd1> 0 && dsdp->pstep<1.0 && dsdp->goty0==DSDP_FALSE){
      info=DSDPResetY0(dsdp);DSDPCHKERR(info); continue;
      dsdp->goty0=DSDP_FALSE;
    }
    
  } /* End of Dual Scaling Algorithm */
  
  
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DSDPChooseBarrierParameter"
/*!
\fn int DSDPChooseBarrierParameter(DSDP dsdp, double mutarget, double *ppstep, double *nextmutarget);
\param dsdp solver
\param mutarget current barrier parameter
\param ppstep set to primal step length
\param nextmutarget set to new target barrier parameter
\brief DSDP barrier heuristic choses the smalles value of mu such that X>0

This routine implements a dynamic strategy for reducing the barrier parameter.  Basically,
it looks for the smallest barrier parameter for which the primal matrix X is psd.  Lower
and upper bounds to this parameter also apply.

*/
int DSDPChooseBarrierParameter(DSDP dsdp, double mutarget, double *ppstep, double *nextmutarget){
  int attempt,info,count=0;
  double pnorm,pstep=*ppstep,pmumin,mur, dmury1, mutarget2;
  DSDPTruth psdefinite=DSDP_FALSE;

  DSDPFunctionBegin;

  *nextmutarget=mutarget;
  /* Compute a feasible primal matrix with the smallest possible value of \mu  */
  /* Start by finding a psd feasible matrix */
  if (*ppstep >=1){
    pstep=1.0;
    /* PS should alread be formed and factored */
  } else {
    mur=-1.0/mutarget;
    info=DSDPComputePDY(dsdp,mutarget,dsdp->dy,&pnorm); DSDPCHKERR(info);
    info=DSDPComputeMaxStepLength(dsdp,dsdp->dy,DUAL_FACTOR,&pstep);DSDPCHKERR(info);
    /*  pstep=DSDPMin(0.97*pstep,1.0);  */
    if (pstep<1.0) {pstep=DSDPMin(0.97*pstep,1.0);} else {pstep=DSDPMin(1.0*pstep,1.0);}
    while (psdefinite==DSDP_FALSE){
      if (count > 2 && pstep <1e-8){pstep=0;break;}
      info=DSDPComputePY(dsdp,pstep,dsdp->ytemp);DSDPCHKERR(info);
      info=DSDPComputeSS(dsdp,dsdp->ytemp,PRIMAL_FACTOR,&psdefinite);DSDPCHKERR(info);
      if (psdefinite==DSDP_FALSE){ 
	if (count>1) pstep=0.5*pstep; else pstep=0.97*pstep;
	DSDPLogInfo(0,2,"Reducing pstep: %8.8e\n",pstep);
	count++;
      }
    }
    *ppstep=pstep;
    if (pstep > dsdp->xmaker[0].pstep || mutarget < dsdp->xmaker[0].mu * 1.0e-8){
      info=DSDPSaveYForX(dsdp,mutarget,pstep);DSDPCHKERR(info);
    }
    if (pstep==0){
      DSDPFunctionReturn(0);
    }
  }

  /* Now determine how much smaller of mu can be used */
  mur=pstep/mutarget;
  info=DSDPComputePDY1(dsdp,mur,dsdp->rhstemp); DSDPCHKERR(info);

  /* Smallest value of mu that gives a positive definite matrix */
  info=DSDPComputeMaxStepLength(dsdp,dsdp->rhstemp,PRIMAL_FACTOR,&dmury1);DSDPCHKERR(info);
  dmury1 = DSDPMin(1000,0.97*dmury1);
  
  /*  We could test the point, My tests say its not necessary -- its good!  */
  attempt=0;psdefinite=DSDP_FALSE;
  pmumin=mutarget / (1.0 + 1.0 * dmury1); /* This should be positive definite */
  while ( 0  && psdefinite==DSDP_FALSE){
    pmumin=mutarget / (1.0 + 1.0 * dmury1); /* This should be positive definite */
    if (attempt>2){pmumin=mutarget;}/* We have actually factored this one.  It is PSD. */
    info=DSDPComputePDY(dsdp,pmumin,dsdp->dy,&pnorm); DSDPCHKERR(info);
    info=DSDPComputePY(dsdp,pstep,dsdp->ytemp); DSDPCHKERR(info);
    info=DSDPComputeSS(dsdp,dsdp->ytemp,PRIMAL_FACTOR,&psdefinite);DSDPCHKERR(info);
    if (psdefinite==DSDP_FALSE){ dmury1*=0.9; /* printf("NO GOOD \n"); */ }
    else { /* printf("ITS GOOD \n"); */}
    attempt++;
    DSDPLogInfo(0,2,"GOT X: Smallest Mu for feasible X: %4.4e.\n",pmumin);
  }
  
  DSDPLogInfo(0,6,"GOT X: Smallest Mu for feasible X: %4.4e \n",pmumin);
  
  mutarget2=mutarget;
  if (dsdp->pstep==1){
    mutarget2 = pmumin;
  } else {
    /*      printf("PMUMIN: %4.4e MUTARGET: %4.4e \n",pmumin,dsdp->mutarget); */
    mutarget2 = mutarget;
    mutarget2 = (pstep)*mutarget + (1.0-pstep)*dsdp->mu;
    mutarget2 = (pstep)*pmumin + (1.0-pstep)*dsdp->mu;
  }

  mutarget2=DSDPMax(dsdp->mu/dsdp->rhon,mutarget2);
  if (dsdp->mu0>0){mutarget2=DSDPMin(mutarget2,dsdp->mu0);}

  *nextmutarget=mutarget2;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPResetY0"
/*!
\fn int DSDPResetY0(DSDP dsdp);
\param dsdp solver

\brief After 1 iteration, consider increasing the variable r 
*/
int DSDPResetY0(DSDP dsdp){
  int info;
  double r,rr,dd;
  DSDPTruth psdefinite;
  DSDPFunctionBegin;
  info=DSDPComputeDY(dsdp,dsdp->mutarget,dsdp->dy,&dsdp->pnorm); DSDPCHKERR(info);
  info=DSDPVecCopy(dsdp->y0,dsdp->y);DSDPCHKERR(info);
  info=DSDPGetRR(dsdp,&r);DSDPCHKERR(info);
  rr=DSDPMax(r*10000,1e12);
  info=DSDPSetRR(dsdp,rr);DSDPCHKERR(info);
  info=DSDPComputeSS(dsdp,dsdp->y,DUAL_FACTOR,&psdefinite);DSDPCHKERR(info);
  info=DSDPComputeLogSDeterminant(dsdp,&dsdp->logdet);DSDPCHKERR(info);
  info=DSDPSetY(dsdp,1.0,dsdp->logdet,dsdp->y);DSDPCHKERR(info);
  info=DSDPVecGetR(dsdp->b,&dd);DSDPCHKERR(info);

  dsdp->mutarget=fabs(rr*dd);
  dsdp->mu=fabs(rr*dd);

  /*  dsdp->par.mu0=mutarget; */
  dsdp->goty0=DSDP_TRUE;
  DSDPLogInfo(0,2,"Restart algorithm\n");
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPComputeDualStepDirections(DSDP dsdp);

\param dsdp the solver

\brief Compute the step direction by computing a linear system and
solving it.

DSDP first attempts unpreconditioned CG to the matrix.  Once the number
of iterations becomes too large, it swithes a CG preconditioned
by the Cholesky factorization.  Usually only one iteration of the
preconditioned CG is necessary, but solutions with large norms and
very precise solutions may require additional iterations.

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeDualStepDirections"
int DSDPComputeDualStepDirections(DSDP dsdp){
  int info,computem=1;
  double madd,ymax,cgtol=1e-7;
  DSDPTruth cg1,cg2,psdefinite;
  DSDPFunctionBegin;

  if (dsdp->itnow>30) dsdp->slestype=3;
  if (dsdp->rgap<1e-3) dsdp->slestype=3;
  if (dsdp->m<40) dsdp->slestype=3;
  if (0 && dsdp->itnow>20 && dsdp->m<500) dsdp->slestype=3;
  info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
  if (dsdp->slestype==1){
    cg1=DSDP_TRUE; cg2=DSDP_TRUE;
    info=DSDPInvertS(dsdp);DSDPCHKERR(info);
    info=DSDPComputeG(dsdp,dsdp->rhstemp,dsdp->rhs1,dsdp->rhs2);DSDPCHKERR(info);
    info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs1,dsdp->dy1,cgtol,&cg1);DSDPCHKERR(info);
    if (cg1==DSDP_TRUE){info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs2,dsdp->dy2,cgtol,&cg2);DSDPCHKERR(info);}
    if (cg1==DSDP_FALSE || cg2==DSDP_FALSE) dsdp->slestype=2;
  }
  if (dsdp->slestype==2){
    cg1=DSDP_TRUE; cg2=DSDP_TRUE;
    DSDPLogInfo(0,9,"Compute Hessian\n");
    info=DSDPInvertS(dsdp);DSDPCHKERR(info);
    info=DSDPComputeHessian(dsdp,dsdp->M,dsdp->rhs1,dsdp->rhs2);DSDPCHKERR(info);
    computem=0;
    DSDPLogInfo(0,9,"Apply CG\n");
    info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs1,dsdp->dy1,cgtol,&cg1);DSDPCHKERR(info);
    if (cg1==DSDP_TRUE){info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs2,dsdp->dy2,cgtol,&cg2);DSDPCHKERR(info);}
    if (cg1==DSDP_FALSE || cg2==DSDP_FALSE) dsdp->slestype=3;
    
  }
  if (dsdp->slestype==3){
    DSDPLogInfo(0,9,"Factor Hessian\n");
    psdefinite=DSDP_FALSE;
    if (dsdp->Mshift < 1e-12 || dsdp->rgap<0.1 || dsdp->Mshift > 1e-6){
      madd=dsdp->Mshift;
    } else {
      madd=1e-13;
    }
    if (computem){
      info=DSDPInvertS(dsdp);DSDPCHKERR(info);
    }
    while (psdefinite==DSDP_FALSE){
      if (0==1 && dsdp->Mshift>dsdp->maxschurshift){ 
	info = DSDPSetConvergenceFlag(dsdp,DSDP_INDEFINITE_SCHUR_MATRIX); DSDPCHKERR(info);  
	break;
      }
      if (0 && dsdp->Mshift*ymax>dsdp->pinfeastol/10){ 
	info = DSDPSetConvergenceFlag(dsdp,DSDP_INDEFINITE_SCHUR_MATRIX); DSDPCHKERR(info);  
	break;
      }
      if (madd*ymax>dsdp->pinfeastol*1000){ 
	info = DSDPSetConvergenceFlag(dsdp,DSDP_INDEFINITE_SCHUR_MATRIX); DSDPCHKERR(info);  
	break;
      }
      if (computem){
	info=DSDPComputeHessian(dsdp,dsdp->M,dsdp->rhs1,dsdp->rhs2);DSDPCHKERR(info);}
      if (0==1){info=DSDPSchurMatView(dsdp->M);DSDPCHKERR(info);}
      info = DSDPSchurMatShiftDiagonal(dsdp->M,madd);DSDPCHKERR(info);
      info = DSDPSchurMatFactor(dsdp->M,&psdefinite); DSDPCHKERR(info);
      computem=1;
      if (psdefinite==DSDP_FALSE){ 
	madd=madd*4 + 1.0e-13;
      }
    }
    dsdp->Mshift=madd;
    if (psdefinite==DSDP_TRUE){
      info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs1,dsdp->dy1,cgtol,&cg1);DSDPCHKERR(info);
      info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs2,dsdp->dy2,cgtol,&cg2);DSDPCHKERR(info);
    }
  }
  
  DSDPFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeDualStepDirections"
int DSDPRefineStepDirection(DSDP dsdp,DSDPVec rhs, DSDPVec dy){
  int info;
  double cgtol=1e-20;
  DSDPTruth cg1;
  DSDPFunctionBegin;
  
  if (dsdp->reason==DSDP_INDEFINITE_SCHUR_MATRIX){
  } else if (dsdp->reason==DSDP_SMALL_STEPS){
  } else if (dsdp->pstep<1){
  } else {
    dsdp->slestype=4;
    info=DSDPCGSolve(dsdp,dsdp->M,rhs,dy,cgtol,&cg1);DSDPCHKERR(info);    
    dsdp->slestype=3;
  }
  
  DSDPFunctionReturn(0);
}

#include "dsdp5.h"


#undef __FUNCT__  
#define __FUNCT__ "DSDPInitializeVariables"
/*!
\fn int DSDPInitializeVariables( DSDP dsdp );
\param dsdp the solver

\brief Initialize variables and factor S
 */
int DSDPInitializeVariables( DSDP dsdp ){

  int info;
  double r0,mutarget=dsdp->mutarget,penalty;
  DSDPTruth psdefinite=DSDP_FALSE;

  DSDPFunctionBegin;
  info=DSDPGetRR(dsdp,&r0);DSDPCHKERR(info);
  dsdp->rho=dsdp->np*dsdp->rhon;
  if (r0>=0) { /* Use the specified starting point */    
    info=DSDPComputeObjective(dsdp,dsdp->y,&dsdp->ddobj);DSDPCHKERR(info);
    info=DSDPComputeSS(dsdp,dsdp->y,DUAL_FACTOR,&psdefinite);DSDPCHKERR(info);
    if (mutarget<0) mutarget=(dsdp->ppobj-dsdp->ddobj)/(dsdp->rho);
  } else {
    info=DSDPGetPenalty(dsdp,&penalty);DSDPCHKERR(info);
    r0=0.1/(1.0+dsdp->cnorm);
    while (psdefinite==DSDP_FALSE ){
      r0=r0*100;
      DSDPLogInfo(0,9,"Set Initial R0 %4.2e\n",r0);
      info=DSDPSetRR(dsdp,r0);DSDPCHKERR(info);
      info=DSDPComputeSS(dsdp,dsdp->y,DUAL_FACTOR,&psdefinite);DSDPCHKERR(info);
    }
    r0=r0*dsdp->np;
    if (dsdp->cnorm>0 && dsdp->anorm>0 && dsdp->cnorm/dsdp->anorm<1){ r0=r0/(dsdp->cnorm/dsdp->anorm);}
    dsdp->mu=r0*penalty;
    if (mutarget<0){
      mutarget=(dsdp->ppobj-dsdp->ddobj)/(dsdp->rho);
      mutarget=(dsdp->ppobj-dsdp->ddobj)/(dsdp->np);
      mutarget=r0*penalty;
    }
    DSDPLogInfo(0,9,"Set Initial R0 %4.2e\n",r0);
    info=DSDPSetRR(dsdp,r0);DSDPCHKERR(info);
    info=DSDPComputeObjective(dsdp,dsdp->y,&dsdp->ddobj);DSDPCHKERR(info);
    info=DSDPComputeSS(dsdp,dsdp->y,DUAL_FACTOR,&psdefinite);DSDPCHKERR(info);
  }
  info=DSDPComputeObjective(dsdp,dsdp->y,&dsdp->ddobj);DSDPCHKERR(info);
  if (psdefinite==DSDP_FALSE){
    info=DSDPSetConvergenceFlag(dsdp,DSDP_INFEASIBLE_START);DSDPCHKERR(info);
  } else {
    info=DSDPComputeLogSDeterminant(dsdp,&dsdp->logdet);DSDPCHKERR(info);
    info=DSDPComputePotential(dsdp,dsdp->y,dsdp->logdet,&dsdp->potential);DSDPCHKERR(info);
  }
  /* Tough guess, as all these rules suggest */
  info=DSDPSetY(dsdp,0,dsdp->logdet,dsdp->y);DSDPCHKERR(info);
  info=DSDPSaveYForX(dsdp,dsdp->xmaker[0].mu,0);DSDPCHKERR(info);
  dsdp->mutarget=mutarget;
  dsdp->pstep=0.0;
  dsdp->pnorm=0;
  /*  dsdp->par.mu0=mutarget; */
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputeAndFactorS(DSDP dsdp, DSDPTruth *psdefinite)
\brief Compute and factor the dual matrix variables.

This routine may be necessary after computing the X variables.

\param dsdp is the solver
\param psdefinite is DSDP_TRUE if the S variables are positive definite.
\sa DSDPGetY()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeAndFactorS"
int DSDPComputeAndFactorS(DSDP dsdp,DSDPTruth *psdefinite){
  int info;
  DSDPFunctionBegin;
  info=DSDPComputeSS(dsdp,dsdp->y,DUAL_FACTOR,psdefinite);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
