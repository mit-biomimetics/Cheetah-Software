#include "dsdp5.h"
/*!
  \file dsdpprintout.c
  \brief Print iteration statistics.
*/

static int dsdpprintlevel=0;
static int dsdpprintlevel2=0;

#undef __FUNCT__  
#define __FUNCT__ "DSDPPrintStats"
int DSDPPrintStatsFile(DSDP dsdp, void *dummy){
  
  double ppobj,ddobj,pstp,dstp,mu,res,pinfeas,pnorm;
  int iter,info;
  int printlevel=dsdpprintlevel2;
  DSDPTerminationReason reason;

  if(printlevel<=0) return(0);
  if(!dsdpoutputfile) return(0);
  
  info = DSDPStopReason(dsdp,&reason);DSDPCHKERR(info);
  info = DSDPGetIts(dsdp,&iter);DSDPCHKERR(info);
  
  if( (reason!=CONTINUE_ITERATING) || ((iter % printlevel)==0)){
    info = DSDPGetDDObjective(dsdp,&ddobj); DSDPCHKERR(info);
    info = DSDPGetPPObjective(dsdp,&ppobj); DSDPCHKERR(info);
    info = DSDPGetR(dsdp,&res); DSDPCHKERR(info);
    info = DSDPGetPInfeasibility(dsdp,&pinfeas); DSDPCHKERR(info);
    info = DSDPGetStepLengths(dsdp,&pstp,&dstp); DSDPCHKERR(info);
    info = DSDPGetBarrierParameter(dsdp,&mu); DSDPCHKERR(info);
    info = DSDPGetPnorm(dsdp,&pnorm); DSDPCHKERR(info);
    if (reason==CONTINUE_ITERATING && iter>100 && iter%10!=0) return 0;

    if (iter==0){
      fprintf(dsdpoutputfile,"Iter   PP Objective      DD Objective    PInfeas  DInfeas     Mu     StepLength   Pnrm\n");
      fprintf(dsdpoutputfile,"--------------------------------------------------------------------------------------\n");
    }
    fprintf(dsdpoutputfile,"%-3d %16.8e  %16.8e %9.1e %9.1e %9.1e",iter,ppobj,ddobj,pinfeas,res,mu);
    fprintf(dsdpoutputfile,"  %4.2f  %4.2f",pstp,dstp);
    if (pnorm>1.0e3){
      fprintf(dsdpoutputfile,"  %1.0e \n",pnorm);
    } else {
      fprintf(dsdpoutputfile,"  %5.2f \n",pnorm);
    }
    
  }
  return 0;
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPSetStandardMonitor"
int DSDPSetFileMonitor(DSDP dsdp, int printlevel){
  int info;
  dsdpprintlevel2=printlevel;
  info=DSDPSetMonitor(dsdp,DSDPPrintStatsFile,0); DSDPCHKERR(info);
  return (0);
}

/*!
\fn int DSDPPrintStats(DSDP dsdp, void *ctx)
\brief Print statistics about the current solution to standard output.

\param dsdp is the solver
\param ctx is a pointer to a structure (NULL in this case)
\sa DSDPSetStandardMonitor()
\ingroup DSDPSolver
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPPrintStats"
int DSDPPrintStats(DSDP dsdp, void *dummy){
  
  double ppobj,ddobj,pstp,dstp,mu,res,pinfeas,pnorm;
  int iter,info;
  int printlevel=dsdpprintlevel;
  DSDPTerminationReason reason;

  if(printlevel<=0) return(0);

  info = DSDPStopReason(dsdp,&reason);DSDPCHKERR(info);
  info = DSDPGetIts(dsdp,&iter);DSDPCHKERR(info);

  if( (reason!=CONTINUE_ITERATING) || ((iter % printlevel)==0)){
    info = DSDPGetDDObjective(dsdp,&ddobj); DSDPCHKERR(info);
    info = DSDPGetPPObjective(dsdp,&ppobj); DSDPCHKERR(info);
    info = DSDPGetR(dsdp,&res); DSDPCHKERR(info);
    info = DSDPGetPInfeasibility(dsdp,&pinfeas); DSDPCHKERR(info);
    info = DSDPGetStepLengths(dsdp,&pstp,&dstp); DSDPCHKERR(info);
    info = DSDPGetBarrierParameter(dsdp,&mu); DSDPCHKERR(info);
    info = DSDPGetPnorm(dsdp,&pnorm); DSDPCHKERR(info);
    if (0 && reason==CONTINUE_ITERATING && iter>100 && iter%10!=0) return 0;

    if (iter==0){
      printf("Iter   PP Objective      DD Objective    PInfeas   DInfeas     Nu     StepLength   Pnrm\n")
	;
      printf("---------------------------------------------------------------------------------------\n")
	;
    }
    printf("%-3d %16.8e  %16.8e %9.1e %9.1e %9.1e",iter,ppobj,ddobj,pinfeas,res,mu);
    printf("  %4.2f  %4.2f",pstp,dstp);
    if (pnorm>1.0e3){
      printf("  %1.0e \n",pnorm);
    } else {
      printf("  %5.2f \n",pnorm);
    }
    fflush(NULL);
  }
  return 0;
}

/*!
\fn int DSDPSetStandardMonitor(DSDP dsdp, int k)
\brief Print at every kth iteration.

\param dsdp is the solver
\param k is the frequency to print information.
\ingroup DSDPBasic
\verbatim
Iter   PP Objective      DD Objective    PInfeas   DInfeas     Nu     StepLength   Pnrm
---------------------------------------------------------------------------------------
0     1.00000000e+02   -1.13743137e+05   2.2e+00   3.8e+02   1.1e+05  0.00  0.00   0.00
1     1.36503342e+06   -6.65779055e+04   5.1e+00   2.2e+02   1.1e+04  1.00  0.33   4.06
2     1.36631922e+05   -6.21604409e+03   5.4e+00   1.9e+01   4.5e+02  1.00  1.00   7.85
3     5.45799174e+03   -3.18292092e+03   1.5e-03   9.1e+00   7.5e+01  1.00  1.00  17.63
4     1.02930559e+03   -5.39166166e+02   1.1e-05   5.3e-01   2.7e+01  1.00  1.00   7.58
5     4.30074471e+02   -3.02460061e+01   3.3e-09   0.0e+00   5.6e+00  1.00  1.00  11.36
...
11    8.99999824e+00    8.99999617e+00   1.1e-16   0.0e+00   1.7e-08  1.00  1.00   7.03
12    8.99999668e+00    8.99999629e+00   2.9e-19   0.0e+00   3.4e-09  1.00  1.00  14.19
\endverbatim

- \c Iter - the current iteration number,
- \c PP \c Objective - the current objective value in (PP),
- \c DD \c Objective - the current objective value in (DD),
- \c PInfeas - is the largest number \f$(x^u - x^l)_i\f$ in (PP),
- \c DInfeas - the variable r in (DD) that corresponds to the infeasibility of y and S in (D).
- \c Nu - the current barrier parameter \f$\frac{\bar{z} - b^Ty}{\rho}\f$.  This parameter decreases
to zero as the points get closer to the solution,
- \c StepLength - the multiple of the step-directions in (PP) and (DD),
- \c Pnrm - the proximity to a point on the central path: \f$\sqrt{-\nabla \psi^T(y^k,\bar{z}^k ) \Delta y}\f$.

\sa DSDPGetIts()
\sa DSDPGetDDObjective()
\sa DSDPGetPPObjective()
\sa DSDPGetR()
\sa DSDPGetPInfeasibility()
\sa DSDPGetBarrierParameter()
\sa DSDPGetStepLengths()
\sa DSDPGetPnorm()
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetStandardMonitor"
int DSDPSetStandardMonitor(DSDP dsdp, int k){
  int info;
  info=DSDPSetMonitor(dsdp,DSDPPrintStats,0); DSDPCHKERR(info);
  dsdpprintlevel=k;
  return (0);
}

