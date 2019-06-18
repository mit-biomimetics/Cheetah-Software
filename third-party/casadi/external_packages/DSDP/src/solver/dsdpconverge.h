#ifndef __DSDP_CONVERGE_H
#define __DSDP_CONVERGE_H
/*! 
\file dsdpconverge.h
\brief Detect convergence of the solver from the duality gap and step sizes.
*/
#define DSDPHistory 200
typedef struct {
  int    history;    /*  length of records for iter, alpha, relgap    */    
  double alpha[DSDPHistory];    /* History of stepsize                */
  double gaphist[DSDPHistory];  /* History of duality gap             */
  double infhist[DSDPHistory];  /* History of dual infeasiblity       */

  double steptol;
  double rgaptol;
  double pnormtol;
  double dualbound;  
} ConvergenceMonitor;
#endif
