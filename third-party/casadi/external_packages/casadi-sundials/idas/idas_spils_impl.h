/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the common header file (private version) for the Scaled
 * Preconditioned Iterative Linear Solver modules.
 * -----------------------------------------------------------------
 */

#ifndef _IDASSPILS_IMPL_H
#define _IDASSPILS_IMPL_H

#include <idas/idas_spils.h>
#include "idas_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Types of iterative linear solvers */

#define SPILS_SPGMR   1
#define SPILS_SPBCG   2
#define SPILS_SPTFQMR 3

/* Constants */

#define IDA_SPILS_MAXL    5
#define IDA_SPILS_MAXRS   5

/*
 * -----------------------------------------------------------------
 * Types : IDASpilsMemRec, IDASpilsMem                             
 * -----------------------------------------------------------------
 */

typedef struct IDASpilsMemRec {

  int s_type;          /* type of scaled preconditioned iterative LS   */

  int  s_gstype;       /* type of Gram-Schmidt orthogonalization       */
  realtype s_sqrtN;    /* sqrt(N)                                      */
  int  s_maxl;         /* maxl = maximum dimension of the Krylov space */
  int  s_maxrs;        /* maxrs = max. number of GMRES restarts        */
  realtype s_eplifac;  /* eplifac = linear convergence factor          */
  realtype s_dqincfac; /* dqincfac = optional increment factor in Jv   */
  realtype s_epslin;   /* SpgrmSolve tolerance parameter               */

  long int s_npe;      /* npe = total number of precond calls          */   
  long int s_nli;      /* nli = total number of linear iterations      */
  long int s_nps;      /* nps = total number of psolve calls           */
  long int s_ncfl;     /* ncfl = total number of convergence failures  */
  long int s_nres;     /* nres = total number of calls to res          */
  long int s_njtimes;  /* njtimes = total number of calls to jtimes    */

  long int s_nst0;     /* nst0 = saved nst (for performance monitor)   */   
  long int s_nni0;     /* nni0 = saved nni (for performance monitor)   */   
  long int s_nli0;     /* nli0 = saved nli (for performance monitor)   */   
  long int s_ncfn0;    /* ncfn0 = saved ncfn (for performance monitor) */   
  long int s_ncfl0;    /* ncfl0 = saved ncfl (for performance monitor) */   
  long int s_nwarn;    /* nwarn = no. of warnings (for perf. monitor)  */   

  N_Vector s_ytemp;    /* temp vector used by IDAAtimesDQ              */ 
  N_Vector s_yptemp;   /* temp vector used by IDAAtimesDQ              */ 
  N_Vector s_xx;       /* temp vector used by the solve function       */
  N_Vector s_ycur;     /* current y vector in Newton iteration         */
  N_Vector s_ypcur;    /* current yp vector in Newton iteration        */
  N_Vector s_rcur;     /* rcur = F(tn, ycur, ypcur)                    */

  void *s_spils_mem;   /* memory used by the generic solver            */

  long int s_last_flag; /* last error return flag                      */

  /* Preconditioner computation
   * (a) user-provided:
   *     - pdata == user_data
   *     - pfree == NULL (the user dealocates memory for f_data)
   * (b) internal preconditioner module
   *     - pdata == ida_mem
   *     - pfree == set by the prec. module and called in IDASpilsFree
   */

  IDASpilsPrecSetupFn s_pset;
  IDASpilsPrecSolveFn s_psolve;
  void (*s_pfree)(IDAMem IDA_mem);
  void *s_pdata;
  
  /* Jacobian times vector compuation
   * (a) jtimes function provided by the user:
   *     - jdata == user_data
   *     - jtimesDQ == FALSE
   * (b) internal jtimes
   *     - jdata == ida_mem
   *     - jtimesDQ == TRUE
   */

  booleantype s_jtimesDQ;
  IDASpilsJacTimesVecFn s_jtimes;
  void *s_jdata;

} *IDASpilsMem;


/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

/* Atimes and PSolve routines called by generic solver */

int IDASpilsAtimes(void *ida_mem, N_Vector v, N_Vector z);

int IDASpilsPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximation for Jac times vector */

int IDASpilsDQJtimes(realtype tt,
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     N_Vector v, N_Vector Jv, 
                     realtype c_j, void *data, 
                     N_Vector work1, N_Vector work2);



/*
 * -----------------------------------------------------------------
 * Error and Warning Messages
 * -----------------------------------------------------------------
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSGS_TIME "at t = %Lg, "
#define MSGS_FRMT "%Le."

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSGS_TIME "at t = %lg, "
#define MSGS_FRMT "%le."

#else

#define MSGS_TIME "at t = %g, "
#define MSGS_FRMT "%e."

#endif


/* Error Messages */

#define MSGS_IDAMEM_NULL   "Integrator memory is NULL."
#define MSGS_MEM_FAIL      "A memory request failed."
#define MSGS_BAD_NVECTOR   "A required vector operation is not implemented."
#define MSGS_BAD_LSTYPE    "Incompatible linear solver type."
#define MSGS_LMEM_NULL     "Linear solver memory is NULL."
#define MSGS_BAD_GSTYPE    "gstype has an illegal value."
#define MSGS_NEG_MAXRS     "maxrs < 0 illegal."
#define MSGS_NEG_EPLIFAC   "eplifac < 0.0 illegal."
#define MSGS_NEG_DQINCFAC  "dqincfac < 0.0 illegal."

#define MSGS_PSET_FAILED "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."

/* Warning Messages */

#define MSGS_WARN  "Warning: " MSGS_TIME "poor iterative algorithm performance. "

#define MSGS_AVD_WARN  MSGS_WARN "Average number of linear iterations is " MSGS_FRMT
#define MSGS_CFN_WARN  MSGS_WARN "Nonlinear convergence failure rate is " MSGS_FRMT
#define MSGS_CFL_WARN  MSGS_WARN "Linear convergence failure rate is " MSGS_FRMT

/* 
 * -----------------------------------------------------------------
 * PART II - backward problems
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Types : IDASpilsMemRecB, IDASpilsMemB       
 * -----------------------------------------------------------------
 * IDASpgmrB, IDASpbcgB, and IDASptfqmr attach such a structure to the 
 * lmemB filed of IDAadjMem
 * -----------------------------------------------------------------
 */

typedef struct IDASpilsMemRecB {

  IDASpilsJacTimesVecFnB s_jtimesB;
  IDASpilsJacTimesVecFnBS s_jtimesBS;
  IDASpilsPrecSetupFnB s_psetB;
  IDASpilsPrecSetupFnBS s_psetBS;
  IDASpilsPrecSolveFnB s_psolveB;
  IDASpilsPrecSolveFnBS s_psolveBS;
  void *s_P_dataB;

} *IDASpilsMemB;

/*
 * -----------------------------------------------------------------
 * Error Messages 
 * -----------------------------------------------------------------
 */

#define MSGS_LMEMB_NULL "Linear solver memory is NULL for the backward integration."
#define MSGS_BAD_T      "Bad t for interpolation."
#define MSGS_BAD_WHICH  "Illegal value for which."
#define MSGS_NO_ADJ     "Illegal attempt to call before calling IDAAdjInit."

#define MSGS_LMEMB_NULL  "Linear solver memory is NULL for the backward integration."


#ifdef __cplusplus
}
#endif

#endif
