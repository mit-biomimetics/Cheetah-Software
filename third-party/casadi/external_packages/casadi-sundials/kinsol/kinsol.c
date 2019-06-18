/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, Carol Woodward,
 *                and Aaron Collier @ LLNL
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
 * This is the implementation file for the main KINSol solver.
 * It is independent of the KINSol linear solver in use.
 * -----------------------------------------------------------------
 *
 * EXPORTED FUNCTIONS
 * ------------------
 *   Creation and allocation functions
 *     KINCreate
 *     KINInit
 *   Main solver function
 *     KINSol
 *   Deallocation function
 *     KINFree
 *
 * PRIVATE FUNCTIONS
 * -----------------
 *     KINCheckNvector
 *   Memory allocation/deallocation
 *     KINAllocVectors
 *     KINFreeVectors
 *   Initial setup
 *     KINSolInit
 *   Step functions
 *     KINLinSolDrv
 *     KINFullNewton
 *     KINLineSearch
 *     KINConstraint
 *     KINFP
 *     KINPicardAA
 *   Stopping tests
 *     KINStop
 *     KINForcingTerm
 *   Norm functions
 *     KINScFNorm
 *     KINScSNorm
 *   KINSOL Verbose output functions
 *     KINPrintInfo
 *     KINInfoHandler
 *   KINSOL Error Handling functions
 *     KINProcessError
 *     KINErrHandler
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <math.h>

#include "kinsol_impl.h"
#include "kinsol_spils_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * MACRO DEFINITIONS
 * =================================================================
 */

/* Macro: loop */
#define loop for(;;)

/* 
 * =================================================================
 * KINSOL PRIVATE CONSTANTS
 * =================================================================
 */

#define HALF      RCONST(0.5)
#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define ONEPT5    RCONST(1.5)
#define TWO       RCONST(2.0)
#define THREE     RCONST(3.0)
#define FIVE      RCONST(5.0)
#define TWELVE    RCONST(12.0)
#define POINT1    RCONST(0.1)
#define POINT01   RCONST(0.01)
#define POINT99   RCONST(0.99)
#define THOUSAND  RCONST(1000.0)
#define ONETHIRD  RCONST(0.3333333333333333)
#define TWOTHIRDS RCONST(0.6666666666666667)
#define POINT9    RCONST(0.9)
#define POINT0001 RCONST(0.0001)

/* 
 * =================================================================
 * KINSOL ROUTINE-SPECIFIC CONSTANTS
 * =================================================================
 */

/* 
 * Control constants for lower-level functions used by KINSol 
 * ----------------------------------------------------------
 *
 * KINStop return value requesting more iterations
 *    RETRY_ITERATION
 *    CONTINUE_ITERATIONS
 *
 * KINFullNewton, KINLineSearch, KINFP, and KINPicardAA return values:
 *    KIN_SUCCESS
 *    KIN_SYSFUNC_FAIL
 *    STEP_TOO_SMALL
 *
 * KINConstraint return values:
 *    KIN_SUCCESS
 *    CONSTR_VIOLATED
 */

#define RETRY_ITERATION     -998
#define CONTINUE_ITERATIONS -999
#define STEP_TOO_SMALL      -997
#define CONSTR_VIOLATED     -996

/*
 * Algorithmic constants
 * ---------------------
 *
 * MAX_RECVR   max. no. of attempts to correct a recoverable func error
 */

#define MAX_RECVR      5

/*
 * Keys for KINPrintInfo
 * ---------------------
 */

#define PRNT_RETVAL     1
#define PRNT_NNI        2
#define PRNT_TOL        3
#define PRNT_FMAX       4
#define PRNT_PNORM      5
#define PRNT_PNORM1     6
#define PRNT_FNORM      7
#define PRNT_LAM        8
#define PRNT_ALPHA      9
#define PRNT_BETA      10
#define PRNT_ALPHABETA 11
#define PRNT_ADJ       12

/* 
 * =================================================================
 * PRIVATE FUNCTION PROTOTYPES
 * =================================================================
 */

static booleantype KINCheckNvector(N_Vector tmpl);
static booleantype KINAllocVectors(KINMem kin_mem, N_Vector tmpl);
static int KINSolInit(KINMem kin_mem);
static int KINConstraint(KINMem kin_mem );
static void KINForcingTerm(KINMem kin_mem, realtype fnormp);
static void KINFreeVectors(KINMem kin_mem);

static int  KINFullNewton(KINMem kin_mem, realtype *fnormp, 
                          realtype *f1normp, booleantype *maxStepTaken);
static int  KINLineSearch(KINMem kin_mem, realtype *fnormp, 
                          realtype *f1normp, booleantype *maxStepTaken);
static int  KINPicardAA(KINMem kin_mem, long int *iter, realtype *R, 
			realtype *gamma, realtype *fmax);
static int  KINFP(KINMem kin_mem, long int *iter, realtype *R, 
		  realtype *gamma, realtype *fmax);

static int  KINLinSolDrv(KINMem kinmem);
static int  KINPicardFcnEval(KINMem kin_mem, N_Vector gval, N_Vector uval, 
			     N_Vector fval1);
static realtype KINScFNorm(KINMem kin_mem, N_Vector v, N_Vector scale);
static realtype KINScSNorm(KINMem kin_mem, N_Vector v, N_Vector u);
static int KINStop(KINMem kin_mem, booleantype maxStepTaken, 
		   int sflag);
static int AndersenAcc(KINMem kin_mem, N_Vector gval, N_Vector fv, N_Vector x, 
		       N_Vector x_old, int iter, realtype *R, realtype *gamma);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * Creation and allocation functions
 * -----------------------------------------------------------------
 */

/*
 * Function : KINCreate
 *
 * KINCreate creates an internal memory block for a problem to 
 * be solved by KINSOL. If successful, KINCreate returns a pointer
 * to the problem memory. This pointer should be passed to
 * KINInit. If an initialization error occurs, KINCreate prints
 * an error message to standard error and returns NULL. 
 */

void *KINCreate(void)
{
  KINMem kin_mem;
  realtype uround;

  kin_mem = NULL;
  kin_mem = (KINMem) malloc(sizeof(struct KINMemRec));
  if (kin_mem == NULL) {
    KINProcessError(kin_mem, 0, "KINSOL", "KINCreate", MSG_MEM_FAIL);
    return(NULL);
  }

  /* Zero out kin_mem */
  memset(kin_mem, 0, sizeof(struct KINMemRec));

  /* set uround (unit roundoff) */

  kin_mem->kin_uround = uround = UNIT_ROUNDOFF;
  
  /* set default values for solver optional inputs */

  kin_mem->kin_func             = NULL;
  kin_mem->kin_user_data        = NULL;
  kin_mem->kin_constraints      = NULL;
  kin_mem->kin_uscale           = NULL;
  kin_mem->kin_fscale           = NULL;
  kin_mem->kin_fold_aa          = NULL;
  kin_mem->kin_gold_aa          = NULL;
  kin_mem->kin_df_aa            = NULL;
  kin_mem->kin_dg_aa            = NULL;
  kin_mem->kin_q_aa             = NULL;
  kin_mem->kin_qtmp_aa          = NULL;
  kin_mem->kin_gamma_aa         = NULL;
  kin_mem->kin_R_aa             = NULL;
  kin_mem->kin_m_aa             = ZERO;
  kin_mem->kin_aamem_aa         = 0;
  kin_mem->kin_setstop_aa       = 0;
  kin_mem->kin_constraintsSet   = FALSE;
  kin_mem->kin_ehfun            = KINErrHandler;
  kin_mem->kin_eh_data          = kin_mem;
  kin_mem->kin_errfp            = stderr;
  kin_mem->kin_ihfun            = KINInfoHandler;
  kin_mem->kin_ih_data          = kin_mem;
  kin_mem->kin_infofp           = stdout;
  kin_mem->kin_printfl          = PRINTFL_DEFAULT;
  kin_mem->kin_mxiter           = MXITER_DEFAULT;
  kin_mem->kin_noInitSetup      = FALSE;
  kin_mem->kin_msbset           = MSBSET_DEFAULT;
  kin_mem->kin_noResMon         = FALSE;
  kin_mem->kin_msbset_sub       = MSBSET_SUB_DEFAULT;
  kin_mem->kin_update_fnorm_sub = FALSE;
  kin_mem->kin_mxnbcf           = MXNBCF_DEFAULT;
  kin_mem->kin_sthrsh           = TWO;
  kin_mem->kin_noMinEps         = FALSE;
  kin_mem->kin_mxnstepin        = ZERO;
  kin_mem->kin_sqrt_relfunc     = SUNRsqrt(uround);
  kin_mem->kin_scsteptol        = SUNRpowerR(uround,TWOTHIRDS);
  kin_mem->kin_fnormtol         = SUNRpowerR(uround,ONETHIRD);
  kin_mem->kin_etaflag          = KIN_ETACHOICE1;
  kin_mem->kin_eta              = POINT1;     /* default for KIN_ETACONSTANT */
  kin_mem->kin_eta_alpha        = TWO;        /* default for KIN_ETACHOICE2  */
  kin_mem->kin_eta_gamma        = POINT9;     /* default for KIN_ETACHOICE2  */
  kin_mem->kin_MallocDone       = FALSE;
  kin_mem->kin_setupNonNull     = FALSE;
  kin_mem->kin_eval_omega       = TRUE;
  kin_mem->kin_omega            = ZERO;       /* default to using min/max    */
  kin_mem->kin_omega_min        = OMEGA_MIN;
  kin_mem->kin_omega_max        = OMEGA_MAX;

  /* initialize lrw and liw */

  kin_mem->kin_lrw = 17;
  kin_mem->kin_liw = 22;

  /* NOTE: needed since KINInit could be called after KINSetConstraints */

  kin_mem->kin_lrw1 = 0;
  kin_mem->kin_liw1 = 0;

  return((void *) kin_mem);
}

#define errfp (kin_mem->kin_errfp)
#define liw   (kin_mem->kin_liw)
#define lrw   (kin_mem->kin_lrw)

/*
 * Function : KINInit
 *
 * KINInit allocates memory for a problem or execution of KINSol. 
 * If memory is successfully allocated, KIN_SUCCESS is returned.
 * Otherwise, an error message is printed and an error flag
 * returned.
 */

int KINInit(void *kinmem, KINSysFn func, N_Vector tmpl)
{
  long int liw1, lrw1;
  KINMem kin_mem;
  booleantype allocOK, nvectorOK;
  
  /* check kinmem */

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINInit", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (func == NULL) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINInit", MSG_FUNC_NULL);
    return(KIN_ILL_INPUT);
  }

  /* check if all required vector operations are implemented */

  nvectorOK = KINCheckNvector(tmpl);
  if (!nvectorOK) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINInit", MSG_BAD_NVECTOR);
    return(KIN_ILL_INPUT);
  }

  /* set space requirements for one N_Vector */

  if (tmpl->ops->nvspace != NULL) {
    N_VSpace(tmpl, &lrw1, &liw1);
    kin_mem->kin_lrw1 = lrw1;
    kin_mem->kin_liw1 = liw1;
  }
  else {
    kin_mem->kin_lrw1 = 0;
    kin_mem->kin_liw1 = 0;
  }

  /* allocate necessary vectors */

  allocOK = KINAllocVectors(kin_mem, tmpl);
  if (!allocOK) {
    KINProcessError(kin_mem, KIN_MEM_FAIL, "KINSOL", "KINInit", MSG_MEM_FAIL);
    free(kin_mem); kin_mem = NULL;
    return(KIN_MEM_FAIL);
  }

  /* copy the input parameter into KINSol state */

  kin_mem->kin_func = func;

  /* set the linear solver addresses to NULL */

  kin_mem->kin_linit  = NULL;
  kin_mem->kin_lsetup = NULL;
  kin_mem->kin_lsolve = NULL;
  kin_mem->kin_lfree  = NULL;
  kin_mem->kin_lmem   = NULL;
  
  /* problem memory has been successfully allocated */

  kin_mem->kin_MallocDone = TRUE;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Readability constants
 * -----------------------------------------------------------------
 */

#define func             (kin_mem->kin_func)
#define user_data        (kin_mem->kin_user_data)
#define printfl          (kin_mem->kin_printfl)
#define mxiter           (kin_mem->kin_mxiter)
#define noInitSetup      (kin_mem->kin_noInitSetup)
#define retry_nni        (kin_mem->kin_retry_nni)
#define msbset           (kin_mem->kin_msbset)
#define etaflag          (kin_mem->kin_etaflag)
#define eta              (kin_mem->kin_eta)
#define ealpha           (kin_mem->kin_eta_alpha)
#define egamma           (kin_mem->kin_eta_gamma)
#define noMinEps         (kin_mem->kin_noMinEps)
#define mxnewtstep       (kin_mem->kin_mxnewtstep)
#define mxnstepin        (kin_mem->kin_mxnstepin)
#define mxnbcf           (kin_mem->kin_mxnbcf)
#define relfunc          (kin_mem->kin_sqrt_relfunc)
#define fnormtol         (kin_mem->kin_fnormtol)
#define scsteptol        (kin_mem->kin_scsteptol)
#define constraints      (kin_mem->kin_constraints)

#define uround           (kin_mem->kin_uround)
#define nni              (kin_mem->kin_nni)
#define nfe              (kin_mem->kin_nfe)
#define nbcf             (kin_mem->kin_nbcf)  
#define nbktrk           (kin_mem->kin_nbktrk)
#define ncscmx           (kin_mem->kin_ncscmx)
#define stepl            (kin_mem->kin_stepl)
#define stepmul          (kin_mem->kin_stepmul)
#define sthrsh           (kin_mem->kin_sthrsh)
#define linit            (kin_mem->kin_linit)
#define lsetup           (kin_mem->kin_lsetup)
#define lsolve           (kin_mem->kin_lsolve) 
#define lfree            (kin_mem->kin_lfree)
#define constraintsSet   (kin_mem->kin_constraintsSet) 
#define jacCurrent       (kin_mem->kin_jacCurrent)          
#define nnilset          (kin_mem->kin_nnilset)
#define lmem             (kin_mem->kin_lmem)        
#define inexact_ls       (kin_mem->kin_inexact_ls)
#define setupNonNull     (kin_mem->kin_setupNonNull)
#define fval             (kin_mem->kin_fval)      
#define fnorm            (kin_mem->kin_fnorm)
#define f1norm           (kin_mem->kin_f1norm)
#define etaflag          (kin_mem->kin_etaflag)
#define callForcingTerm  (kin_mem->kin_callForcingTerm)
#define uu               (kin_mem->kin_uu)
#define uscale           (kin_mem->kin_uscale)
#define fscale           (kin_mem->kin_fscale)
#define sJpnorm          (kin_mem->kin_sJpnorm)
#define sFdotJp          (kin_mem->kin_sFdotJp)
#define unew             (kin_mem->kin_unew)
#define pp               (kin_mem->kin_pp)
#define vtemp1           (kin_mem->kin_vtemp1)
#define vtemp2           (kin_mem->kin_vtemp2)
#define eps              (kin_mem->kin_eps)
#define liw1             (kin_mem->kin_liw1)
#define lrw1             (kin_mem->kin_lrw1)

#define noResMon         (kin_mem->kin_noResMon)
#define fnorm_sub        (kin_mem->kin_fnorm_sub)
#define msbset_sub       (kin_mem->kin_msbset_sub)
#define nnilset_sub      (kin_mem->kin_nnilset_sub)
#define update_fnorm_sub (kin_mem->kin_update_fnorm_sub)
#define eval_omega       (kin_mem->kin_eval_omega)
#define omega            (kin_mem->kin_omega)
#define omega_min        (kin_mem->kin_omega_min)
#define omega_max        (kin_mem->kin_omega_max)

#define fold             (kin_mem->kin_fold_aa)
#define gold             (kin_mem->kin_gold_aa)
#define df               (kin_mem->kin_df_aa)
#define dg               (kin_mem->kin_dg_aa)
#define Q                (kin_mem->kin_q_aa)
#define qtmp             (kin_mem->kin_qtmp_aa)
#define maa              (kin_mem->kin_m_aa)
#define aamem            (kin_mem->kin_aamem_aa)
#define setstop          (kin_mem->kin_setstop_aa)
#define strategy         (kin_mem->kin_globalstrategy)

/* 
 * -----------------------------------------------------------------
 * Main solver function
 * -----------------------------------------------------------------
 */

/*
 * Function : KINSol
 *
 * KINSol (main KINSOL driver routine) manages the computational
 * process of computing an approximate solution of the nonlinear
 * system F(uu) = 0. The KINSol routine calls the following
 * subroutines:
 *
 *  KINSolInit    checks if initial guess satisfies user-supplied
 *                constraints and initializes linear solver
 *
 *  KINLinSolDrv  interfaces with linear solver to find a
 *                solution of the system J(uu)*x = b (calculate
 *                Newton step)
 *
 *  KINFullNewton/KINLineSearch  implement the global strategy
 *
 *  KINForcingTerm  computes the forcing term (eta)
 *
 *  KINStop  determines if an approximate solution has been found
 */

int KINSol(void *kinmem, N_Vector u, int strategy_in,  
           N_Vector u_scale, N_Vector f_scale)
{
  realtype fnormp, f1normp, epsmin, fmax=ZERO;
  KINMem kin_mem;
  int ret, sflag;
  booleantype maxStepTaken;

  /* intialize to avoid compiler warning messages */

  maxStepTaken = FALSE;
  f1normp = fnormp = -ONE;

  /* initialize epsmin to avoid compiler warning message */

  epsmin = ZERO;

  /* check for kinmem non-NULL */

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSol", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if(kin_mem->kin_MallocDone == FALSE) {
    KINProcessError(NULL, KIN_NO_MALLOC, "KINSOL", "KINSol", MSG_NO_MALLOC);    
    return(KIN_NO_MALLOC);
  }

  /* load input arguments */

  uu = u;
  uscale = u_scale;
  fscale = f_scale;
  strategy = strategy_in;

  /* CSW:  
     Call fixed point solver if requested.  Note that this should probably
     be forked off to a FPSOL solver instead of kinsol in the future. */
  if ( strategy == KIN_FP ) {
    if (uu == NULL) {
      KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSol", MSG_UU_NULL);
      return(KIN_ILL_INPUT);
    }

    if (kin_mem->kin_constraintsSet != FALSE) {
      KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSol", MSG_CONSTRAINTS_NOTOK);
      return(KIN_ILL_INPUT);
    }

    if (printfl > 0)
      KINPrintInfo(kin_mem, PRNT_TOL, "KINSOL", "KINSol", INFO_TOL, scsteptol, fnormtol);

    nfe = nnilset = nnilset_sub = nni = nbcf = nbktrk = 0;
    ret = KINFP(kin_mem, &nni, kin_mem->kin_R_aa, kin_mem->kin_gamma_aa, &fmax);

    switch(ret) {
    case KIN_SYSFUNC_FAIL:
      KINProcessError(kin_mem, KIN_SYSFUNC_FAIL, "KINSOL", "KINSol", MSG_SYSFUNC_FAILED);
      break;
    case KIN_MAXITER_REACHED:
      KINProcessError(kin_mem, KIN_MAXITER_REACHED, "KINSOL", "KINSol", MSG_MAXITER_REACHED);
      break;
    }

    return(ret);
  }

  /* initialize solver */
  ret = KINSolInit(kin_mem);
  if (ret != KIN_SUCCESS) return(ret);

  ncscmx = 0;

  /* Note: The following logic allows the choice of whether or not
     to force a call to the linear solver setup upon a given call to
     KINSol */

  if (noInitSetup) sthrsh = ONE;
  else             sthrsh = TWO;

  /* if eps is to be bounded from below, set the bound */

  if (inexact_ls && !noMinEps) epsmin = POINT01 * fnormtol;


  /* if omega is zero at this point, make sure it will be evaluated
     at each iteration based on the provided min/max bounds and the
     current function norm. */
  if (omega == ZERO) eval_omega = TRUE;
  else               eval_omega = FALSE;
 

  /* CSW:  
     Call fixed point solver for Picard method if requested.  Note that this should probably
     be forked off to a part of an FPSOL solver instead of kinsol in the future. */
  if ( strategy == KIN_PICARD ) {

    kin_mem->kin_gval = N_VClone(unew);
    lrw += lrw1;
    ret = KINPicardAA(kin_mem, &nni, kin_mem->kin_R_aa, kin_mem->kin_gamma_aa, &fmax);

    return(ret);
  }


  loop{

    retry_nni = FALSE;

    nni++;

    /* calculate the epsilon (stopping criteria for iterative linear solver)
       for this iteration based on eta from the routine KINForcingTerm */

    if (inexact_ls) {
      eps = (eta + uround) * fnorm;
      if(!noMinEps) eps = SUNMAX(epsmin, eps);
    }

    repeat_nni:

    /* call the appropriate routine to calculate an acceptable step pp */

    sflag = 0;

    if (strategy == KIN_NONE) {

      /* Full Newton Step*/

      /* call KINLinSolDrv to calculate the (approximate) Newton step, pp */ 
      ret = KINLinSolDrv(kin_mem);
      if (ret != KIN_SUCCESS) break;

      sflag = KINFullNewton(kin_mem, &fnormp, &f1normp, &maxStepTaken);

      /* if sysfunc failed unrecoverably, stop */
      if ((sflag == KIN_SYSFUNC_FAIL) || (sflag == KIN_REPTD_SYSFUNC_ERR)) {
        ret = sflag;
        break;
      }

    } else if (strategy == KIN_LINESEARCH) {

      /* Line Search */

      /* call KINLinSolDrv to calculate the (approximate) Newton step, pp */ 
      ret = KINLinSolDrv(kin_mem);
      if (ret != KIN_SUCCESS) break;

      sflag = KINLineSearch(kin_mem, &fnormp, &f1normp, &maxStepTaken);

      /* if sysfunc failed unrecoverably, stop */
      if ((sflag == KIN_SYSFUNC_FAIL) || (sflag == KIN_REPTD_SYSFUNC_ERR)) {
        ret = sflag;
        break;
      }

      /* if too many beta condition failures, then stop iteration */
      if (nbcf > mxnbcf) {
        ret = KIN_LINESEARCH_BCFAIL;
        break;
      }

    }

    if ( (strategy != KIN_PICARD) && (strategy != KIN_FP) ) {
      
      /* evaluate eta by calling the forcing term routine */
      if (callForcingTerm) KINForcingTerm(kin_mem, fnormp);

      fnorm = fnormp;

      /* call KINStop to check if tolerances where met by this iteration */
      ret = KINStop(kin_mem, maxStepTaken, sflag); 

      if (ret == RETRY_ITERATION) {
	retry_nni = TRUE;
	goto repeat_nni;
      }
    }

    /* update uu after the iteration */
    N_VScale(ONE, unew, uu);

    f1norm = f1normp;

    /* print the current nni, fnorm, and nfe values if printfl > 0 */

    if (printfl>0)
      KINPrintInfo(kin_mem, PRNT_NNI, "KINSOL", "KINSol", INFO_NNI, nni, nfe, fnorm);

    if (ret != CONTINUE_ITERATIONS) break; 

    fflush(errfp);
    
  }  /* end of loop; return */



  if (printfl > 0)
    KINPrintInfo(kin_mem, PRNT_RETVAL, "KINSOL", "KINSol", INFO_RETVAL, ret);

  switch(ret) {
  case KIN_SYSFUNC_FAIL:
    KINProcessError(kin_mem, KIN_SYSFUNC_FAIL, "KINSOL", "KINSol", MSG_SYSFUNC_FAILED);
    break;
  case KIN_REPTD_SYSFUNC_ERR:
    KINProcessError(kin_mem, KIN_REPTD_SYSFUNC_ERR, "KINSOL", "KINSol", MSG_SYSFUNC_REPTD);
    break;
  case KIN_LSETUP_FAIL:
    KINProcessError(kin_mem, KIN_LSETUP_FAIL, "KINSOL", "KINSol", MSG_LSETUP_FAILED);
    break;
  case KIN_LSOLVE_FAIL:
    KINProcessError(kin_mem, KIN_LSOLVE_FAIL, "KINSOL", "KINSol", MSG_LSOLVE_FAILED);
    break;
  case KIN_LINSOLV_NO_RECOVERY:
    KINProcessError(kin_mem, KIN_LINSOLV_NO_RECOVERY, "KINSOL", "KINSol", MSG_LINSOLV_NO_RECOVERY);
    break;
  case KIN_LINESEARCH_NONCONV:
    KINProcessError(kin_mem, KIN_LINESEARCH_NONCONV, "KINSOL", "KINSol", MSG_LINESEARCH_NONCONV);
    break;
  case KIN_LINESEARCH_BCFAIL:
    KINProcessError(kin_mem, KIN_LINESEARCH_BCFAIL, "KINSOL", "KINSol", MSG_LINESEARCH_BCFAIL);
    break;
  case KIN_MAXITER_REACHED:
    KINProcessError(kin_mem, KIN_MAXITER_REACHED, "KINSOL", "KINSol", MSG_MAXITER_REACHED);
    break;
  case KIN_MXNEWT_5X_EXCEEDED:
    KINProcessError(kin_mem, KIN_MXNEWT_5X_EXCEEDED, "KINSOL", "KINSol", MSG_MXNEWT_5X_EXCEEDED);
    break;
  }
  
  return(ret);
}

/* 
 * -----------------------------------------------------------------
 * Deallocation function
 * -----------------------------------------------------------------
 */

/*
 * Function : KINFree
 *
 * This routine frees the problem memory allocated by KINInit.
 * Such memory includes all the vectors allocated by
 * KINAllocVectors, and the memory lmem for the linear solver
 * (deallocated by a call to lfree).
 */

void KINFree(void **kinmem)
{
  KINMem kin_mem;

  if (*kinmem == NULL) return;

  kin_mem = (KINMem) (*kinmem);
  KINFreeVectors(kin_mem);

  /* call lfree if non-NULL */

  if (lfree != NULL) lfree(kin_mem);

  free(*kinmem);
  *kinmem = NULL;
}

/* 
 * =================================================================
 * PRIVATE FUNCTIONS
 * =================================================================
 */

/*
 * Function : KINCheckNvector
 *
 * This routine checks if all required vector operations are
 * implemented (excluding those required by KINConstraint). If all
 * necessary operations are present, then KINCheckNvector returns
 * TRUE. Otherwise, FALSE is returned.
 */

static booleantype KINCheckNvector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone     == NULL) ||
      (tmpl->ops->nvdestroy   == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) ||
      (tmpl->ops->nvprod      == NULL) ||
      (tmpl->ops->nvdiv       == NULL) ||
      (tmpl->ops->nvscale     == NULL) ||
      (tmpl->ops->nvabs       == NULL) ||
      (tmpl->ops->nvinv       == NULL) ||
      (tmpl->ops->nvmaxnorm   == NULL) ||
      (tmpl->ops->nvmin       == NULL) ||
      (tmpl->ops->nvwl2norm   == NULL)) return(FALSE);
  else return(TRUE);
}

/* 
 * -----------------------------------------------------------------
 * Memory allocation/deallocation
 * -----------------------------------------------------------------
 */

/*
 * Function : KINAllocVectors
 *
 * This routine allocates the KINSol vectors. If all memory
 * allocations are successful, KINAllocVectors returns TRUE.
 * Otherwise all allocated memory is freed and KINAllocVectors
 * returns FALSE.
 */

static booleantype KINAllocVectors(KINMem kin_mem, N_Vector tmpl)
{
  /* allocate unew, fval, pp, vtemp1 and vtemp2. */  
  /* allocate df, dg, q, for Andersen Acceleration, Broyden and EN */
 
  unew = N_VClone(tmpl);
  if (unew == NULL) return(FALSE);

  fval = N_VClone(tmpl);
  if (fval == NULL) {
    N_VDestroy(unew);
    return(FALSE);
  }

  pp = N_VClone(tmpl);
  if (pp == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    return(FALSE);
  }

  vtemp1 = N_VClone(tmpl);
  if (vtemp1 == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    N_VDestroy(pp);
    return(FALSE);
  }

  vtemp2 = N_VClone(tmpl);
  if (vtemp2 == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    N_VDestroy(pp);
    N_VDestroy(vtemp1);
    return(FALSE);
  }

  /* update solver workspace lengths */

  liw += 5*liw1;
  lrw += 5*lrw1;

  if (maa) {
    kin_mem->kin_R_aa = (realtype *) malloc((maa*maa) * sizeof(realtype));
    if (kin_mem->kin_R_aa == NULL) {
      KINProcessError(kin_mem, 0, "KINSOL", "KINAllocVectors", MSG_MEM_FAIL);
      return(KIN_MEM_FAIL);
    }
    kin_mem->kin_gamma_aa = (realtype *)malloc(maa * sizeof(realtype));
    if (kin_mem->kin_gamma_aa == NULL) {
      KINProcessError(kin_mem, 0, "KINSOL", "KINAllocVectors", MSG_MEM_FAIL);
      return(KIN_MEM_FAIL);
    }
  } 

  if (maa) {
    fold = N_VClone(tmpl);
    if (fold == NULL) {
      N_VDestroy(unew);
      N_VDestroy(fval);
      N_VDestroy(pp);
      N_VDestroy(vtemp1);
      N_VDestroy(vtemp2);
      return(FALSE);
    }
    gold = N_VClone(tmpl);
    if (gold == NULL) {
      N_VDestroy(unew);
      N_VDestroy(fval);
      N_VDestroy(pp);
      N_VDestroy(vtemp1);
      N_VDestroy(vtemp2);
      N_VDestroy(fold);
      return(FALSE);
    }
    df = N_VCloneVectorArray(maa,tmpl);
    if (df == NULL) {
      N_VDestroy(unew);
      N_VDestroy(fval);
      N_VDestroy(pp);
      N_VDestroy(vtemp1);
      N_VDestroy(vtemp2);
      N_VDestroy(fold);
      N_VDestroy(gold);
      return(FALSE);
    }
    dg = N_VCloneVectorArray(maa,tmpl);
    if (dg == NULL) {
      N_VDestroy(unew);
      N_VDestroy(fval);
      N_VDestroy(pp);
      N_VDestroy(vtemp1);
      N_VDestroy(vtemp2);
      N_VDestroy(fold);
      N_VDestroy(gold);
      N_VDestroyVectorArray(df, maa);
      return(FALSE);
    }
    /* update solver workspace lengths */

    liw += 2*maa*liw1+2;
    lrw += 2*maa*lrw1+2;

    if (aamem) {
      Q = N_VCloneVectorArray(maa,tmpl);
      if (Q == NULL) {
	N_VDestroy(unew);
	N_VDestroy(fval);
	N_VDestroy(pp);
	N_VDestroy(vtemp1);
	N_VDestroy(vtemp2);
	N_VDestroy(fold);
	N_VDestroy(gold);
	N_VDestroyVectorArray(df, maa);
	N_VDestroyVectorArray(dg, maa);
	return(FALSE);
      }
      qtmp = N_VCloneVectorArray(maa,tmpl);
      if (qtmp == NULL) {
	N_VDestroy(unew);
	N_VDestroy(fval);
	N_VDestroy(pp);
	N_VDestroy(vtemp1);
	N_VDestroy(vtemp2);
	N_VDestroy(fold);
	N_VDestroy(gold);
	N_VDestroyVectorArray(df, maa);
	N_VDestroyVectorArray(dg, maa);
	N_VDestroyVectorArray(Q, maa);
	return(FALSE);
      }
      liw += 2*maa*liw1;
      lrw += 2*maa*lrw1;
    }
  }
  return(TRUE);
}

/*
 * KINFreeVectors
 *
 * This routine frees the KINSol vectors allocated by
 * KINAllocVectors.
 */

static void KINFreeVectors(KINMem kin_mem)
{
  if (unew != NULL)   N_VDestroy(unew);
  if (fval != NULL)   N_VDestroy(fval);
  if (pp != NULL)     N_VDestroy(pp);
  if (vtemp1 != NULL) N_VDestroy(vtemp1);
  if (vtemp2 != NULL) N_VDestroy(vtemp2);

  if ( (kin_mem->kin_globalstrategy == KIN_PICARD) && (kin_mem->kin_gval != NULL) ) 
    N_VDestroy(kin_mem->kin_gval);

  if ( ((strategy == KIN_PICARD) || (strategy == KIN_FP)) && (maa > 0) ) {
    free(kin_mem->kin_R_aa);
    free(kin_mem->kin_gamma_aa);
  }

  if (maa)
  {
     if (fold != NULL) N_VDestroy(fold);
     if (gold != NULL) N_VDestroy(gold);
     N_VDestroyVectorArray(df,maa);
     N_VDestroyVectorArray(dg,maa);
     lrw -= (2*maa*lrw1+2);
     liw -= (2*maa*liw1+2);
     if (aamem)
     {
        N_VDestroyVectorArray(Q,maa);
        N_VDestroyVectorArray(qtmp,maa);
        lrw -= 2*maa*lrw1;
        liw -= 2*maa*liw1;
     }
  }

  lrw -= 5*lrw1;
  liw -= 5*liw1;

  if (kin_mem->kin_constraintsSet) {
    if (constraints != NULL) N_VDestroy(constraints);
    lrw -= lrw1;
    liw -= liw1;
  }

  return;
}

/* 
 * -----------------------------------------------------------------
 * Initial setup
 * -----------------------------------------------------------------
 */

/*
 * KINSolInit
 *
 * KINSolInit initializes the problem for the specific input
 * received in this call to KINSol (which calls KINSolInit). All
 * problem specification inputs are checked for errors. If any error
 * occurs during initialization, it is reported to the file whose
 * file pointer is errfp.
 *
 * The possible return values for KINSolInit are:
 *   KIN_SUCCESS : indicates a normal initialization
 *
 *   KIN_ILL_INPUT : indicates that an input error has been found
 *
 *   KIN_INITIAL_GUESS_OK : indicates that the guess uu
 *                          satisfied the system func(uu) = 0
 *                          within the tolerances specified
 */

static int KINSolInit(KINMem kin_mem)
{
  int retval;
  realtype fmax;
  
  /* check for illegal input parameters */

  if (uu == NULL) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_UU_NULL);
    return(KIN_ILL_INPUT);
  }

  if ( (strategy != KIN_NONE) && (strategy != KIN_LINESEARCH) && 
       (strategy != KIN_PICARD) && (strategy != KIN_FP) ) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_BAD_GLSTRAT);
    return(KIN_ILL_INPUT);
  }

  if (uscale == NULL)  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_BAD_USCALE);
    return(KIN_ILL_INPUT);
  }

  if (N_VMin(uscale) <= ZERO){
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_USCALE_NONPOSITIVE);
    return(KIN_ILL_INPUT);
  }

  if (fscale == NULL)  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_BAD_FSCALE);
    return(KIN_ILL_INPUT);
  }

  if (N_VMin(fscale) <= ZERO){
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_FSCALE_NONPOSITIVE);
    return(KIN_ILL_INPUT);
  }

  if ( (constraints != NULL) && ( (strategy == KIN_PICARD) || (strategy == KIN_FP) ) ) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_CONSTRAINTS_NOTOK);
    return(KIN_ILL_INPUT);
  }
    

  /* set the constraints flag */

  if (constraints == NULL) 
    constraintsSet = FALSE;
  else {
    constraintsSet = TRUE;
    if ((constraints->ops->nvconstrmask  == NULL) ||
	(constraints->ops->nvminquotient == NULL)) {
      KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_BAD_NVECTOR);
      return(KIN_ILL_INPUT);
    }
  }

  /* check the initial guess uu against the constraints */

  if (constraintsSet) {
    if (!N_VConstrMask(constraints, uu, vtemp1)) {
      KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSolInit", MSG_INITIAL_CNSTRNT);
      return(KIN_ILL_INPUT);
    }
  }
  
  /* all error checking is complete at this point */

  if (printfl > 0)
    KINPrintInfo(kin_mem, PRNT_TOL, "KINSOL", "KINSolInit", INFO_TOL, scsteptol, fnormtol);

  /* calculate the default value for mxnewtstep (maximum Newton step) */

  if (mxnstepin == ZERO) mxnewtstep = THOUSAND * N_VWL2Norm(uu, uscale);
  else                   mxnewtstep = mxnstepin;
  if (mxnewtstep < ONE) mxnewtstep = ONE;


  /* additional set-up for inexact linear solvers */

  if (inexact_ls) {

    /* set up the coefficients for the eta calculation */

    callForcingTerm = (etaflag != KIN_ETACONSTANT);

    /* this value is always used for choice #1 */

    if (etaflag == KIN_ETACHOICE1) ealpha = (ONE + SUNRsqrt(FIVE)) * HALF;

    /* initial value for eta set to 0.5 for other than the 
       KIN_ETACONSTANT option */

    if (etaflag != KIN_ETACONSTANT) eta = HALF;

    /* disable residual monitoring if using an inexact linear solver */

    noResMon = TRUE;

  } else {

    callForcingTerm = FALSE;

  }

  /* initialize counters */

  nfe = nnilset = nnilset_sub = nni = nbcf = nbktrk = 0;

  /* see if the initial guess uu satisfies the nonlinear system */
  retval = func(uu, fval, user_data); nfe++;

  if (retval < 0) {
    KINProcessError(kin_mem, KIN_SYSFUNC_FAIL, "KINSOL", "KINSolInit", 
		    MSG_SYSFUNC_FAILED);
    return(KIN_SYSFUNC_FAIL);
  } else if (retval > 0) {
    KINProcessError(kin_mem, KIN_FIRST_SYSFUNC_ERR, "KINSOL", "KINSolInit", 
		    MSG_SYSFUNC_FIRST);
    return(KIN_FIRST_SYSFUNC_ERR);
  }

  fmax = KINScFNorm(kin_mem, fval, fscale);
  if (fmax <= (POINT01 * fnormtol)) {
    kin_mem->kin_fnorm = N_VWL2Norm(fval, fscale);
    return(KIN_INITIAL_GUESS_OK);
  }

  if (printfl > 1)
    KINPrintInfo(kin_mem, PRNT_FMAX, "KINSOL", "KINSolInit", INFO_FMAX, fmax);
  
  /* initialize the linear solver if linit != NULL */

  if (linit != NULL) {
    retval = linit(kin_mem);
    if (retval != 0) {
      KINProcessError(kin_mem, KIN_LINIT_FAIL, "KINSOL", "KINSolInit", MSG_LINIT_FAIL);
      return(KIN_LINIT_FAIL);
    }
  }

  /* initialize the L2 (Euclidean) norms of f for the linear iteration steps */

  fnorm = N_VWL2Norm(fval, fscale);
  f1norm = HALF * fnorm * fnorm;
  fnorm_sub = fnorm;

  if (printfl > 0)
    KINPrintInfo(kin_mem, PRNT_NNI, "KINSOL", "KINSolInit", 
		 INFO_NNI, nni, nfe, fnorm);

  /* problem has now been successfully initialized */

  return(KIN_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Step functions
 * -----------------------------------------------------------------
 */

/*
 * KINLinSolDrv
 *
 * This routine handles the process of solving for the approximate
 * solution of the Newton equations in the Newton iteration.
 * Subsequent routines handle the nonlinear aspects of its
 * application. 
 */

static int KINLinSolDrv(KINMem kin_mem)
{
  N_Vector x, b;
  int retval;

  if ((nni - nnilset) >= msbset) {
    sthrsh = TWO;
    update_fnorm_sub = TRUE;
  }

  loop{

    jacCurrent = FALSE;

    if ((sthrsh > ONEPT5) && setupNonNull) {
      retval = lsetup(kin_mem);
      jacCurrent = TRUE;
      nnilset = nni;
      nnilset_sub = nni;
      if (retval != 0) return(KIN_LSETUP_FAIL);
    }

    /* rename vectors for readability */

    b = unew;
    x = pp;

    /* load b with the current value of -fval */

    N_VScale(-ONE, fval, b);

    /* call the generic 'lsolve' routine to solve the system Jx = b */

    retval = lsolve(kin_mem, x, b, &sJpnorm, &sFdotJp);

    if (retval == 0)                          return(KIN_SUCCESS);
    else if (retval < 0)                      return(KIN_LSOLVE_FAIL);
    else if ((!setupNonNull) || (jacCurrent)) return(KIN_LINSOLV_NO_RECOVERY);

    /* loop back only if the linear solver setup is in use 
       and Jacobian information is not current */

    sthrsh = TWO;

  }
}

/*
 * KINFullNewton
 *
 * This routine is the main driver for the Full Newton
 * algorithm. Its purpose is to compute unew = uu + pp in the
 * direction pp from uu, taking the full Newton step. The
 * step may be constrained if the constraint conditions are
 * violated, or if the norm of pp is greater than mxnewtstep. 
 */

static int KINFullNewton(KINMem kin_mem, realtype *fnormp, realtype *f1normp,
                         booleantype *maxStepTaken)
{
  realtype pnorm, ratio;
  booleantype fOK;
  int ircvr, retval;

  *maxStepTaken = FALSE;
  pnorm = N_VWL2Norm(pp, uscale);
  ratio = ONE;
  if (pnorm > mxnewtstep) {
    ratio = mxnewtstep / pnorm;
    N_VScale(ratio, pp, pp);
    pnorm = mxnewtstep;
  }

  if (printfl > 0)
    KINPrintInfo(kin_mem, PRNT_PNORM, "KINSOL", "KINFullNewton", INFO_PNORM, pnorm);

  /* If constraints are active, then constrain the step accordingly */

  stepl = pnorm;
  stepmul = ONE;
  if (constraintsSet) {
    retval = KINConstraint(kin_mem);
    if (retval == CONSTR_VIOLATED) {
      /* Apply stepmul set in KINConstraint */
      ratio *= stepmul;
      N_VScale(stepmul, pp, pp);
      pnorm *= stepmul;
      stepl = pnorm;
      if (printfl > 0)
        KINPrintInfo(kin_mem, PRNT_PNORM, "KINSOL", "KINFullNewton", INFO_PNORM, pnorm);
      if (pnorm <= scsteptol) {
        N_VLinearSum(ONE, uu, ONE, pp, unew);
        return(STEP_TOO_SMALL);}
    }
  }
 
  /* Attempt (at most MAX_RECVR times) to evaluate function at the new iterate */
  
  fOK = FALSE;

  for (ircvr = 1; ircvr <= MAX_RECVR; ircvr++) {

    /* compute the iterate unew = uu + pp */
    N_VLinearSum(ONE, uu, ONE, pp, unew);

    /* evaluate func(unew) and its norm, and return */
    retval = func(unew, fval, user_data); nfe++;

    /* if func was successful, accept pp */
    if (retval == 0) {fOK = TRUE; break;}

    /* if func failed unrecoverably, give up */
    else if (retval < 0) return(KIN_SYSFUNC_FAIL);

    /* func failed recoverably; cut step in half and try again */
    ratio *= HALF;
    N_VScale(HALF, pp, pp);
    pnorm *= HALF;
    stepl = pnorm;
  }

  /* If func() failed recoverably MAX_RECVR times, give up */

  if (!fOK) return(KIN_REPTD_SYSFUNC_ERR);

  /* Evaluate function norms */

  *fnormp = N_VWL2Norm(fval,fscale);
  *f1normp = HALF * (*fnormp) * (*fnormp);

  /* scale sFdotJp and sJpnorm by ratio for later use in KINForcingTerm */

  sFdotJp *= ratio;
  sJpnorm *= ratio;
 
  if (printfl > 1) 
    KINPrintInfo(kin_mem, PRNT_FNORM, "KINSOL", "KINFullNewton", INFO_FNORM, *fnormp);

  if (pnorm > (POINT99 * mxnewtstep)) *maxStepTaken = TRUE; 

  return(KIN_SUCCESS);
}

/*
 * KINLineSearch
 *
 * The routine KINLineSearch implements the LineSearch algorithm.
 * Its purpose is to find unew = uu + rl * pp in the direction pp
 * from uu so that:
 *                                    t
 *  func(unew) <= func(uu) + alpha * g  (unew - uu) (alpha = 1.e-4)
 *
 *    and
 *                                   t
 *  func(unew) >= func(uu) + beta * g  (unew - uu) (beta = 0.9)
 *
 * where 0 < rlmin <= rl <= rlmax.
 *
 * Note:
 *             mxnewtstep
 *  rlmax = ----------------   if uu+pp is feasible
 *          ||uscale*pp||_L2
 *
 *  rlmax = 1   otherwise
 *
 *    and
 *
 *                 scsteptol
 *  rlmin = --------------------------
 *          ||           pp         ||
 *          || -------------------- ||_L-infinity
 *          || (1/uscale + SUNRabs(uu)) ||
 *
 *
 * If the system function fails unrecoverably at any time, KINLineSearch 
 * returns KIN_SYSFUNC_FAIL which will halt the solver.
 *
 * We attempt to corect recoverable system function failures only before 
 * the alpha-condition loop; i.e. when the solution is updated with the 
 * full Newton step (possibly reduced due to constraint violations). 
 * Once we find a feasible pp, we assume that any update up to pp is
 * feasible.
 * 
 * If the step size is limited due to constraint violations and/or 
 * recoverable system function failures, we set rlmax=1 to ensure
 * that the update remains feasible during the attempts to enforce 
 * the beta-condition (this is not an isse while enforcing the alpha
 * condition, as rl can only decrease from 1 at that stage)
 */

static int KINLineSearch(KINMem kin_mem, realtype *fnormp, realtype *f1normp,
                         booleantype *maxStepTaken)
{
  realtype pnorm, ratio, slpi, rlmin, rlength, rl, rlmax, rldiff;
  realtype rltmp, rlprev, pt1trl, f1nprv, rllo, rlinc, alpha, beta;
  realtype alpha_cond, beta_cond, rl_a, tmp1, rl_b, tmp2, disc;
  int ircvr, nbktrk_l, retval;
  booleantype firstBacktrack, fOK;

  /* Initializations */

  nbktrk_l = 0;          /* local backtracking counter */
  ratio    = ONE;        /* step change ratio          */
  alpha    = POINT0001;
  beta     = POINT9;

  firstBacktrack = TRUE;
  *maxStepTaken = FALSE;

  rlprev = f1nprv = ZERO;

  /* Compute length of Newton step */

  pnorm = N_VWL2Norm(pp, uscale);
  rlmax = mxnewtstep / pnorm;
  stepl = pnorm;

  /* If the full Newton step is too large, set it to the maximum allowable value */

  if(pnorm > mxnewtstep ) {
    ratio = mxnewtstep / pnorm;
    N_VScale(ratio, pp, pp);
    pnorm = mxnewtstep;
    rlmax = ONE;
    stepl = pnorm;
  }

  /* If constraint checking is activated, check and correct violations */

  stepmul = ONE;

  if(constraintsSet){
    retval = KINConstraint(kin_mem);
    if(retval == CONSTR_VIOLATED){
      /* Apply stepmul set in KINConstraint */
      N_VScale(stepmul, pp, pp);
      ratio *= stepmul;
      pnorm *= stepmul;
      rlmax = ONE;
      stepl = pnorm;
      if (printfl > 0) KINPrintInfo(kin_mem, PRNT_PNORM1, "KINSOL", "KINLineSearch", INFO_PNORM1, pnorm);
      if (pnorm <= scsteptol) {
        N_VLinearSum(ONE, uu, ONE, pp, unew);
        return(STEP_TOO_SMALL);}
    }
  }

  /* Attempt (at most MAX_RECVR times) to evaluate function at the new iterate */
  
  fOK = FALSE;

  for (ircvr = 1; ircvr <= MAX_RECVR; ircvr++) {

    /* compute the iterate unew = uu + pp */
    N_VLinearSum(ONE, uu, ONE, pp, unew);

    /* evaluate func(unew) and its norm, and return */
    retval = func(unew, fval, user_data); nfe++;

    /* if func was successful, accept pp */
    if (retval == 0) {fOK = TRUE; break;}

    /* if func failed unrecoverably, give up */
    else if (retval < 0) return(KIN_SYSFUNC_FAIL);

    /* func failed recoverably; cut step in half and try again */
    N_VScale(HALF, pp, pp);
    ratio *= HALF;
    pnorm *= HALF;
    rlmax = ONE;
    stepl = pnorm;

  }

  /* If func() failed recoverably MAX_RECVR times, give up */

  if (!fOK) return(KIN_REPTD_SYSFUNC_ERR);

  /* Evaluate function norms */

  *fnormp = N_VWL2Norm(fval, fscale);
  *f1normp = HALF * (*fnormp) * (*fnormp) ;

  /* Estimate the line search value rl (lambda) to satisfy both ALPHA and BETA conditions */

  slpi = sFdotJp * ratio;
  rlength = KINScSNorm(kin_mem, pp, uu);
  rlmin = scsteptol / rlength;
  rl = ONE;

  if (printfl > 2)
    KINPrintInfo(kin_mem, PRNT_LAM, "KINSOL", "KINLineSearch", INFO_LAM, rlmin, f1norm, pnorm);

  /* Loop until the ALPHA condition is satisfied. Terminate if rl becomes too small */

  loop {
    
    /* Evaluate test quantity */

    alpha_cond = f1norm + (alpha * slpi * rl);

    if (printfl > 2)
      KINPrintInfo(kin_mem, PRNT_ALPHA, "KINSOL", "KINLinesearch", 
                   INFO_ALPHA, *fnormp, *f1normp, alpha_cond, rl);

    /* If ALPHA condition is satisfied, break out from loop */

    if ((*f1normp) <= alpha_cond) break;

    /* Backtracking. Use quadratic fit the first time and cubic fit afterwards. */

    if (firstBacktrack) {

      rltmp = -slpi / (TWO * ((*f1normp) - f1norm - slpi));
      firstBacktrack = FALSE;

    } else {

      tmp1 = (*f1normp) - f1norm - (rl * slpi);
      tmp2 = f1nprv - f1norm - (rlprev * slpi);
      rl_a = ((ONE / (rl * rl)) * tmp1) - ((ONE / (rlprev * rlprev)) * tmp2);
      rl_b = ((-rlprev / (rl * rl)) * tmp1) + ((rl / (rlprev * rlprev)) * tmp2);
      tmp1 = ONE / (rl - rlprev);
      rl_a *= tmp1;
      rl_b *= tmp1;
      disc = (rl_b * rl_b) - (THREE * rl_a * slpi);

      if (SUNRabs(rl_a) < uround) {        /* cubic is actually just a quadratic (rl_a ~ 0) */
        rltmp = -slpi / (TWO * rl_b);
      } else {                         /* real cubic */
        rltmp = (-rl_b + SUNRsqrt(disc)) / (THREE * rl_a);
      }
    }
      if (rltmp > (HALF * rl)) rltmp = HALF * rl;

    /* Set new rl (do not allow a reduction by a factor larger than 10) */

    rlprev = rl;
    f1nprv = (*f1normp);
    pt1trl = POINT1 * rl;
    rl = SUNMAX(pt1trl, rltmp);
    nbktrk_l++;

    /* Update unew and re-evaluate function */

    N_VLinearSum(ONE, uu, rl, pp, unew);

    retval = func(unew, fval, user_data); nfe++;
    if (retval != 0) return(KIN_SYSFUNC_FAIL);

    *fnormp = N_VWL2Norm(fval, fscale);
    *f1normp = HALF * (*fnormp) * (*fnormp) ;

    /* Check if rl (lambda) is too small */

    if (rl < rlmin) {
      /* unew sufficiently distinct from uu cannot be found.
         copy uu into unew (step remains unchanged) and 
         return STEP_TOO_SMALL */
      N_VScale(ONE, uu, unew);
      return(STEP_TOO_SMALL);
    }

  } /* end ALPHA condition loop */


  /* ALPHA condition is satisfied. Now check the BETA condition */

  beta_cond = f1norm + (beta * slpi * rl);

  if ((*f1normp) < beta_cond) {

    /* BETA condition not satisfied */

    if ((rl == ONE) && (pnorm < mxnewtstep)) {

      do {

        rlprev = rl;
        f1nprv = *f1normp;
        rl = SUNMIN((TWO * rl), rlmax);
        nbktrk_l++;

        N_VLinearSum(ONE, uu, rl, pp, unew);
        retval = func(unew, fval, user_data); nfe++;
        if (retval != 0) return(KIN_SYSFUNC_FAIL);
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);

        alpha_cond = f1norm + (alpha * slpi * rl);
        beta_cond = f1norm + (beta * slpi * rl);

        if (printfl > 2)
          KINPrintInfo(kin_mem, PRNT_BETA, "KINSOL", "KINLineSearch", 
                       INFO_BETA, *f1normp, beta_cond, rl);

      } while (((*f1normp) <= alpha_cond) && 
	       ((*f1normp) < beta_cond) && (rl < rlmax));

    } /* enf if (rl == ONE) block */

    if ((rl < ONE) || ((rl > ONE) && (*f1normp > alpha_cond))) {

      rllo = SUNMIN(rl, rlprev);
      rldiff = SUNRabs(rlprev - rl);

      do {

        rlinc = HALF * rldiff;
        rl = rllo + rlinc;
        nbktrk_l++;

        N_VLinearSum(ONE, uu, rl, pp, unew);
        retval = func(unew, fval, user_data); nfe++;
        if (retval != 0) return(KIN_SYSFUNC_FAIL);
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);

        alpha_cond = f1norm + (alpha * slpi * rl);
        beta_cond = f1norm + (beta * slpi * rl);

        if (printfl > 2)
          KINPrintInfo(kin_mem, PRNT_ALPHABETA, "KINSOL", "KINLineSearch", 
                       INFO_ALPHABETA, *f1normp, alpha_cond, beta_cond, rl);

        if ((*f1normp) > alpha_cond) rldiff = rlinc;
        else if (*f1normp < beta_cond) {
          rllo = rl;
          rldiff = rldiff - rlinc;
        }

      } while ((*f1normp > alpha_cond) ||
	       ((*f1normp < beta_cond) && (rldiff >= rlmin)));

      if ((*f1normp) < beta_cond) {

	/* beta condition could not be satisfied so set unew to last u value
	   that satisfied the alpha condition and continue */

        N_VLinearSum(ONE, uu, rllo, pp, unew);
        retval = func(unew, fval, user_data); nfe++;
        if (retval != 0) return(KIN_SYSFUNC_FAIL);
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);   

	/* increment beta-condition failures counter */

        nbcf++;

      }

    }  /* end of if (rl < ONE) block */

  }  /* end of if (f1normp < beta_cond) block */

  /* Update number of backtracking operations */

  nbktrk += nbktrk_l;

  if (printfl > 1)
    KINPrintInfo(kin_mem, PRNT_ADJ, "KINSOL", "KINLineSearch", INFO_ADJ, nbktrk_l);

  /* scale sFdotJp and sJpnorm by rl * ratio for later use in KINForcingTerm */

  sFdotJp = sFdotJp * rl * ratio;
  sJpnorm = sJpnorm * rl * ratio;

  if ((rl * pnorm) > (POINT99 * mxnewtstep)) *maxStepTaken = TRUE;

  return(KIN_SUCCESS);
}

/*
 * Function : KINConstraint
 *
 * This routine checks if the proposed solution vector uu + pp
 * violates any constraints. If a constraint is violated, then the
 * scalar stepmul is determined such that uu + stepmul * pp does
 * not violate any constraints.
 *
 * Note: This routine is called by the functions
 *       KINLineSearch and KINFullNewton.
 */

static int KINConstraint(KINMem kin_mem)
{
  N_VLinearSum(ONE, uu, ONE, pp, vtemp1);

  /* if vtemp1[i] violates constraint[i] then vtemp2[i] = 1
     else vtemp2[i] = 0 (vtemp2 is the mask vector) */

  if(N_VConstrMask(constraints, vtemp1, vtemp2)) return(KIN_SUCCESS);

  /* vtemp1[i] = SUNRabs(pp[i]) */

  N_VAbs(pp, vtemp1);

  /* consider vtemp1[i] only if vtemp2[i] = 1 (constraint violated) */

  N_VProd(vtemp2, vtemp1, vtemp1);

  N_VAbs(uu, vtemp2);
  stepmul = POINT9 * N_VMinQuotient(vtemp2, vtemp1);

  return(CONSTR_VIOLATED);
}

/*
 * -----------------------------------------------------------------
 * Stopping tests
 * -----------------------------------------------------------------
 */

/*
 * KINStop
 *
 * This routine checks the current iterate unew to see if the
 * system func(unew) = 0 is satisfied by a variety of tests.
 *
 * strategy is one of KIN_NONE or KIN_LINESEARCH
 * sflag    is one of KIN_SUCCESS, STEP_TOO_SMALL
 */

static int KINStop(KINMem kin_mem, booleantype maxStepTaken, int sflag)
{
  realtype fmax, rlength, omexp;
  N_Vector delta;

  /* Check for too small a step */

  if (sflag == STEP_TOO_SMALL) {

    if (setupNonNull && !jacCurrent) {
      /* If the Jacobian is out of date, update it and retry */
      sthrsh = TWO;
      return(RETRY_ITERATION);
    } else {
      /* Give up */
      if (strategy == KIN_NONE)  return(KIN_STEP_LT_STPTOL);
      else                       return(KIN_LINESEARCH_NONCONV);
    }

  }

  /* Check tolerance on scaled function norm at the current iterate */

  fmax = KINScFNorm(kin_mem, fval, fscale);

  if (printfl > 1) 
    KINPrintInfo(kin_mem, PRNT_FMAX, "KINSOL", "KINStop", INFO_FMAX, fmax);

  if (fmax <= fnormtol) return(KIN_SUCCESS);

  /* Check if the scaled distance between the last two steps is too small */
  /* NOTE: pp used as work space to store this distance */

  delta = pp;
  N_VLinearSum(ONE, unew, -ONE, uu, delta);
  rlength = KINScSNorm(kin_mem, delta, unew);

  if (rlength <= scsteptol) {

    if (setupNonNull && !jacCurrent) {
      /* If the Jacobian is out of date, update it and retry */
      sthrsh = TWO;
      return(CONTINUE_ITERATIONS);
    } else {
      /* give up */
      return(KIN_STEP_LT_STPTOL);
    }

  }

  /* Check if the maximum number of iterations is reached */

  if (nni >= mxiter) return(KIN_MAXITER_REACHED);

  /* Check for consecutive number of steps taken of size mxnewtstep
     and if not maxStepTaken, then set ncscmx to 0 */
 
  if (maxStepTaken) ncscmx++;
  else              ncscmx = 0;
 
  if (ncscmx == 5) return(KIN_MXNEWT_5X_EXCEEDED);

  /* Proceed according to the type of linear solver used */

  if (inexact_ls) {

    /* We're doing inexact Newton.
       Load threshold for reevaluating the Jacobian. */

    sthrsh = rlength;

  } else if (!noResMon) {

    /* We're doing modified Newton and the user did not disable residual monitoring.
       Check if it is time to monitor residual. */

    if ((nni - nnilset_sub) >= msbset_sub) {

      /* Residual monitoring needed */

      nnilset_sub = nni;

      /* If indicated, estimate new OMEGA value */
      if (eval_omega) {
        omexp = SUNMAX(ZERO,(fnorm/fnormtol)-ONE);
        omega = (omexp > TWELVE)? omega_max : SUNMIN(omega_min*SUNRexp(omexp), omega_max);
      }   
      /* Check if making satisfactory progress */

      if (fnorm > omega*fnorm_sub) {
        /* Insufficient progress */
	if (setupNonNull && !jacCurrent) {
          /* If the Jacobian is out of date, update it and retry */
	  sthrsh = TWO;
          return(CONTINUE_ITERATIONS);
	} else {
          /* Otherwise, we cannot do anything, so just return. */
        }
      } else {
        /* Sufficient progress */
	fnorm_sub = fnorm;
	sthrsh = ONE;
      }

    } else {

      /* Residual monitoring not needed */

      /* Reset sthrsh */
      if (retry_nni || update_fnorm_sub) fnorm_sub = fnorm;
      if (update_fnorm_sub) update_fnorm_sub = FALSE;
      sthrsh = ONE;

    }

  }

  /* if made it to here, then the iteration process is not finished
     so return CONTINUE_ITERATIONS flag */

  return(CONTINUE_ITERATIONS);
}

/*
 * KINForcingTerm
 *
 * This routine computes eta, the scaling factor in the linear
 * convergence stopping tolerance eps when choice #1 or choice #2
 * forcing terms are used. Eta is computed here for all but the
 * first iterative step, which is set to the default in routine
 * KINSolInit.
 *
 * This routine was written by Homer Walker of Utah State
 * University with subsequent modifications by Allan Taylor @ LLNL.
 *
 * It is based on the concepts of the paper 'Choosing the forcing
 * terms in an inexact Newton method', SIAM J Sci Comput, 17
 * (1996), pp 16 - 32, or Utah State University Research Report
 * 6/94/75 of the same title.
 */

static void KINForcingTerm(KINMem kin_mem, realtype fnormp)
{
  realtype eta_max, eta_min, eta_safe, linmodel_norm;

  eta_max  = POINT9;
  eta_min  = POINT0001;
  eta_safe = HALF;

  /* choice #1 forcing term */

  if (etaflag == KIN_ETACHOICE1) {

    /* compute the norm of f + Jp , scaled L2 norm */

    linmodel_norm = SUNRsqrt((fnorm * fnorm) + (TWO * sFdotJp) + (sJpnorm * sJpnorm));

    /* form the safeguarded for choice #1 */ 

    eta_safe = SUNRpowerR(eta, ealpha); 
    eta = SUNRabs(fnormp - linmodel_norm) / fnorm;
  }

  /* choice #2 forcing term */

  if (etaflag == KIN_ETACHOICE2) {
    eta_safe = egamma * SUNRpowerR(eta, ealpha); 
    eta = egamma * SUNRpowerR((fnormp / fnorm), ealpha); 
  }

  /* apply safeguards */
 
  if(eta_safe < POINT1) eta_safe = ZERO;
  eta = SUNMAX(eta, eta_safe);
  eta = SUNMAX(eta, eta_min);
  eta = SUNMIN(eta, eta_max);

  return; 
}


/*
 * -----------------------------------------------------------------
 * Norm functions
 * -----------------------------------------------------------------
 */

/*
 * Function : KINScFNorm
 *
 * This routine computes the max norm for scaled vectors. The
 * scaling vector is scale, and the vector of which the norm is to
 * be determined is vv. The returned value, fnormval, is the
 * resulting scaled vector norm. This routine uses N_Vector
 * functions from the vector module.
 */

static realtype KINScFNorm(KINMem kin_mem, N_Vector v, N_Vector scale)
{
  N_VProd(scale, v, vtemp1);
  return(N_VMaxNorm(vtemp1));
}

/*
 * Function : KINScSNorm
 *
 * This routine computes the max norm of the scaled steplength, ss.
 * Here ucur is the current step and usc is the u scale factor.
 */

static realtype KINScSNorm(KINMem kin_mem, N_Vector v, N_Vector u)
{
  realtype length;

  N_VInv(uscale, vtemp1);
  N_VAbs(u, vtemp2);
  N_VLinearSum(ONE, vtemp1, ONE, vtemp2, vtemp1);
  N_VDiv(v, vtemp1, vtemp1);

  length = N_VMaxNorm(vtemp1);

  return(length);
}

/* 
 * =================================================================
 * KINSOL Verbose output functions
 * =================================================================
 */

/* 
 * KINPrintInfo
 *
 * KINPrintInfo is a high level error handling function
 * Based on the value info_code, it composes the info message and
 * passes it to the info handler function.
 */

#define ihfun    (kin_mem->kin_ihfun)
#define ih_data  (kin_mem->kin_ih_data)

void KINPrintInfo(KINMem kin_mem, 
                  int info_code, const char *module, const char *fname, 
                  const char *msgfmt, ...)
{
  va_list ap;
  char msg[256], msg1[40];
  char retstr[30];
  int ret;

  /* Initialize argument processing 
   (msgfrmt is the last required argument) */

  va_start(ap, msgfmt); 

  if (info_code == PRNT_RETVAL) {

    /* If info_code = PRNT_RETVAL, decode the numeric value */

    ret = va_arg(ap, int);

    switch(ret) {
    case KIN_SUCCESS:
      sprintf(retstr, "KIN_SUCCESS");
      break;
    case KIN_SYSFUNC_FAIL:
      sprintf(retstr, "KIN_SYSFUNC_FAIL");
      break;
    case KIN_STEP_LT_STPTOL:
      sprintf(retstr, "KIN_STEP_LT_STPTOL");
      break;
    case KIN_LINESEARCH_NONCONV:
      sprintf(retstr, "KIN_LINESEARCH_NONCONV");
      break;
    case KIN_LINESEARCH_BCFAIL:
      sprintf(retstr, "KIN_LINESEARCH_BCFAIL");
      break;
    case KIN_MAXITER_REACHED:
      sprintf(retstr, "KIN_MAXITER_REACHED");
      break;
    case KIN_MXNEWT_5X_EXCEEDED:
      sprintf(retstr, "KIN_MXNEWT_5X_EXCEEDED");
      break;
    case KIN_LINSOLV_NO_RECOVERY:
      sprintf(retstr, "KIN_LINSOLV_NO_RECOVERY");
      break;
    case KIN_LSETUP_FAIL:
      sprintf(retstr, "KIN_PRECONDSET_FAILURE");
      break;
    case KIN_LSOLVE_FAIL:
      sprintf(retstr, "KIN_PRECONDSOLVE_FAILURE");
      break;
    }

    /* Compose the message */

    sprintf(msg1, msgfmt, ret);
    sprintf(msg,"%s (%s)",msg1,retstr);


  } else {
  
    /* Compose the message */

    vsprintf(msg, msgfmt, ap);

  }

  /* call the info message handler */

  ihfun(module, fname, msg, ih_data);

  /* finalize argument processing */

  va_end(ap);

  return;
}


/*
 * KINInfoHandler 
 *
 * This is the default KINSOL info handling function.
 * It sends the info message to the stream pointed to by kin_infofp 
 */

#define infofp (kin_mem->kin_infofp)

void KINInfoHandler(const char *module, const char *function, 
                    char *msg, void *data)
{
  KINMem kin_mem;

  /* data points to kin_mem here */

  kin_mem = (KINMem) data;
  
#ifndef NO_FPRINTF_OUTPUT
  if (infofp != NULL) {
    fprintf(infofp,"\n[%s] %s\n",module, function);
    fprintf(infofp,"   %s\n",msg);
  }
#endif  

}

/* 
 * =================================================================
 * KINSOL Error Handling functions
 * =================================================================
 */

/*
 * KINProcessError 
 *
 * KINProcessError is a high level error handling function.
 * - If cv_mem==NULL it prints the error message to stderr.
 * - Otherwise, it sets up and calls the error handling function 
 *   pointed to by cv_ehfun.
 */

#define ehfun    (kin_mem->kin_ehfun)
#define eh_data  (kin_mem->kin_eh_data)

void KINProcessError(KINMem kin_mem, 
                    int error_code, const char *module, const char *fname, 
                    const char *msgfmt, ...)
{
  va_list ap;
  char msg[256];

  /* Initialize the argument pointer variable 
     (msgfmt is the last required argument to KINProcessError) */

  va_start(ap, msgfmt);

  /* Compose the message */

  vsprintf(msg, msgfmt, ap);

  if (kin_mem == NULL) {    /* We write to stderr */
#ifndef NO_FPRINTF_OUTPUT
    fprintf(stderr, "\n[%s ERROR]  %s\n  ", module, fname);
    fprintf(stderr, "%s\n\n", msg);
#endif

  } else {                 /* We can call ehfun */
    ehfun(error_code, module, fname, msg, eh_data);
  }

  /* Finalize argument processing */
  va_end(ap);

  return;
}

/* 
 * KINErrHandler 
 *
 * This is the default error handling function.
 * It sends the error message to the stream pointed to by kin_errfp 
 */

#define errfp    (kin_mem->kin_errfp)

void KINErrHandler(int error_code, const char *module,
                   const char *function, char *msg, void *data)
{
  KINMem kin_mem;
  char err_type[10];

  /* data points to kin_mem here */

  kin_mem = (KINMem) data;

  if (error_code == KIN_WARNING)
    sprintf(err_type,"WARNING");
  else
    sprintf(err_type,"ERROR");

#ifndef NO_FPRINTF_OUTPUT
  if (errfp != NULL) {
    fprintf(errfp,"\n[%s %s]  %s\n",module,err_type,function);
    fprintf(errfp,"  %s\n\n",msg);
  }
#endif

  return;
}


/*
 * =======================================================================
 * Picard and fixed point solvers
 * =======================================================================
 */

/*
 * KINPicardAA
 *
 * This routine is the main driver for the Picard iteration with acclerated fixed point. 
 */

static int KINPicardAA(KINMem kin_mem, long int *iterp, realtype *R, realtype *gamma, 
		       realtype *fmaxptr)
{
  booleantype fOK;
  int retval, ret; 
  long int iter;
  realtype fmax, epsmin, fnormp;
  N_Vector delta, gval;
  
  delta = kin_mem->kin_vtemp1;
  gval = kin_mem->kin_gval;
  ret = CONTINUE_ITERATIONS;
  fOK = TRUE;
  fmax = fnormtol + ONE;
  iter = 0;
  epsmin = ZERO;
  fnormp = -ONE;

  N_VConst(ZERO, gval);

  /* if eps is to be bounded from below, set the bound */
  if (inexact_ls && !noMinEps) epsmin = POINT01 * fnormtol;

  while (ret == CONTINUE_ITERATIONS) {

    iter++;

    /* Update the forcing term for the inexact linear solves */
    if (inexact_ls) {
      eps = (eta + uround) * fnorm;
      if(!noMinEps) eps = SUNMAX(epsmin, eps);
    }

    /* evaluate g = uu - L^{-1}func(u) and return if failed.  
       For Picard, assume that the fval vector has been filled 
       with an eval of the nonlinear residual prior to this call. */
    retval = KINPicardFcnEval(kin_mem, gval, uu, fval);
    if (retval == 0) {
      fOK = TRUE;
    } else if (retval < 0) {
      fOK = FALSE;
      ret = KIN_SYSFUNC_FAIL;
      break;
    }

    if (maa == 0) { 
      N_VScale(ONE, gval, unew);
    }
    else {  /* use Anderson, if desired */
      N_VScale(ONE, uu, unew);
      AndersenAcc(kin_mem, gval, delta, unew, uu, (int)(iter-1), R, gamma);
    }

    /* Fill the Newton residual based on the new solution iterate */
    retval = func(unew, fval, user_data); nfe++;
    if (retval == 0) {
      fOK = TRUE;
    }
    else if (retval < 0) {
      ret = KIN_SYSFUNC_FAIL;
      break;
    }

    /* Evaluate function norms */
    fnormp = N_VWL2Norm(fval,fscale);
    fmax = KINScFNorm(kin_mem, fval, fscale); /* measure  || F(x) ||_max */
    fnorm = fmax;
    *fmaxptr = fmax;
    
    if (printfl > 1) 
      KINPrintInfo(kin_mem, PRNT_FMAX, "KINSOL", "KINPicardAA", INFO_FMAX, fmax);

    /* print the current iter, fnorm, and nfe values if printfl > 0 */
    if (printfl>0)
      KINPrintInfo(kin_mem, PRNT_NNI, "KINSOL", "KINPicardAA", INFO_NNI, iter, nfe, fnorm);

    /* Check if the maximum number of iterations is reached */
    if (iter >= mxiter) {
      ret = KIN_MAXITER_REACHED;
    }
    if (fmax <= fnormtol) { 
      ret = KIN_SUCCESS;
    }
    
    if (ret == CONTINUE_ITERATIONS) { 
      /* Only update solution if taking a next iteration.  */
      N_VScale(ONE, unew, uu);
      /* evaluate eta by calling the forcing term routine */
      if (callForcingTerm) KINForcingTerm(kin_mem, fnormp);
    }

    fflush(errfp);
    
  }  /* end of loop; return */

  *iterp = iter;

  if (printfl > 0)
    KINPrintInfo(kin_mem, PRNT_RETVAL, "KINSOL", "KINPicardAA", INFO_RETVAL, ret);

  return(ret); 
}

/*
 * KINPicardFcnEval
 *
 * This routine evaluates the Picard fixed point function
 * using the linear solver, gval = u - L^{-1}F(u).
 * The function assumes the user has defined L either through
 * a user-supplied matvec if using a SPILS solver or through
 * a supplied matrix if using a dense solver.  This assumption is
 * tested by a check on the strategy and the requisite functionality
 * within the linear solve routines.
 *
 * This routine fills gval = uu - L^{-1}F(uu) given uu and fval = F(uu).
 */ 

static int KINPicardFcnEval(KINMem kin_mem, N_Vector gval, N_Vector uval, N_Vector fval1)
{
  int retval;

  if ((nni - nnilset) >= msbset) {
    sthrsh = TWO;
    update_fnorm_sub = TRUE;
  }

  loop{

    jacCurrent = FALSE;

    if ((sthrsh > ONEPT5) && setupNonNull) {
      retval = lsetup(kin_mem);
      jacCurrent = TRUE;
      nnilset = nni;
      nnilset_sub = nni;
      if (retval != 0) return(KIN_LSETUP_FAIL);
    }

    /* call the generic 'lsolve' routine to solve the system Lx = -fval
       Note that we are using gval to hold x. */
    N_VScale(-ONE, fval1, fval1);
    retval = lsolve(kin_mem, gval, fval1, &sJpnorm, &sFdotJp);

    if (retval == 0) {
      /* Update gval = uval + gval since gval = -L^{-1}F(uu)  */
      N_VLinearSum(ONE, uval, ONE, gval, gval);
      return(KIN_SUCCESS);
    }
    else if (retval < 0)                      return(KIN_LSOLVE_FAIL);
    else if ((!setupNonNull) || (jacCurrent)) return(KIN_LINSOLV_NO_RECOVERY);

    /* loop back only if the linear solver setup is in use
       and matrix information is not current */

    sthrsh = TWO;
  }

}


/*
 * KINFP
 *
 * This routine is the main driver for the fixed point iteration with
 * Anderson Acceleration.
 */

static int KINFP(KINMem kin_mem, long int *iterp, 
		 realtype *R, realtype *gamma, 
		 realtype *fmaxptr)
{
  booleantype fOK;
  int retval, ret; 
  long int iter;
  realtype fmax;
  N_Vector delta;
  
  delta = kin_mem->kin_vtemp1;
  ret = CONTINUE_ITERATIONS;
  fOK = TRUE;
  fmax = fnormtol + ONE;
  iter = 0;

  while (ret == CONTINUE_ITERATIONS) {

    iter++;

    /* evaluate func(uu) and return if failed */
    retval = func(uu, fval, user_data); nfe++;

    if (retval < 0) {
      fOK = FALSE;
      ret = KIN_SYSFUNC_FAIL;
      break;
    }

    if (maa == 0) { 
      N_VScale(ONE, fval, unew);
    }
    else {  /* use Anderson, if desired */
      N_VScale(ONE, uu, unew);
      AndersenAcc(kin_mem, fval, delta, unew, uu, (int)(iter-1), R, gamma);
    }

    N_VLinearSum(ONE, unew, -ONE, uu, delta);
    fmax = KINScFNorm(kin_mem, delta, fscale); /* measure  || g(x)-x || */
    
    if (printfl > 1) 
      KINPrintInfo(kin_mem, PRNT_FMAX, "KINSOL", "KINFP", INFO_FMAX, fmax);
    
    fnorm = fmax;
    *fmaxptr = fmax;

    /* print the current iter, fnorm, and nfe values if printfl > 0 */
    if (printfl>0)
      KINPrintInfo(kin_mem, PRNT_NNI, "KINSOL", "KINFP", INFO_NNI, iter, nfe, fnorm);

    /* Check if the maximum number of iterations is reached */
    if (iter >= mxiter) {
      ret = KIN_MAXITER_REACHED;
    }
    if (fmax <= fnormtol) { 
      ret = KIN_SUCCESS;
    }
    
    if (ret == CONTINUE_ITERATIONS) { 
      /* Only update solution if taking a next iteration.  */
      /* CSW  Should put in a conditional to send back the newest iterate or
	 the one consistent with the fval */
      N_VScale(ONE, unew, uu);
    }

    fflush(errfp);
    
  }  /* end of loop; return */

  *iterp = iter;

  if (printfl > 0)
    KINPrintInfo(kin_mem, PRNT_RETVAL, "KINSOL", "KINFP", INFO_RETVAL, ret);

  return(ret); 
} 


 /* -----------------------------------------------------------------
 * Stopping tests
 * -----------------------------------------------------------------
 */


/*
 * ========================================================================
 * Andersen Acceleration
 * ========================================================================
 */

static int AndersenAcc(KINMem kin_mem, N_Vector gval, N_Vector fv, 
		       N_Vector x, N_Vector xold, 
		       int iter, realtype *R, realtype *gamma)
{
    
  int i_pt, i, j, lAA, imap, jmap;
  int *ipt_map;
  realtype alfa;
  
  ipt_map = (int *) malloc(maa * sizeof(int));
  i_pt = iter-1 - ((iter-1)/maa)*maa;
  N_VLinearSum(ONE, gval, -1.0, xold, fv);
  if (iter > 0) {
    /* compute dg_new = gval -gval_old*/
    N_VLinearSum(ONE, gval, -1.0, gold, dg[i_pt]);
    /* compute df_new = fval - fval_old */
    N_VLinearSum(ONE, fv, -1.0, fold, df[i_pt]);
  }
    
  N_VScale(ONE, gval, gold);
  N_VScale(ONE, fv, fold);
  
  if (iter == 0) {
    N_VScale(ONE, gval, x);
  }
  else {
    if (iter == 1) {
      N_VScale(ONE,df[i_pt],qtmp[i_pt]);
      R[0] = sqrt(N_VDotProd(df[i_pt],df[i_pt])); 
      alfa = 1/R[0];
      N_VScale(alfa,df[i_pt],Q[i_pt]);
      ipt_map[0] = 0;
    }
    else if (iter < maa) {
      N_VScale(ONE,df[i_pt],qtmp[i_pt]);
      for (j=0; j < (iter-1); j++) {
	ipt_map[j] = j;
	R[(iter-1)*maa+j] = N_VDotProd(Q[j],qtmp[i_pt]);
	N_VLinearSum(ONE,qtmp[i_pt],-R[(iter-1)*maa+j],Q[j],qtmp[i_pt]);
      }
      R[(iter-1)*maa+iter-1] = sqrt(N_VDotProd(qtmp[i_pt],qtmp[i_pt])); 
      N_VScale((1/R[(iter-1)*maa+iter-1]),qtmp[i_pt],Q[i_pt]);
      ipt_map[iter-1] = iter-1;
    }
    else {
      j = 0;
      for (i=i_pt+1; i < maa; i++)
	ipt_map[j++] = i;
      for (i=0; i < (i_pt+1); i++)
	ipt_map[j++] = i;
      
      for (i=0; i < maa; i++)
	N_VScale(ONE,df[i],qtmp[i]);
      for (i=0; i < maa; i++) {
	imap = ipt_map[i];
	R[i*maa+i] = sqrt(N_VDotProd(qtmp[imap],qtmp[imap]));
	N_VScale((1/R[i*maa+i]),qtmp[imap],Q[imap]);
	for (j = i+1; j < maa; j++) {
	  jmap = ipt_map[j];
	  R[j*maa+i] = N_VDotProd(qtmp[jmap],Q[imap]);
	  N_VLinearSum(ONE,qtmp[jmap],-R[j*maa+i],Q[imap],qtmp[jmap]);
	}            
      }
    }
    /* Solve least squares problem and update solution */
    lAA = iter;
    if (maa < iter) lAA = maa;
    N_VScale(ONE, gval, x);
    for (i=0; i < lAA; i++)
      gamma[i] = N_VDotProd(fv,Q[ipt_map[i]]);
    for (i=lAA-1; i > -1; i--) {
      for (j=i+1; j < lAA; j++) {
	gamma[i] = gamma[i]-R[j*maa+i]*gamma[j]; 
      }
      gamma[i] = gamma[i]/R[i*maa+i];
      N_VLinearSum(ONE,x,-gamma[i],dg[ipt_map[i]],x);
    }
  }
  free(ipt_map);

  return 0;
}
