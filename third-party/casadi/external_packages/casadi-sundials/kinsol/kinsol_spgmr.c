/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
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
 * This is the implementation file for the KINSOL scaled,
 * preconditioned GMRES linear solver, KINSpgmr.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "kinsol_impl.h"
#include <kinsol/kinsol_spgmr.h>
#include "kinsol_spils_impl.h"

#include <sundials/sundials_spgmr.h>
#include <sundials/sundials_math.h>

/*
 * -----------------------------------------------------------------
 * private constants
 * -----------------------------------------------------------------
 */

#define ZERO RCONST(0.0)

/*
 * -----------------------------------------------------------------
 * function prototypes
 * -----------------------------------------------------------------
 */

/* KINSpgmr linit, lsetup, lsolve, and lfree routines */

static int KINSpgmrInit(KINMem kin_mem);
static int KINSpgmrSetup(KINMem kin_mem);
static int KINSpgmrSolve(KINMem kin_mem, N_Vector xx, 
                         N_Vector bb, realtype *sJpnorm, realtype *sFdotJp);
static void KINSpgmrFree(KINMem kin_mem);

/*
 * -----------------------------------------------------------------
 * readability replacements
 * -----------------------------------------------------------------
 */

#define nni            (kin_mem->kin_nni)
#define nnilset        (kin_mem->kin_nnilset)
#define func           (kin_mem->kin_func)
#define user_data      (kin_mem->kin_user_data)
#define printfl        (kin_mem->kin_printfl)
#define linit          (kin_mem->kin_linit)
#define lsetup         (kin_mem->kin_lsetup)
#define lsolve         (kin_mem->kin_lsolve)
#define lfree          (kin_mem->kin_lfree)
#define lmem           (kin_mem->kin_lmem)
#define inexact_ls     (kin_mem->kin_inexact_ls)
#define uu             (kin_mem->kin_uu)
#define fval           (kin_mem->kin_fval)
#define uscale         (kin_mem->kin_uscale)
#define fscale         (kin_mem->kin_fscale)
#define sqrt_relfunc   (kin_mem->kin_sqrt_relfunc)
#define jacCurrent     (kin_mem->kin_jacCurrent)
#define eps            (kin_mem->kin_eps)
#define errfp          (kin_mem->kin_errfp)
#define infofp         (kin_mem->kin_infofp)
#define setupNonNull   (kin_mem->kin_setupNonNull)
#define vtemp1         (kin_mem->kin_vtemp1)
#define vec_tmpl       (kin_mem->kin_vtemp1)
#define vtemp2         (kin_mem->kin_vtemp2)
#define strategy       (kin_mem->kin_globalstrategy)

#define pretype   (kinspils_mem->s_pretype)
#define gstype    (kinspils_mem->s_gstype)
#define nli       (kinspils_mem->s_nli)
#define npe       (kinspils_mem->s_npe)
#define nps       (kinspils_mem->s_nps)
#define ncfl      (kinspils_mem->s_ncfl)
#define njtimes   (kinspils_mem->s_njtimes)
#define nfes      (kinspils_mem->s_nfes)
#define new_uu    (kinspils_mem->s_new_uu)
#define spils_mem (kinspils_mem->s_spils_mem)

#define jtimesDQ  (kinspils_mem->s_jtimesDQ)
#define jtimes    (kinspils_mem->s_jtimes)
#define J_data    (kinspils_mem->s_J_data)

#define last_flag (kinspils_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : KINSpgmr
 * -----------------------------------------------------------------
 * This routine allocates and initializes the memory record and
 * sets function fields specific to the SPGMR linear solver module.
 * KINSpgmr sets the kin_linit, kin_lsetup, kin_lsolve, and
 * kin_lfree fields in *kinmem to be KINSpgmrInit, KINSpgmrSetup,
 * KINSpgmrSolve, and KINSpgmrFree, respectively. It allocates
 * memory for a structure of type KINSpilsMemRec and sets the
 * kin_lmem field in *kinmem to the address of this structure. It
 * also calls SpgmrMalloc to allocate memory for the module
 * SPGMR. In summary, KINSpgmr sets various fields in the
 * KINSpilsMemRec structure.
 * -----------------------------------------------------------------
 */

int KINSpgmr(void *kinmem, int maxl)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  SpgmrMem spgmr_mem;
  int maxl1;

  if (kinmem == NULL){
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS", "KINSpgmr", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);  
  }
  kin_mem = (KINMem) kinmem;

  /* check for required vector operations */

  /* Note: do NOT need to check for N_VLinearSum, N_VProd, N_VScale, N_VDiv, 
     or N_VWL2Norm because they are required by KINSOL */

  if ((vec_tmpl->ops->nvconst == NULL) ||
      (vec_tmpl->ops->nvdotprod == NULL) ||
      (vec_tmpl->ops->nvl1norm == NULL)) {
    KINProcessError(NULL, KINSPILS_ILL_INPUT, "KINSPILS", "KINSpgmr", MSGS_BAD_NVECTOR);
    return(KINSPILS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* set four main function fields in kin_mem */

  linit  = KINSpgmrInit; 
  lsetup = KINSpgmrSetup;
  lsolve = KINSpgmrSolve;
  lfree  = KINSpgmrFree;

  /* get memory for KINSpilsMemRec */
  kinspils_mem = NULL;
  kinspils_mem = (KINSpilsMem) malloc(sizeof(struct KINSpilsMemRec));
  if (kinspils_mem == NULL){
    KINProcessError(NULL, KINSPILS_MEM_FAIL, "KINSPILS", "KINSpgmr", MSGS_MEM_FAIL);
    return(KINSPILS_MEM_FAIL);  
  }

  /* Set ILS type */
  kinspils_mem->s_type = SPILS_SPGMR;

  /* set SPGMR parameters that were passed in call sequence */

  maxl1 = (maxl <= 0) ? KINSPILS_MAXL : maxl;
  kinspils_mem->s_maxl = maxl1;  

  /* Set defaults for Jacobian-related fileds */

  jtimesDQ = TRUE;
  jtimes   = NULL;
  J_data   = NULL;

  /* Set defaults for preconditioner-related fields */

  kinspils_mem->s_pset   = NULL;
  kinspils_mem->s_psolve = NULL;
  kinspils_mem->s_pfree  = NULL;
  kinspils_mem->s_P_data = kin_mem->kin_user_data;

  /* Set default values for the rest of the SPGMR parameters */

  kinspils_mem->s_pretype   = PREC_NONE;
  kinspils_mem->s_gstype    = MODIFIED_GS;
  kinspils_mem->s_maxlrst   = 0;
  kinspils_mem->s_last_flag = KINSPILS_SUCCESS;

  /* Call SpgmrMalloc to allocate workspace for SPGMR */

  /* vec_tmpl passed as template vector */
  spgmr_mem = NULL;
  spgmr_mem = SpgmrMalloc(maxl1, vec_tmpl);
  if (spgmr_mem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_FAIL, "KINSPILS", "KINSpgmr", MSGS_MEM_FAIL);
    free(kinspils_mem); kinspils_mem = NULL;
    return(KINSPILS_MEM_FAIL);
  }

  /* This is an iterative linear solver */

  inexact_ls = TRUE;

  /* Attach SPGMR memory to spils memory structure */
  spils_mem = (void *) spgmr_mem;

  /* attach linear solver memory to KINSOL memory */
  lmem = kinspils_mem;

  return(KINSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * additional readability replacements
 * -----------------------------------------------------------------
 */

#define maxl    (kinspils_mem->s_maxl)
#define maxlrst (kinspils_mem->s_maxlrst)
#define pset    (kinspils_mem->s_pset)
#define psolve  (kinspils_mem->s_psolve)
#define P_data  (kinspils_mem->s_P_data)

/*
 * -----------------------------------------------------------------
 * Function : KINSpgmrInit
 * -----------------------------------------------------------------
 * This routine initializes variables associated with the GMRES
 * linear solver. Memory allocation was done previously in
 * KINSpgmr.
 * -----------------------------------------------------------------
 */

static int KINSpgmrInit(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;

  kinspils_mem = (KINSpilsMem) lmem;

  /* initialize counters */

  npe = nli = nps = ncfl = 0;
  njtimes = nfes = 0;

  /* set preconditioner type */

  if (psolve != NULL) {
    pretype = PREC_RIGHT;
  } else {
    pretype = PREC_NONE;
  }
  
  /* set setupNonNull to TRUE iff there is preconditioning with setup */

  setupNonNull = (psolve != NULL) && (pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */

  if (jtimesDQ) {
    jtimes = KINSpilsDQJtimes;
    J_data = kin_mem;
  } else {
    J_data = user_data;
  }

  if ( (strategy == KIN_PICARD) && jtimesDQ ) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSpgmrInit", 
		    MSG_NOL_FAIL);
    return(KIN_ILL_INPUT);
  }

  last_flag = KINSPILS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpgmrSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the SPGMR linear
 * solver, that is, it is an interface to the user-supplied
 * routine pset.
 * -----------------------------------------------------------------
 */

static int KINSpgmrSetup(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;
  int ret;

  kinspils_mem = (KINSpilsMem) lmem;

  /* call pset routine */

  ret = pset(uu, uscale, fval, fscale, P_data, vtemp1, vtemp2); 

  last_flag = ret;

  npe++;
  nnilset = nni; 

  /* return the same value ret that pset returned */

  return(ret);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpgmrSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic SPGMR solver
 * SpgmrSolve for the solution of the linear system Ax = b.
 *
 * Appropriate variables are passed to SpgmrSolve and the counters
 * nli, nps, and ncfl are incremented, and the return value is set
 * according to the success of SpgmrSolve. The success flag is
 * returned if SpgmrSolve converged, or if the residual was reduced.
 * Of the other error conditions, only preconditioner solver
 * failure is specifically returned. Otherwise a generic flag is
 * returned to denote failure of this routine.
 * -----------------------------------------------------------------
 */

static int KINSpgmrSolve(KINMem kin_mem, N_Vector xx, N_Vector bb, 
                         realtype *sJpnorm, realtype *sFdotJp)
{
  KINSpilsMem kinspils_mem;
  SpgmrMem spgmr_mem;
  int ret, nli_inc, nps_inc;
  realtype res_norm;
  
  kinspils_mem = (KINSpilsMem) lmem;

  spgmr_mem = (SpgmrMem) spils_mem;

  /* Set initial guess to xx = 0. bb is set, by the routine
     calling KINSpgmrSolve, to the RHS vector for the system
     to be solved. */ 
 
  N_VConst(ZERO, xx);

  new_uu = TRUE;  /* set flag required for user Jacobian routine */

  /* call SpgmrSolve */

  ret = SpgmrSolve(spgmr_mem, kin_mem, xx, bb, pretype, gstype, eps, 
                   maxlrst, kin_mem, fscale, fscale, KINSpilsAtimes,
                   KINSpilsPSolve, &res_norm, &nli_inc, &nps_inc);

  /* increment counters nli, nps, and ncfl 
     (nni is updated in the KINSol main iteration loop) */

  nli = nli + (long int) nli_inc;
  nps = nps + (long int) nps_inc;

  if (printfl > 2) 
    KINPrintInfo(kin_mem, PRNT_NLI, "KINSPGMR", "KINSpgmrSolve", INFO_NLI, nli_inc);

  if (ret != 0) ncfl++;
  last_flag = ret;

  if ( (ret != 0) && (ret != SPGMR_RES_REDUCED) ) {

    /* Handle all failure returns from SpgmrSolve */

    switch(ret) {
    case SPGMR_PSOLVE_FAIL_REC:
    case SPGMR_ATIMES_FAIL_REC:
      return(1);
      break;
    case SPGMR_CONV_FAIL:
    case SPGMR_QRFACT_FAIL:
    case SPGMR_MEM_NULL:
    case SPGMR_GS_FAIL:
    case SPGMR_QRSOL_FAIL:
    case SPGMR_ATIMES_FAIL_UNREC:
    case SPGMR_PSOLVE_FAIL_UNREC:
      return(-1);
      break;
    }
  }

  /*  SpgmrSolve returned either SPGMR_SUCCESS or SPGMR_RES_REDUCED.

     Compute the terms sJpnorm and sFdotJp for use in the linesearch
     routine and in KINForcingTerm.  Both of these terms are subsequently
     corrected if the step is reduced by constraints or the linesearch.

     sJpnorm is the norm of the scaled product (scaled by fscale) of the
     current Jacobian matrix J and the step vector p (= solution vector xx).

     sFdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale.                            */

  ret = KINSpilsAtimes(kin_mem, xx, bb);
  if (ret > 0) {
    last_flag = SPGMR_ATIMES_FAIL_REC;
    return(1);
  }      
  else if (ret < 0) {
    last_flag = SPGMR_ATIMES_FAIL_UNREC;
    return(-1);
  }

  *sJpnorm = N_VWL2Norm(bb, fscale);
  N_VProd(bb, fscale, bb);
  N_VProd(bb, fscale, bb);
  *sFdotJp = N_VDotProd(fval, bb);

  if (printfl > 2) KINPrintInfo(kin_mem, PRNT_EPS, "KINSPGMR",
                     "KINSpgmrSolve", INFO_EPS, res_norm, eps);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpgmrFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the SPGMR linear solver.
 * -----------------------------------------------------------------
 */

static void KINSpgmrFree(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;
  SpgmrMem spgmr_mem;

  kinspils_mem = (KINSpilsMem) lmem;

  spgmr_mem = (SpgmrMem) spils_mem;
  SpgmrFree(spgmr_mem);

  if (kinspils_mem->s_pfree != NULL) (kinspils_mem->s_pfree)(kin_mem);

  free(kinspils_mem); kinspils_mem = NULL;
}
