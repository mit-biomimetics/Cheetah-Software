/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * -----------------------------------------------------------------
 * Programmer(s): Carol S. Woodward @ LLNL
 *    Based on kinsol_spgmr.c
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
 * preconditioned Flexible GMRES linear solver, KINSpfgmr.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "kinsol_impl.h"
#include <kinsol/kinsol_spfgmr.h>
#include "kinsol_spils_impl.h"

#include <sundials/sundials_spfgmr.h>
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

/* KINSpfgmr linit, lsetup, lsolve, and lfree routines */

static int KINSpfgmrInit(KINMem kin_mem);
static int KINSpfgmrSetup(KINMem kin_mem);
static int KINSpfgmrSolve(KINMem kin_mem, N_Vector xx, 
			  N_Vector bb, realtype *sJpnorm, realtype *sFdotJp);
static void KINSpfgmrFree(KINMem kin_mem);

/*
 * -----------------------------------------------------------------
 * readability replacements
 * -----------------------------------------------------------------
 */
/*
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
*/

/*
 * -----------------------------------------------------------------
 * Function : KINSpfgmr
 * -----------------------------------------------------------------
 * This routine allocates and initializes the memory record and
 * sets function fields specific to the SPFGMR linear solver module.
 * KINSpfgmr sets the kin_linit, kin_lsetup, kin_lsolve, and
 * kin_lfree fields in *kinmem to be KINSpfgmrInit, KINSpfgmrSetup,
 * KINSpfgmrSolve, and KINSpfgmrFree, respectively. It allocates
 * memory for a structure of type KINSpilsMemRec and sets the
 * kin_lmem field in *kinmem to the address of this structure. It
 * also calls SpfgmrMalloc to allocate memory for the module
 * SPFGMR. In summary, KINSpfgmr sets various fields in the
 * KINSpilsMemRec structure.
 * -----------------------------------------------------------------
 */

int KINSpfgmr(void *kinmem, int maxl)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  SpfgmrMem spfgmr_mem;
  int maxl1;

  if (kinmem == NULL){
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS", "KINSpfgmr", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);  
  }
  kin_mem = (KINMem) kinmem;

  /* check for required vector operations */

  /* Note: do NOT need to check for N_VLinearSum, N_VProd, N_VScale, N_VDiv, 
     or N_VWL2Norm because they are required by KINSOL */

  if ((kin_mem->kin_vtemp1->ops->nvconst == NULL) ||
      (kin_mem->kin_vtemp1->ops->nvdotprod == NULL) ||
      (kin_mem->kin_vtemp1->ops->nvl1norm == NULL)) {
    KINProcessError(NULL, KINSPILS_ILL_INPUT, "KINSPILS", "KINSpfgmr", MSGS_BAD_NVECTOR);
    return(KINSPILS_ILL_INPUT);
  }

  if (kin_mem->kin_lfree != NULL) kin_mem->kin_lfree(kin_mem);

  /* set four main function fields in kin_mem */

  kin_mem->kin_linit  = KINSpfgmrInit; 
  kin_mem->kin_lsetup = KINSpfgmrSetup;
  kin_mem->kin_lsolve = KINSpfgmrSolve;
  kin_mem->kin_lfree  = KINSpfgmrFree;

  /* get memory for KINSpilsMemRec */
  kinspils_mem = NULL;
  kinspils_mem = (KINSpilsMem) malloc(sizeof(struct KINSpilsMemRec));
  if (kinspils_mem == NULL){
    KINProcessError(NULL, KINSPILS_MEM_FAIL, "KINSPILS", "KINSpfgmr", MSGS_MEM_FAIL);
    return(KINSPILS_MEM_FAIL);  
  }

  /* Set ILS type */
  kinspils_mem->s_type = SPILS_SPFGMR;

  /* set SPFGMR parameters that were passed in call sequence */

  maxl1 = (maxl <= 0) ? KINSPILS_MAXL : maxl;
  kinspils_mem->s_maxl = maxl1;  

  /* Set defaults for Jacobian-related fields */

  kinspils_mem->s_jtimesDQ = TRUE;
  kinspils_mem->s_jtimes   = NULL;
  kinspils_mem->s_J_data   = NULL;

  /* Set defaults for preconditioner-related fields */

  kinspils_mem->s_pset   = NULL;
  kinspils_mem->s_psolve = NULL;
  kinspils_mem->s_pfree  = NULL;
  kinspils_mem->s_P_data = kin_mem->kin_user_data;

  /* Set default values for the rest of the SPFGMR parameters */

  kinspils_mem->s_pretype   = PREC_NONE;
  kinspils_mem->s_gstype    = MODIFIED_GS;
  kinspils_mem->s_maxlrst   = 0;
  kinspils_mem->s_last_flag = KINSPILS_SUCCESS;

  /* Call SpfgmrMalloc to allocate workspace for SPFGMR */

  /* vtemp1 passed as template vector */
  spfgmr_mem = NULL;
  spfgmr_mem = SpfgmrMalloc(maxl1, kin_mem->kin_vtemp1);
  if (spfgmr_mem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_FAIL, "KINSPILS", "KINSpfgmr", MSGS_MEM_FAIL);
    free(kinspils_mem); kinspils_mem = NULL;
    return(KINSPILS_MEM_FAIL);
  }

  /* This is an iterative linear solver */

  kin_mem->kin_inexact_ls = TRUE;

  /* Attach SPFGMR memory to spils memory structure */
  kinspils_mem->s_spils_mem = (void *) spfgmr_mem;

  /* attach linear solver memory to KINSOL memory */
  kin_mem->kin_lmem = kinspils_mem;

  return(KINSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * additional readability replacements
 * -----------------------------------------------------------------
 */

/*
#define maxl    (kinspils_mem->s_maxl)
#define maxlrst (kinspils_mem->s_maxlrst)
#define pset    (kinspils_mem->s_pset)
#define psolve  (kinspils_mem->s_psolve)
#define P_data  (kinspils_mem->s_P_data)
*/

/*
 * -----------------------------------------------------------------
 * Function : KINSpfgmrInit
 * -----------------------------------------------------------------
 * This routine initializes variables associated with the FGMRES
 * linear solver. Memory allocation was done previously in
 * KINSpfgmr.
 * -----------------------------------------------------------------
 */

static int KINSpfgmrInit(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;

  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* initialize counters */
  
  kinspils_mem->s_npe = 0;
  kinspils_mem->s_nli = 0;
  kinspils_mem->s_nps = 0;
  kinspils_mem->s_ncfl = 0;
  kinspils_mem->s_njtimes = 0;
  kinspils_mem->s_nfes = 0;

  /* set preconditioner type */

  if (kinspils_mem->s_psolve != NULL) {
    kinspils_mem->s_pretype = PREC_RIGHT;
  } else {
    kinspils_mem->s_pretype = PREC_NONE;
  }
  
  /* set setupNonNull to TRUE iff there is preconditioning with setup */

  kin_mem->kin_setupNonNull = 
    (kinspils_mem->s_psolve != NULL) && (kinspils_mem->s_pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (kinspils_mem->s_jtimesDQ) {
    kinspils_mem->s_jtimes = KINSpilsDQJtimes;
    kinspils_mem->s_J_data = kin_mem;
  } else {
    kinspils_mem->s_J_data = kin_mem->kin_user_data;
  }

  if ( (kin_mem->kin_globalstrategy == KIN_PICARD) 
       && (kinspils_mem->s_jtimesDQ) ) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSpfgmrInit", 
		    MSG_NOL_FAIL);
    return(KIN_ILL_INPUT);
  }

  kinspils_mem->s_last_flag = KINSPILS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpfgmrSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the SPFGMR linear
 * solver, that is, it is an interface to the user-supplied
 * routine pset.
 * -----------------------------------------------------------------
 */

static int KINSpfgmrSetup(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;
  int ret;

  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* call pset routine */

  ret = kinspils_mem->s_pset(kin_mem->kin_uu, kin_mem->kin_uscale, 
			     kin_mem->kin_fval, kin_mem->kin_fscale, 
			     kinspils_mem->s_P_data, 
			     kin_mem->kin_vtemp1, kin_mem->kin_vtemp2); 

  kinspils_mem->s_last_flag = ret;

  (kinspils_mem->s_npe)++;
  kin_mem->kin_nnilset = kin_mem->kin_nni; 

  /* return the same value ret that pset returned */
  return(ret);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpfgmrSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic SPFGMR solver
 * SpfgmrSolve for the solution of the linear system Ax = b.
 *
 * Appropriate variables are passed to SpfgmrSolve and the counters
 * nli, nps, and ncfl are incremented, and the return value is set
 * according to the success of SpfgmrSolve. The success flag is
 * returned if SpfgmrSolve converged, or if the residual was reduced.
 * Of the other error conditions, only preconditioner solver
 * failure is specifically returned. Otherwise a generic flag is
 * returned to denote failure of this routine.
 * -----------------------------------------------------------------
 */

static int KINSpfgmrSolve(KINMem kin_mem, N_Vector xx, N_Vector bb, 
                         realtype *sJpnorm, realtype *sFdotJp)
{
  KINSpilsMem kinspils_mem;
  SpfgmrMem spfgmr_mem;
  int ret, nli_inc, nps_inc;
  realtype res_norm;
  
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  spfgmr_mem = (SpfgmrMem) kinspils_mem->s_spils_mem;

  /* Set initial guess to xx = 0. bb is set, by the routine
     calling KINSpfgmrSolve, to the RHS vector for the system
     to be solved. */ 
  N_VConst(ZERO, xx);

  kinspils_mem->s_new_uu = TRUE; /* set flag required for user Jacobian routine */

  /* call SpfgmrSolve */
  ret = SpfgmrSolve(spfgmr_mem, kin_mem, xx, bb, kinspils_mem->s_pretype, 
		    kinspils_mem->s_gstype, kin_mem->kin_eps, 
		    kinspils_mem->s_maxlrst, kinspils_mem->s_maxl, 
		    kin_mem, kin_mem->kin_fscale, kin_mem->kin_fscale, 
		    KINSpilsAtimes, KINSpilsPSolve, &res_norm, &nli_inc, &nps_inc);

  /* increment counters nli, nps, and ncfl 
     (nni is updated in the KINSol main iteration loop) */

  kinspils_mem->s_nli = kinspils_mem->s_nli + (long int) nli_inc;
  kinspils_mem->s_nps = kinspils_mem->s_nps + (long int) nps_inc;

  if (kin_mem->kin_printfl > 2) 
    KINPrintInfo(kin_mem, PRNT_NLI, "KINSPFGMR", "KINSpfgmrSolve", INFO_NLI, nli_inc);

  if (ret != 0) (kinspils_mem->s_ncfl)++;
  kinspils_mem->s_last_flag = ret;

  if ( (ret != 0) && (ret != SPFGMR_RES_REDUCED) ) {

    /* Handle all failure returns from SpfgmrSolve */

    switch(ret) {
    case SPFGMR_PSOLVE_FAIL_REC:
    case SPFGMR_ATIMES_FAIL_REC:
      return(1);
      break;
    case SPFGMR_CONV_FAIL:
    case SPFGMR_QRFACT_FAIL:
    case SPFGMR_MEM_NULL:
    case SPFGMR_GS_FAIL:
    case SPFGMR_QRSOL_FAIL:
    case SPFGMR_ATIMES_FAIL_UNREC:
    case SPFGMR_PSOLVE_FAIL_UNREC:
      return(-1);
      break;
    }
  }

  /*  SpfgmrSolve returned either SPFGMR_SUCCESS or SPFGMR_RES_REDUCED.

     Compute the terms sJpnorm and sFdotJp for use in the linesearch
     routine and in KINForcingTerm.  Both of these terms are subsequently
     corrected if the step is reduced by constraints or the linesearch.

     sJpnorm is the norm of the scaled product (scaled by fscale) of the
     current Jacobian matrix J and the step vector p (= solution vector xx).

     sFdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale.                            */

  ret = KINSpilsAtimes(kin_mem, xx, bb);
  if (ret > 0) {
    kinspils_mem->s_last_flag = SPFGMR_ATIMES_FAIL_REC;
    return(1);
  }      
  else if (ret < 0) {
    kinspils_mem->s_last_flag = SPFGMR_ATIMES_FAIL_UNREC;
    return(-1);
  }

  *sJpnorm = N_VWL2Norm(bb, kin_mem->kin_fscale);
  N_VProd(bb, kin_mem->kin_fscale, bb);
  N_VProd(bb, kin_mem->kin_fscale, bb);
  *sFdotJp = N_VDotProd(kin_mem->kin_fval, bb);

  if (kin_mem->kin_printfl > 2) KINPrintInfo(kin_mem, PRNT_EPS, "KINSPFGMR",
					     "KINSpfgmrSolve", INFO_EPS, 
					     res_norm, kin_mem->kin_eps);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpfgmrFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the SPFGMR linear solver.
 * -----------------------------------------------------------------
 */

static void KINSpfgmrFree(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;
  SpfgmrMem spfgmr_mem;

  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  spfgmr_mem = (SpfgmrMem) kinspils_mem->s_spils_mem;
  SpfgmrFree(spfgmr_mem);

  if (kinspils_mem->s_pfree != NULL) (kinspils_mem->s_pfree)(kin_mem);

  free(kinspils_mem); kinspils_mem = NULL;
}
