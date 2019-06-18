/*
 * -----------------------------------------------------------------
 * $Revision: 4397 $
 * $Date: 2015-02-28 14:03:10 -0800 (Sat, 28 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
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
 * This is the implementation file for the KINSOL interface to the
 * scaled, preconditioned Bi-CGSTAB (SPBCG) iterative linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "kinsol_impl.h"
#include <kinsol/kinsol_spbcgs.h>
#include "kinsol_spils_impl.h"

#include <sundials/sundials_spbcgs.h>
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

/* KINSpbcg linit, lsetup, lsolve, and lfree routines */

static int KINSpbcgInit(KINMem kin_mem);
static int KINSpbcgSetup(KINMem kin_mem);
static int KINSpbcgSolve(KINMem kin_mem, N_Vector xx,
			 N_Vector bb, realtype *sJpnorm, realtype *sFdotJp);
static void KINSpbcgFree(KINMem kin_mem);

/*
 * -----------------------------------------------------------------
 * readability replacements
 * -----------------------------------------------------------------
 */

#define nni          (kin_mem->kin_nni)
#define nnilset      (kin_mem->kin_nnilset)
#define func         (kin_mem->kin_func)
#define user_data    (kin_mem->kin_user_data)
#define printfl      (kin_mem->kin_printfl)
#define linit        (kin_mem->kin_linit)
#define lsetup       (kin_mem->kin_lsetup)
#define lsolve       (kin_mem->kin_lsolve)
#define lfree        (kin_mem->kin_lfree)
#define lmem         (kin_mem->kin_lmem)
#define inexact_ls   (kin_mem->kin_inexact_ls)
#define uu           (kin_mem->kin_uu)
#define fval         (kin_mem->kin_fval)
#define uscale       (kin_mem->kin_uscale)
#define fscale       (kin_mem->kin_fscale)
#define sqrt_relfunc (kin_mem->kin_sqrt_relfunc)
#define eps          (kin_mem->kin_eps)
#define errfp        (kin_mem->kin_errfp)
#define infofp       (kin_mem->kin_infofp)
#define setupNonNull (kin_mem->kin_setupNonNull)
#define vtemp1       (kin_mem->kin_vtemp1)
#define vec_tmpl     (kin_mem->kin_vtemp1)
#define vtemp2       (kin_mem->kin_vtemp2)
#define strategy     (kin_mem->kin_globalstrategy)

#define pretype   (kinspils_mem->s_pretype)
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
 * Function : KINSpbcg
 * -----------------------------------------------------------------
 * This routine allocates and initializes the memory record and
 * sets function fields specific to the SPBCG linear solver module.
 * KINSpbcg sets the kin_linit, kin_lsetup, kin_lsolve, and
 * kin_lfree fields in *kinmem to be KINSpbcgInit, KINSpbcgSetup,
 * KINSpbcgSolve, and KINSpbcgFree, respectively. It allocates
 * memory for a structure of type KINSpilsMemRec and sets the
 * kin_lmem field in *kinmem to the address of this structure. It
 * also calls SpbcgMalloc to allocate memory for the module
 * SPBCG. It sets setupNonNull in (*kin_mem) and sets various
 * fields in the KINSpilsMemRec structure.
 * Finally, KINSpbcg allocates memory for local vectors, and calls
 * SpbcgMalloc to allocate memory for the Spbcg solver.
 * -----------------------------------------------------------------
 */

int KINSpbcg(void *kinmem, int maxl)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  SpbcgMem spbcg_mem;
  int maxl1;

  if (kinmem == NULL){
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS", "KINSpbcg", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);  
  }
  kin_mem = (KINMem) kinmem;

  /* check for required vector operations */

  /* Note: do NOT need to check for N_VLinearSum, N_VProd, N_VScale, N_VDiv, 
     or N_VWL2Norm because they are required by KINSOL */

  if ((vec_tmpl->ops->nvconst == NULL) ||
      (vec_tmpl->ops->nvdotprod == NULL) ||
      (vec_tmpl->ops->nvl1norm == NULL)) {
    KINProcessError(NULL, KINSPILS_ILL_INPUT, "KINSPILS", "KINSpbcg", MSGS_BAD_NVECTOR);
    return(KINSPILS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* set four main function fields in kin_mem */

  linit  = KINSpbcgInit; 
  lsetup = KINSpbcgSetup;
  lsolve = KINSpbcgSolve;
  lfree  = KINSpbcgFree;

  /* get memory for KINSpilsMemRec */
  kinspils_mem = NULL;
  kinspils_mem = (KINSpilsMem) malloc(sizeof(struct KINSpilsMemRec));
  if (kinspils_mem == NULL){
    KINProcessError(NULL, KINSPILS_MEM_FAIL, "KINSPILS", "KINSpbcg", MSGS_MEM_FAIL);
    return(KINSPILS_MEM_FAIL);  
  }

  /* Set ILS type */
  kinspils_mem->s_type = SPILS_SPBCG;

  /* set SPBCG parameters that were passed in call sequence */

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

  /* Set default values for the rest of the SPBCG parameters */

  kinspils_mem->s_pretype   = PREC_NONE;
  kinspils_mem->s_last_flag = KINSPILS_SUCCESS;

  /* Call SpbcgMalloc to allocate workspace for SPBCG */

  /* vec_tmpl passed as template vector */
  spbcg_mem = NULL;
  spbcg_mem = SpbcgMalloc(maxl1, vec_tmpl);
  if (spbcg_mem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_FAIL, "KINSPILS", "KINSpbcg", MSGS_MEM_FAIL);
    free(kinspils_mem); kinspils_mem = NULL;
    return(KINSPILS_MEM_FAIL);
  }

  /* This is an iterative linear solver */

  inexact_ls = TRUE;

  /* Attach SPBCG memory to spils memory structure */
  spils_mem = (void *) spbcg_mem;

  /* attach linear solver memory to KINSOL memory */
  lmem = kinspils_mem;

  return(KINSPILS_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * additional readability replacements
 * -----------------------------------------------------------------
 */

#define maxl   (kinspils_mem->s_maxl)
#define pset   (kinspils_mem->s_pset)
#define psolve (kinspils_mem->s_psolve)
#define P_data (kinspils_mem->s_P_data)

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgInit
 * -----------------------------------------------------------------
 * This routine initializes variables associated with the SPBCG
 * iterative linear solver. Mmemory allocation was done previously
 * in KINSpbcg.
 * -----------------------------------------------------------------
 */

static int KINSpbcgInit(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;
  SpbcgMem spbcg_mem;

  kinspils_mem = (KINSpilsMem) lmem;
  spbcg_mem = (SpbcgMem) spils_mem;

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

  setupNonNull = ((psolve != NULL) && (pset != NULL));

  /* Set Jacobian-related fields, based on jtimesDQ */

  if (jtimesDQ) {
    jtimes = KINSpilsDQJtimes;
    J_data = kin_mem;
  } else {
    J_data = user_data;
  }

  if ( (strategy == KIN_PICARD) && jtimesDQ ) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL", "KINSpbcgInit", 
		    MSG_NOL_FAIL);
    return(KIN_ILL_INPUT);
  }

 /*  Set maxl in the SPBCG memory in case it was changed by the user */
  spbcg_mem->l_max  = maxl;

  last_flag = KINSPILS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the SPBCG linear
 * solver, that is, it is an interface to the user-supplied
 * routine pset.
 * -----------------------------------------------------------------
 */

static int KINSpbcgSetup(KINMem kin_mem)
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
 * Function : KINSpbcgSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic SPBCG solver routine
 * called SpbcgSolve for the solution of the linear system Ax = b.
 *
 * Appropriate variables are passed to SpbcgSolve and the counters
 * nli, nps, and ncfl are incremented, and the return value is set
 * according to the success of SpbcgSolve. The success flag is
 * returned if SpbcgSolve converged, or if the residual was reduced.
 * Of the other error conditions, only preconditioner solver
 * failure is specifically returned. Otherwise a generic flag is
 * returned to denote failure of this routine.
 * -----------------------------------------------------------------
 */

static int KINSpbcgSolve(KINMem kin_mem, N_Vector xx, N_Vector bb, 
                         realtype *sJpnorm, realtype *sFdotJp)
{
  KINSpilsMem kinspils_mem;
  SpbcgMem spbcg_mem;
  int ret, nli_inc, nps_inc;
  realtype res_norm;
  
  kinspils_mem = (KINSpilsMem) lmem;
  spbcg_mem = (SpbcgMem) spils_mem;

  /* Set initial guess to xx = 0. bb is set, by the routine
     calling KINSpbcgSolve, to the RHS vector for the system
     to be solved. */ 
 
  N_VConst(ZERO, xx);

  new_uu = TRUE;  /* set flag required for user Jacobian routine */

  /* call SpbcgSolve */

  ret = SpbcgSolve(spbcg_mem, kin_mem, xx, bb, pretype, eps,
                   kin_mem, fscale, fscale, KINSpilsAtimes,
                   KINSpilsPSolve, &res_norm, &nli_inc, &nps_inc);

  /* increment counters nli, nps, and ncfl 
     (nni is updated in the KINSol main iteration loop) */

  nli = nli + (long int) nli_inc;
  nps = nps + (long int) nps_inc;

  if (printfl > 2) 
    KINPrintInfo(kin_mem, PRNT_NLI, "KINSPBCG", "KINSpbcgSolve", INFO_NLI, nli_inc);

  if (ret != 0) ncfl++;
  last_flag = ret;

  if ( (ret != 0) && (ret != SPBCG_RES_REDUCED) ) {

    /* Handle all failure returns from SpbcgSolve */

    switch(ret) {
    case SPBCG_PSOLVE_FAIL_REC:
    case SPBCG_ATIMES_FAIL_REC:
      return(1);
      break;
    case SPBCG_CONV_FAIL:
    case SPBCG_MEM_NULL:
    case SPBCG_ATIMES_FAIL_UNREC:
    case SPBCG_PSOLVE_FAIL_UNREC:
      return(-1);
      break;
    }
  }

  /*  SpbcgSolve returned either SPBCG_SUCCESS or SPBCG_RES_REDUCED.

     Compute the terms sJpnorm and sFdotJp for use in the linesearch
     routine and in KINForcingTerm.  Both of these terms are subsequently
     corrected if the step is reduced by constraints or the linesearch.

     sJpnorm is the norm of the scaled product (scaled by fscale) of the
     current Jacobian matrix J and the step vector p (= solution vector xx).

     sFdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale.                            */

  ret = KINSpilsAtimes(kin_mem, xx, bb);
  if (ret > 0) {
    last_flag = SPBCG_ATIMES_FAIL_REC;
    return(1);
  }      
  else if (ret < 0) {
    last_flag = SPBCG_ATIMES_FAIL_UNREC;
    return(-1);
  }

  *sJpnorm = N_VWL2Norm(bb, fscale);
  N_VProd(bb, fscale, bb);
  N_VProd(bb, fscale, bb);
  *sFdotJp = N_VDotProd(fval, bb);

  if (printfl > 2) KINPrintInfo(kin_mem, PRNT_EPS, "KINSPBCG",
                     "KINSpbcgSolve", INFO_EPS, res_norm, eps);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgFree
 * -----------------------------------------------------------------
 * Frees memory specific to the SPBCG linear solver module.
 * -----------------------------------------------------------------
 */

static void KINSpbcgFree(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;
  SpbcgMem spbcg_mem;

  kinspils_mem = (KINSpilsMem) lmem;

  spbcg_mem = (SpbcgMem) spils_mem;
  SpbcgFree(spbcg_mem);

  if (kinspils_mem->s_pfree != NULL) (kinspils_mem->s_pfree)(kin_mem);

  free(kinspils_mem); kinspils_mem = NULL;
}
