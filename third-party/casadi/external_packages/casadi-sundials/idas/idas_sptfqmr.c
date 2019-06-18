/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
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
 * This is the implementation file for the IDAS scaled preconditioned
 * TFQMR linear solver module, IDASPTFQMR.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <idas/idas_sptfqmr.h>
#include "idas_spils_impl.h"
#include "idas_impl.h"

#include <sundials/sundials_sptfqmr.h>
#include <sundials/sundials_math.h>

/* Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define PT9  RCONST(0.9)
#define PT05 RCONST(0.05)

/* IDASPTFQMR linit, lsetup, lsolve, lperf, and lfree routines */

static int IDASptfqmrInit(IDAMem IDA_mem);

static int IDASptfqmrSetup(IDAMem IDA_mem, 
			   N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
			   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int IDASptfqmrSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
			   N_Vector yy_now, N_Vector yp_now, N_Vector rr_now);

static int IDASptfqmrPerf(IDAMem IDA_mem, int perftask);

static int IDASptfqmrFree(IDAMem IDA_mem);

/* IDASPTFQMR lfree function for backward problem. */
static void IDASptfqmrFreeB(IDABMem IDAB_mem);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/* Readability Replacements */

#define nst          (IDA_mem->ida_nst)
#define tn           (IDA_mem->ida_tn)
#define cj           (IDA_mem->ida_cj)
#define epsNewt      (IDA_mem->ida_epsNewt)
#define res          (IDA_mem->ida_res)
#define user_data    (IDA_mem->ida_user_data)
#define ewt          (IDA_mem->ida_ewt)
#define errfp        (IDA_mem->ida_errfp)
#define linit        (IDA_mem->ida_linit)
#define lsetup       (IDA_mem->ida_lsetup)
#define lsolve       (IDA_mem->ida_lsolve)
#define lperf        (IDA_mem->ida_lperf)
#define lfree        (IDA_mem->ida_lfree)
#define lmem         (IDA_mem->ida_lmem)
#define nni          (IDA_mem->ida_nni)
#define ncfn         (IDA_mem->ida_ncfn)
#define setupNonNull (IDA_mem->ida_setupNonNull)
#define vec_tmpl     (IDA_mem->ida_tempv1)

#define sqrtN       (idaspils_mem->s_sqrtN)
#define epslin      (idaspils_mem->s_epslin)
#define ytemp       (idaspils_mem->s_ytemp)
#define yptemp      (idaspils_mem->s_yptemp)
#define xx          (idaspils_mem->s_xx)
#define ycur        (idaspils_mem->s_ycur)
#define ypcur       (idaspils_mem->s_ypcur)
#define rcur        (idaspils_mem->s_rcur)
#define npe         (idaspils_mem->s_npe)
#define nli         (idaspils_mem->s_nli)
#define nps         (idaspils_mem->s_nps)
#define ncfl        (idaspils_mem->s_ncfl)
#define nst0        (idaspils_mem->s_nst0)
#define nni0        (idaspils_mem->s_nni0)
#define nli0        (idaspils_mem->s_nli0)
#define ncfn0       (idaspils_mem->s_ncfn0)
#define ncfl0       (idaspils_mem->s_ncfl0)
#define nwarn       (idaspils_mem->s_nwarn)
#define njtimes     (idaspils_mem->s_njtimes)
#define nres        (idaspils_mem->s_nres)
#define spils_mem   (idaspils_mem->s_spils_mem)

#define jtimesDQ    (idaspils_mem->s_jtimesDQ)
#define jtimes      (idaspils_mem->s_jtimes)
#define jdata       (idaspils_mem->s_jdata)

#define last_flag   (idaspils_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmr
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDASPTFQMR linear solver module.
 *
 * IDASptfqmr first calls the existing lfree routine if this is not NULL.
 * It then sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDASptfqmrInit, IDASptfqmrSetup,
 * IDASptfqmrSolve, IDASptfqmrPerf, and IDASptfqmrFree, respectively.
 * It allocates memory for a structure of type IDASpilsMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem). It then sets various fields
 * in the IDASpilsMemRec structure. Finally, IDASptfqmr allocates 
 * memory for ytemp, yptemp, and xx, and calls SptfqmrMalloc to 
 * allocate memory for the Sptfqmr solver.
 *
 * The return value of IDASptfqmr is:
 *   IDASPILS_SUCCESS   =  0 if successful
 *   IDASPILS_MEM_FAIL  = -1 if IDA_mem is NULL or a memory
 *                             allocation failed
 *   IDASPILS_ILL_INPUT = -2 if a required vector operation is not
 *                             implemented.
 * -----------------------------------------------------------------
 */

int IDASptfqmr(void *ida_mem, int maxl)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  SptfqmrMem sptfqmr_mem;
  int flag, maxl1;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPTFQMR", "IDASptfqmr", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if N_VDotProd is present */
  if (vec_tmpl->ops->nvdotprod == NULL) {
    IDAProcessError(NULL, IDASPILS_ILL_INPUT, "IDASPTFQMR", "IDASptfqmr", MSGS_BAD_NVECTOR);
    return(IDASPILS_ILL_INPUT);
  }

  if (lfree != NULL) flag = lfree((IDAMem) ida_mem);

  /* Set five main function fields in ida_mem */
  linit  = IDASptfqmrInit;
  lsetup = IDASptfqmrSetup;
  lsolve = IDASptfqmrSolve;
  lperf  = IDASptfqmrPerf;
  lfree  = IDASptfqmrFree;

  /* Get memory for IDASpilsMemRec */
  idaspils_mem = NULL;
  idaspils_mem = (IDASpilsMem) malloc(sizeof(struct IDASpilsMemRec));
  if (idaspils_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPTFQMR", "IDASptfqmr", MSGS_MEM_FAIL);
    return(IDASPILS_MEM_FAIL);
  }

  /* Set ILS type */
  idaspils_mem->s_type = SPILS_SPTFQMR;

  /* Set SPTFQMR parameters that were passed in call sequence */
  maxl1 = (maxl <= 0) ? IDA_SPILS_MAXL : maxl;
  idaspils_mem->s_maxl = maxl1;

  /* Set defaults for Jacobian-related fileds */
  jtimesDQ = TRUE;
  jtimes   = NULL;
  jdata    = NULL;

  /* Set defaults for preconditioner-related fields */
  idaspils_mem->s_pset   = NULL;
  idaspils_mem->s_psolve = NULL;
  idaspils_mem->s_pfree  = NULL;
  idaspils_mem->s_pdata  = IDA_mem->ida_user_data;

  /* Set default values for the rest of the Sptfqmr parameters */
  idaspils_mem->s_eplifac   = PT05;
  idaspils_mem->s_dqincfac  = ONE;

  idaspils_mem->s_last_flag = IDASPILS_SUCCESS;

  /* Set setupNonNull to FALSE */
  setupNonNull = FALSE;

  /* Allocate memory for ytemp, yptemp, and xx */

  ytemp = N_VClone(vec_tmpl);
  if (ytemp == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPTFQMR", "IDASptfqmr", MSGS_MEM_FAIL);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  yptemp = N_VClone(vec_tmpl);
  if (yptemp == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPTFQMR", "IDASptfqmr", MSGS_MEM_FAIL);
    N_VDestroy(ytemp);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  xx = N_VClone(vec_tmpl);
  if (xx == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPTFQMR", "IDASptfqmr", MSGS_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = SUNRsqrt(N_VDotProd(ytemp, ytemp));

  /* Call SptfqmrMalloc to allocate workspace for Sptfqmr */
  sptfqmr_mem = NULL;
  sptfqmr_mem = SptfqmrMalloc(maxl1, vec_tmpl);
  if (sptfqmr_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPTFQMR", "IDASptfqmr", MSGS_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    N_VDestroy(xx);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  /* Attach SPTFQMR memory to spils memory structure */
  spils_mem = (void *)sptfqmr_mem;

  /* Attach linear solver memory to the integrator memory */
  lmem = idaspils_mem;

  return(IDASPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDASPTFQMR interface routines
 * -----------------------------------------------------------------
 */

/* Additional readability Replacements */

#define maxl     (idaspils_mem->s_maxl)
#define eplifac  (idaspils_mem->s_eplifac)
#define psolve   (idaspils_mem->s_psolve)
#define pset     (idaspils_mem->s_pset)
#define pdata    (idaspils_mem->s_pdata)

static int IDASptfqmrInit(IDAMem IDA_mem)
{
  IDASpilsMem idaspils_mem;
  SptfqmrMem sptfqmr_mem;

  idaspils_mem = (IDASpilsMem) lmem;
  sptfqmr_mem = (SptfqmrMem) spils_mem;

  /* Initialize counters */
  npe = nli = nps = ncfl = 0;
  njtimes = nres = 0;

  /* Set setupNonNull to TRUE iff there is preconditioning with setup */
  setupNonNull = (psolve != NULL) && (pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (jtimesDQ) {
    jtimes = IDASpilsDQJtimes;
    jdata = IDA_mem;
  } else {
    jdata = user_data;
  }

  /*  Set maxl in the SPTFQMR memory in case it was changed by the user */
  sptfqmr_mem->l_max  = maxl;

  last_flag = IDASPILS_SUCCESS;

  return(0);
}

static int IDASptfqmrSetup(IDAMem IDA_mem, 
			   N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
			   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int retval;
  IDASpilsMem idaspils_mem;

  idaspils_mem = (IDASpilsMem) lmem;

  /* Call user setup routine pset and update counter npe */
  retval = pset(tn, yy_p, yp_p, rr_p, cj, pdata,
                tmp1, tmp2, tmp3);
  npe++;

  if (retval < 0) {
    IDAProcessError(IDA_mem, SPTFQMR_PSET_FAIL_UNREC, "IDASPTFQMR", "IDASptfqmrSetup", MSGS_PSET_FAILED);
    last_flag = SPTFQMR_PSET_FAIL_UNREC;
    return(-1);
  }
  if (retval > 0) {
    last_flag = SPTFQMR_PSET_FAIL_REC;
    return(+1);
  }

  last_flag = SPTFQMR_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmrSolve
 * -----------------------------------------------------------------
 * Note: The x-scaling and b-scaling arrays are both equal to weight.
 *
 * We set the initial guess, x = 0, then call SptfqmrSolve.
 * We copy the solution x into b, and update the counters nli, nps,
 * and ncfl. If SptfqmrSolve returned nli_inc = 0 (hence x = 0), we
 * take the SPTFQMR vtemp vector (= P_inverse F) as the correction
 * vector instead. Finally, we set the return value according to the
 * success of SptfqmrSolve.
 * -----------------------------------------------------------------
 */

static int IDASptfqmrSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
			   N_Vector yy_now, N_Vector yp_now, N_Vector rr_now)
{
  IDASpilsMem idaspils_mem;
  SptfqmrMem sptfqmr_mem;
  int pretype, nli_inc, nps_inc, retval;
  realtype res_norm;

  idaspils_mem = (IDASpilsMem) lmem;

  sptfqmr_mem = (SptfqmrMem)spils_mem;

  /* Set SptfqmrSolve convergence test constant epslin, in terms of the
     Newton convergence test constant epsNewt and safety factors. The factor
     sqrt(Neq) assures that the TFQMR convergence test is applied to the
     WRMS norm of the residual vector, rather than the weighted L2 norm. */
  epslin = sqrtN*eplifac*epsNewt;

  /* Set vectors ycur, ypcur, and rcur for use by the Atimes and Psolve */
  ycur = yy_now;
  ypcur = yp_now;
  rcur = rr_now;

  /* Set SptfqmrSolve inputs pretype and initial guess xx = 0 */  
  pretype = (psolve == NULL) ? PREC_NONE : PREC_LEFT;
  N_VConst(ZERO, xx);
  
  /* Call SptfqmrSolve and copy xx to bb */
  retval = SptfqmrSolve(sptfqmr_mem, IDA_mem, xx, bb, pretype, epslin,
                      IDA_mem, weight, weight, IDASpilsAtimes,
                      IDASpilsPSolve, &res_norm, &nli_inc, &nps_inc);

  if (nli_inc == 0) N_VScale(ONE, SPTFQMR_VTEMP(sptfqmr_mem), bb);
  else N_VScale(ONE, xx, bb);
  
  /* Increment counters nli, nps, and return if successful */
  nli += nli_inc;
  nps += nps_inc;
  if (retval != SPTFQMR_SUCCESS) ncfl++;

    /* Interpret return value from SpgmrSolve */

  last_flag = retval;

  switch(retval) {

  case SPTFQMR_SUCCESS:
    return(0);
    break;
  case SPTFQMR_RES_REDUCED:
    return(1);
    break;
  case SPTFQMR_CONV_FAIL:
    return(1);
    break;
  case SPTFQMR_PSOLVE_FAIL_REC:
    return(1);
    break;
  case SPTFQMR_ATIMES_FAIL_REC:
    return(1);
    break;
  case SPTFQMR_MEM_NULL:
    return(-1);
    break;
  case SPTFQMR_ATIMES_FAIL_UNREC:
    IDAProcessError(IDA_mem, SPTFQMR_ATIMES_FAIL_UNREC, "IDASPTFQMR", "IDASptfqmrSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SPTFQMR_PSOLVE_FAIL_UNREC:
    IDAProcessError(IDA_mem, SPTFQMR_PSOLVE_FAIL_UNREC, "IDASPTFQMR", "IDASptfqmrSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmrPerf
 * -----------------------------------------------------------------
 * This routine handles performance monitoring specific to the
 * IDASPTFQMR linear solver. When perftask = 0, it saves values of
 * various counters. When perftask = 1, it examines difference
 * quotients in these counters, and depending on their values, it
 * prints up to three warning messages. Messages are printed up to
 * a maximum of 10 times.
 * -----------------------------------------------------------------
 */

static int IDASptfqmrPerf(IDAMem IDA_mem, int perftask)
{
  IDASpilsMem idaspils_mem;
  realtype avdim, rcfn, rcfl;
  long int nstd, nnid;
  booleantype lavd, lcfn, lcfl;

  idaspils_mem = (IDASpilsMem) lmem;

  if (perftask == 0) {
    nst0 = nst;  nni0 = nni;  nli0 = nli;
    ncfn0 = ncfn;  ncfl0 = ncfl;  
    nwarn = 0;
    return(0);
  }

  nstd = nst - nst0;  nnid = nni - nni0;
  if (nstd == 0 || nnid == 0) return(0);
  avdim = (realtype) ((nli - nli0)/((realtype) nnid));
  rcfn = (realtype) ((ncfn - ncfn0)/((realtype) nstd));
  rcfl = (realtype) ((ncfl - ncfl0)/((realtype) nnid));
  lavd = (avdim > ((realtype) maxl));
  lcfn = (rcfn > PT9);
  lcfl = (rcfl > PT9);
  if (!(lavd || lcfn || lcfl)) return(0);
  nwarn++;
  if (nwarn > 10) return(1);
  if (lavd) 
    IDAProcessError(IDA_mem, IDA_WARNING, "IDASPTFQMR", "IDASptfqmrPerf", MSGS_AVD_WARN, tn, avdim);
  if (lcfn) 
    IDAProcessError(IDA_mem, IDA_WARNING, "IDASPTFQMR", "IDASptfqmrPerf", MSGS_CFN_WARN, tn, rcfn);
  if (lcfl) 
    IDAProcessError(IDA_mem, IDA_WARNING, "IDASPTFQMR", "IDASptfqmrPerf", MSGS_CFL_WARN, tn, rcfl);

  return(0);
}

static int IDASptfqmrFree(IDAMem IDA_mem)
{
  IDASpilsMem idaspils_mem;
  SptfqmrMem sptfqmr_mem;

  idaspils_mem = (IDASpilsMem) lmem;

  N_VDestroy(ytemp);
  N_VDestroy(yptemp);
  N_VDestroy(xx);
  
  sptfqmr_mem = (SptfqmrMem)spils_mem;
  SptfqmrFree(sptfqmr_mem);

  if (idaspils_mem->s_pfree != NULL) (idaspils_mem->s_pfree)(IDA_mem);

  free(idaspils_mem); idaspils_mem = NULL;

  return(0);
}

/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */
/*
 * IDASptfqmrB
 *
 * Wrapper for the backward phase
 *
 */

int IDASptfqmrB(void *ida_mem, int which, int maxlB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  IDASpilsMemB idaspilsB_mem;
  void *ida_memB;
  int flag;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPTFQMR", "IDASptfqmrB", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASPTFQMR", "IDASptfqmrB",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASPTFQMR", "IDASptfqmrB", MSGS_BAD_WHICH);
    return(IDASPILS_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  /* ida_mem corresponding to 'which' problem. */
  ida_memB = (void *) IDAB_mem->IDA_mem;

  
  /* Get memory for IDASpilsMemRecB */
  idaspilsB_mem = NULL;
  idaspilsB_mem = (IDASpilsMemB) malloc(sizeof(struct IDASpilsMemRecB));
  if (idaspilsB_mem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_MEM_FAIL, "IDASPTFQMR", "IDASptfqmrB", MSGS_MEM_FAIL);
    return(IDASPILS_MEM_FAIL);
  }

  idaspilsB_mem->s_psetB = NULL;
  idaspilsB_mem->s_psolveB = NULL;
  idaspilsB_mem->s_P_dataB = NULL;

  /* initialize Jacobian function */
  idaspilsB_mem->s_jtimesB = NULL;

  /* attach lmem and lfree */
  IDAB_mem->ida_lmem = idaspilsB_mem;
  IDAB_mem->ida_lfree = IDASptfqmrFreeB;

  flag = IDASptfqmr(IDAB_mem->IDA_mem, maxlB);

  if (flag != IDASPILS_SUCCESS) {
    free(idaspilsB_mem);
    idaspilsB_mem = NULL;
  }

  return(flag);
}

/*
 * IDASptfqmrFreeB 
 */

static void IDASptfqmrFreeB(IDABMem IDAB_mem)
{
  IDASpilsMemB idaspilsB_mem;

  idaspilsB_mem = (IDASpilsMemB) IDAB_mem->ida_lmem;

  free(idaspilsB_mem);
}
