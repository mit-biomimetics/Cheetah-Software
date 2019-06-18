/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
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
 * This file contains implementations of the banded difference
 * quotient Jacobian-based preconditioner and solver routines for
 * use with the CVSPILS linear solvers.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_impl.h"
#include "cvodes_bandpre_impl.h"
#include "cvodes_spils_impl.h"

#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_spgmr.h>

#include <sundials/sundials_math.h>

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Prototypes of cvBandPrecSetup and cvBandPrecSolve */
  
static int cvBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int cvBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bp_data, N_Vector tmp);

/* Prototype for cvBandPrecFree */

static void cvBandPrecFree(CVodeMem cv_mem);

/* Prototype for difference quotient Jacobian calculation routine */

static int cvBandPrecDQJac(CVBandPrecData pdata, 
                           realtype t, N_Vector y, N_Vector fy, 
                           N_Vector ftemp, N_Vector ytemp);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/* Redability replacements */

#define vec_tmpl (cv_mem->cv_tempv)

/*
 * -----------------------------------------------------------------
 * Initialization, Free, and Get Functions
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVBandPrecInit will
 *       first test for a compatible N_Vector internal representation
 *       by checking that the function N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */
int CVBandPrecInit(void *cvode_mem, long int N, long int mu, long int ml)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBandPrecData pdata;
  long int mup, mlp, storagemu;
  int flag;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBANDPRE", "CVBandPrecInit", MSGBP_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Test if the NVECTOR package is compatible with the BAND preconditioner */
  if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVBANDPRE", "CVBandPrecInit", MSGBP_BAD_NVECTOR);
    return(CVSPILS_ILL_INPUT);
  }

  pdata = NULL;
  pdata = (CVBandPrecData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Load pointers and bandwidths into pdata block. */
  pdata->cvode_mem = cvode_mem;
  pdata->N = N;
  pdata->mu = mup = SUNMIN(N-1, SUNMAX(0,mu));
  pdata->ml = mlp = SUNMIN(N-1, SUNMAX(0,ml));

  /* Initialize nfeBP counter */
  pdata->nfeBP = 0;

  /* Allocate memory for saved banded Jacobian approximation. */
  pdata->savedJ = NULL;
  pdata->savedJ = NewBandMat(N, mup, mlp, mup);
  if (pdata->savedJ == NULL) {
    free(pdata); pdata = NULL;
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }


  /* Allocate memory for banded preconditioner. */
  storagemu = SUNMIN(N-1, mup+mlp);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Allocate memory for pivot array. */
  pdata->lpivots = NULL;
  pdata->lpivots = NewLintArray(N);
  if (pdata->lpivots == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBANDPRE", "CVBandPrecInit", MSGBP_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Overwrite the P_data field in the SPILS memory */
  cvspils_mem->s_P_data = pdata;

  /* Attach the pfree function */
  cvspils_mem->s_pfree = cvBandPrecFree;

  /* Attach preconditioner solve and setup functions */
  flag = CVSpilsSetPreconditioner(cvode_mem, cvBandPrecSetup, cvBandPrecSolve);

  return(flag);
}

int CVBandPrecGetWorkSpace(void *cvode_mem, long int *lenrwBP, long int *leniwBP)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBandPrecData pdata;
  long int N, ml, mu, smu;

  
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBANDPRE", "CVBandPrecGetWorkSpace", MSGBP_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBANDPRE", "CVBandPrecGetWorkSpace", MSGBP_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  if (cvspils_mem->s_P_data == NULL) {
    cvProcessError(cv_mem, CVSPILS_PMEM_NULL, "CVBANDPRE", "CVBandPrecGetWorkSpace", MSGBP_PMEM_NULL);
    return(CVSPILS_PMEM_NULL);
  } 
  pdata = (CVBandPrecData) cvspils_mem->s_P_data;

  N   = pdata->N;
  mu  = pdata->mu;
  ml  = pdata->ml;
  smu = SUNMIN( N-1, mu + ml);

  *leniwBP = pdata->N;
  *lenrwBP = N * ( 2*ml + smu + mu + 2 );

  return(CVSPILS_SUCCESS);
}

int CVBandPrecGetNumRhsEvals(void *cvode_mem, long int *nfevalsBP)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBandPrecData pdata;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBANDPRE", "CVBandPrecGetNumRhsEvals", MSGBP_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBANDPRE", "CVBandPrecGetNumRhsEvals", MSGBP_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  if (cvspils_mem->s_P_data == NULL) {
    cvProcessError(cv_mem, CVSPILS_PMEM_NULL, "CVBANDPRE", "CVBandPrecGetNumRhsEvals", MSGBP_PMEM_NULL);
    return(CVSPILS_PMEM_NULL);
  } 
  pdata = (CVBandPrecData) cvspils_mem->s_P_data;

  *nfevalsBP = pdata->nfeBP;

  return(CVSPILS_SUCCESS);
}

/* Readability Replacements */

#define N       (pdata->N)
#define mu      (pdata->mu)
#define ml      (pdata->ml)
#define lpivots (pdata->lpivots)
#define savedJ  (pdata->savedJ)
#define savedP  (pdata->savedP)
#define nfeBP   (pdata->nfeBP)

/*
 * -----------------------------------------------------------------
 * cvBandPrecSetup
 * -----------------------------------------------------------------
 * Together cvBandPrecSetup and cvBandPrecSolve use a banded
 * difference quotient Jacobian to create a preconditioner.
 * cvBandPrecSetup calculates a new J, if necessary, then
 * calculates P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of cvBandPrecSetup are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *         namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == FALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == TRUE means that Jacobian data from the
 *                  previous PrecSetup call will be reused
 *                  (with the current value of gamma).
 *         A cvBandPrecSetup call with jok == TRUE should only
 *         occur after a call with jok == FALSE.
 *
 * *jcurPtr is a pointer to an output integer flag which is
 *          set by CVBandPrecond as follows:
 *            *jcurPtr = TRUE if Jacobian data was recomputed.
 *            *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                       but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * bp_data is a pointer to preconditoner data (set by CVBandPrecInit)
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *           for vectors of length N for work space. This
 *           routine uses only tmp1 and tmp2.
 *
 * The value to be returned by the cvBandPrecSetup function is
 *   0  if successful, or
 *   1  if the band factorization failed.
 * -----------------------------------------------------------------
 */

static int cvBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CVBandPrecData pdata;
  CVodeMem cv_mem;
  int retval;
  long int ier;

  /* Assume matrix and lpivots have already been allocated. */
  pdata = (CVBandPrecData) bp_data;

  cv_mem = (CVodeMem) pdata->cvode_mem;

  if (jok) {

    /* If jok = TRUE, use saved copy of J. */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mu, ml);

  } else {

    /* If jok = FALSE, call cvBandPrecDQJac for new J value. */
    *jcurPtr = TRUE;
    SetToZero(savedJ);

    retval = cvBandPrecDQJac(pdata, t, y, fy, tmp1, tmp2);
    if (retval < 0) {
      cvProcessError(cv_mem, -1, "CVBANDPRE", "cvBandPrecSetup", MSGBP_RHSFUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(savedJ, savedP, mu, ml);

  }
  
  /* Scale and add I to get savedP = I - gamma*J. */
  BandScale(-gamma, savedP);
  AddIdentity(savedP);
 
  /* Do LU factorization of matrix. */
  ier = BandGBTRF(savedP, lpivots);
 
  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvBandPrecSolve
 * -----------------------------------------------------------------
 * cvBandPrecSolve solves a linear system P z = r, where P is the
 * matrix computed by CVBandPrecond.
 *
 * The parameters of cvBandPrecSolve used here are as follows:
 *
 * r       is the right-hand side vector of the linear system.
 *
 * bp_data is a pointer to preconditoner data (set by CVBandPrecInit)
 *
 * z       is the output vector computed by cvBandPrecSolve.
 *
 * The value returned by the cvBandPrecSolve function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */ 

static int cvBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bp_data, N_Vector tmp)
{
  CVBandPrecData pdata;
  realtype *zd;

  /* Assume matrix and lpivots have already been allocated. */
  pdata = (CVBandPrecData) bp_data;

  /* Copy r to z. */
  N_VScale(ONE, r, z);

  /* Do band backsolve on the vector z. */
  zd = N_VGetArrayPointer(z);

  BandGBTRS(savedP, lpivots, zd);

  return(0);
}


static void cvBandPrecFree(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;
  CVBandPrecData pdata;

  if (cv_mem->cv_lmem == NULL) return;
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  
  if (cvspils_mem->s_P_data == NULL) return;
  pdata = (CVBandPrecData) cvspils_mem->s_P_data;

  DestroyMat(savedJ);
  DestroyMat(savedP);
  DestroyArray(lpivots);

  free(pdata);
  pdata = NULL;
}


#define ewt       (cv_mem->cv_ewt)
#define uround    (cv_mem->cv_uround)
#define h         (cv_mem->cv_h)
#define f         (cv_mem->cv_f)
#define user_data (cv_mem->cv_user_data)

/*
 * -----------------------------------------------------------------
 * cvBandPrecDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static int cvBandPrecDQJac(CVBandPrecData pdata, 
                           realtype t, N_Vector y, N_Vector fy, 
                           N_Vector ftemp, N_Vector ytemp)
{
  CVodeMem cv_mem;
  realtype fnorm, minInc, inc, inc_inv, srur;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int retval;

  cv_mem = (CVodeMem) pdata->cvode_mem;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp. */
  ewt_data   = N_VGetArrayPointer(ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector. */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f. */
  srur = SUNRsqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing. */
  width = ml + mu + 1;
  ngroups = SUNMIN(width, N);
  
  for (group = 1; group <= ngroups; group++) {
    
    /* Increment all y_j in group. */
    for(j = group-1; j < N; j += width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y. */

    retval = f(t, ytemp, ftemp, user_data);
    nfeBP++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients. */
    for (j = group-1; j < N; j += width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(savedJ,j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mu);
      i2 = SUNMIN(j+ml, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

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
 * CVBandPrecInitB, CVBPSp*B
 *
 * Wrappers for the backward phase around the corresponding 
 * CVODES functions
 */

int CVBandPrecInitB(void *cvode_mem, int which, long int nB, long int muB, long int mlB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBANDPRE", "CVBandPrecInitB", MSGBP_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVBANDPRE", "CVBandPrecInitB", MSGBP_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVBANDPRE", "CVBandPrecInitB", MSGBP_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvB_mem->cv_pfree = NULL;

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVBandPrecInit(cvodeB_mem, nB, muB, mlB);

  return(flag);
}

