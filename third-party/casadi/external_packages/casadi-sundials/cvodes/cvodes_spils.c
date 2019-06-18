/*
 * -----------------------------------------------------------------
 * $Revision: 4233 $
 * $Date: 2014-10-08 16:11:32 -0700 (Wed, 08 Oct 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s):Radu Serban @ LLNL
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
 * This is the implementation file for the CVSPILS linear solvers.
 *
 * Part II contains wrappers for using the CVODES iterative linear 
 * solvers on adjoint (backward) problems.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_impl.h"
#include "cvodes_spils_impl.h"

/* Private constants */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define ONE    RCONST(1.0)

/* Algorithmic constants */

#define MAX_ITERS  3  /* max. number of attempts to recover in DQ J*v */

/* 
 * =================================================================
 * PRIVATE FUNCTION PROTOTYPES
 * =================================================================
 */

/*
 * cvSpilsPrecSetupBWrapper and cvSpilsPrecSetupBSWrapper have type
 * CVSpilsPrecSetupFn, and wrap around user-provided functions of
 * type CVSpilsPrecSetupFnB and CVSpilsPrecSetupFnBS, respectively.
 */

static int cvSpilsPrecSetupBWrapper(realtype t, N_Vector yB,
                                    N_Vector fyB, booleantype jokB,
                                    booleantype *jcurPtrB, realtype gammaB,
                                    void *cvode_mem,
                                    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

static int cvSpilsPrecSetupBSWrapper(realtype t, N_Vector yB,
                                     N_Vector fyB, booleantype jokB,
                                     booleantype *jcurPtrB, realtype gammaB,
                                     void *cvode_mem,
                                     N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/*
 * cvSpilsPrecSolveBWrapper and cvSpilsPrecSolveBSWrapper have type
 * CVSpilsPrecSolveFn, and wrap around user-provided functions of
 * type CVSpilsPrecSolveFnB and CVSpilsPrecSolveFnBS, respectively.
 */

static int cvSpilsPrecSolveBWrapper(realtype t, N_Vector yB, N_Vector fyB,
                                    N_Vector rB, N_Vector zB,
                                    realtype gammaB, realtype deltaB,
                                    int lrB, void *cvode_mem, N_Vector tmpB);
  
static int cvSpilsPrecSolveBSWrapper(realtype t, N_Vector yB, N_Vector fyB,
                                     N_Vector rB, N_Vector zB,
                                     realtype gammaB, realtype deltaB,
                                     int lrB, void *cvode_mem, N_Vector tmpB);
  
/*
 * cvSpilsJacTimesVecBWrapper and cvSpilsJacTimesVecBSWrapper have type
 * CVSpilsJacTimesVecFn, and wrap around user-provided functions of
 * type CVSpilsJacTimesVecFnB and CVSpilsJacTimesVecFnBS, respectively.
 */

static int cvSpilsJacTimesVecBWrapper(N_Vector vB, N_Vector JvB, realtype t, 
                                      N_Vector yB, N_Vector fyB, 
                                      void *cvode_mem, N_Vector tmpB);

static int cvSpilsJacTimesVecBSWrapper(N_Vector vB, N_Vector JvB, realtype t, 
                                       N_Vector yB, N_Vector fyB, 
                                       void *cvode_mem, N_Vector tmpB);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/* Readability Replacements */

#define lrw1      (cv_mem->cv_lrw1)
#define liw1      (cv_mem->cv_liw1)
#define tq        (cv_mem->cv_tq)
#define tn        (cv_mem->cv_tn)
#define h         (cv_mem->cv_h)
#define gamma     (cv_mem->cv_gamma)
#define nfe       (cv_mem->cv_nfe)
#define f         (cv_mem->cv_f)
#define user_data (cv_mem->cv_user_data)
#define ewt       (cv_mem->cv_ewt)
#define lmem      (cv_mem->cv_lmem)

#define ils_type (cvspils_mem->s_type)
#define sqrtN    (cvspils_mem->s_sqrtN)   
#define ytemp    (cvspils_mem->s_ytemp)
#define x        (cvspils_mem->s_x)
#define ycur     (cvspils_mem->s_ycur)
#define fcur     (cvspils_mem->s_fcur)
#define delta    (cvspils_mem->s_delta)
#define npe      (cvspils_mem->s_npe)
#define nli      (cvspils_mem->s_nli)
#define nps      (cvspils_mem->s_nps)
#define ncfl     (cvspils_mem->s_ncfl)
#define njtimes  (cvspils_mem->s_njtimes)
#define nfes     (cvspils_mem->s_nfes)

#define jtimesDQ (cvspils_mem->s_jtimesDQ)
#define jtimes   (cvspils_mem->s_jtimes)
#define j_data   (cvspils_mem->s_j_data)

#define last_flag (cvspils_mem->s_last_flag)


/*
 * -----------------------------------------------------------------
 * OPTIONAL INPUT and OUTPUT FUNCTIONS
 * -----------------------------------------------------------------
 */


/*
 * -----------------------------------------------------------------
 * CVSpilsSetPrecType
 * -----------------------------------------------------------------
 */

int CVSpilsSetPrecType(void *cvode_mem, int pretype)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetPrecType", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetPrecType", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetPrecType", MSGS_BAD_PRETYPE);
    return(CVSPILS_ILL_INPUT);
  }

  cvspils_mem->s_pretype = pretype;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsSetGSType
 * -----------------------------------------------------------------
 */

int CVSpilsSetGSType(void *cvode_mem, int gstype)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetGSType", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetGSType", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  if (ils_type != SPILS_SPGMR) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetGSType", MSGS_BAD_LSTYPE);
    return(CVSPILS_ILL_INPUT);
  }

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetGSType", MSGS_BAD_GSTYPE);
    return(CVSPILS_ILL_INPUT);
  }

  cvspils_mem->s_gstype = gstype;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpilsSetMaxl
 * -----------------------------------------------------------------
 */

int CVSpilsSetMaxl(void *cvode_mem, int maxl)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetMaxl", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(NULL, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetMaxl", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  if (ils_type == SPILS_SPGMR) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetMaxl", MSGS_BAD_LSTYPE);
    return(CVSPILS_ILL_INPUT);
  }

  mxl = (maxl <= 0) ? CVSPILS_MAXL : maxl;
  cvspils_mem->s_maxl = mxl;

  /*  spbcg_mem->l_max  = mxl; */

  return(CVSPILS_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * CVSpilsSetEpsLin
 * -----------------------------------------------------------------
 */

int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetEpsLin", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetEpsLin", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  /* Check for legal eplifac */
  if(eplifac < ZERO) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetEpsLin", MSGS_BAD_EPLIN);
    return(CVSPILS_ILL_INPUT);
  }

  cvspils_mem->s_eplifac = (eplifac == ZERO) ? CVSPILS_EPLIN : eplifac;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsSetPrecSetupFn
 * -----------------------------------------------------------------
 */

int CVSpilsSetPreconditioner(void *cvode_mem,
                             CVSpilsPrecSetupFn pset, CVSpilsPrecSolveFn psolve)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetPreconditioner", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  cvspils_mem->s_pset = pset;
  cvspils_mem->s_psolve = psolve;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsSetJacTimesVecFn
 * -----------------------------------------------------------------
 */

int CVSpilsSetJacTimesVecFn(void *cvode_mem, CVSpilsJacTimesVecFn jtv)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFn", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  if (jtv != NULL) {
    jtimesDQ = FALSE;
    jtimes = jtv;
  } else {
    jtimesDQ = TRUE;
  }

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetWorkSpace
 * -----------------------------------------------------------------
 */

int CVSpilsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  int maxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsGetWorkSpace", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  
  switch(ils_type) {
  case SPILS_SPGMR:
    maxl = cvspils_mem->s_maxl;
    *lenrwLS = lrw1*(maxl + 5) + maxl*(maxl + 4) + 1;
    *leniwLS = liw1*(maxl + 5);
    break;
  case SPILS_SPBCG:
    *lenrwLS = lrw1 * 9;
    *leniwLS = liw1 * 9;
    break;
  case SPILS_SPTFQMR:
    *lenrwLS = lrw1*11;
    *leniwLS = liw1*11;
    break;
  }


  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetNumPrecEvals
 * -----------------------------------------------------------------
 */

int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsGetNumPrecEvals", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  *npevals = npe;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetNumPrecSolves
 * -----------------------------------------------------------------
 */

int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsGetNumPrecSolves", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  *npsolves = nps;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetNumLinIters
 * -----------------------------------------------------------------
 */

int CVSpilsGetNumLinIters(void *cvode_mem, long int *nliters)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsGetNumLinIters", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  *nliters = nli;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetNumConvFails
 * -----------------------------------------------------------------
 */

int CVSpilsGetNumConvFails(void *cvode_mem, long int *nlcfails)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsGetNumConvFails", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  *nlcfails = ncfl;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetNumJtimesEvals
 * -----------------------------------------------------------------
 */

int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsGetNumJtimesEvals", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  *njvevals = njtimes;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetNumRhsEvals
 * -----------------------------------------------------------------
 */

int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsGetNumRhsEvals", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsGetNumRhsEvals", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  *nfevalsLS = nfes;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetLastFlag
 * -----------------------------------------------------------------
 */

int CVSpilsGetLastFlag(void *cvode_mem, long int *flag)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsGetLastFlag", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsGetLastFlag", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) lmem;

  *flag = last_flag;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *CVSpilsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVSPILS_SUCCESS:
    sprintf(name,"CVSPILS_SUCCESS");
    break;    
  case CVSPILS_MEM_NULL:
    sprintf(name,"CVSPILS_MEM_NULL");
    break;
  case CVSPILS_LMEM_NULL:
    sprintf(name,"CVSPILS_LMEM_NULL");
    break;
  case CVSPILS_ILL_INPUT:
    sprintf(name,"CVSPILS_ILL_INPUT");
    break;
  case CVSPILS_MEM_FAIL:
    sprintf(name,"CVSPILS_MEM_FAIL");
    break;
  case CVSPILS_PMEM_NULL:
    sprintf(name,"CVSPILS_PMEM_NULL");
    break;
  case CVSPILS_NO_ADJ:
    sprintf(name,"CVSPILS_NO_ADJ");
    break;
  case CVSPILS_LMEMB_NULL:
    sprintf(name,"CVSPILS_LMEMB_NULL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * CVSPILS private functions
 * -----------------------------------------------------------------
 */


/* Additional readability Replacements */

#define pretype (cvspils_mem->s_pretype)
#define eplifac (cvspils_mem->s_eplifac)
#define maxl    (cvspils_mem->s_maxl)
#define psolve  (cvspils_mem->s_psolve)
#define P_data  (cvspils_mem->s_P_data)

/*
 * -----------------------------------------------------------------
 * CVSpilsAtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Mv, where
 * M = I - gamma*J. The product J*v is obtained by calling the jtimes 
 * routine. It is then scaled by -gamma and added to v to obtain M*v.
 * The return value is the same as the value returned by jtimes --
 * 0 if successful, nonzero otherwise.
 * -----------------------------------------------------------------
 */

int CVSpilsAtimes(void *cvode_mem, N_Vector v, N_Vector z)
{
  CVodeMem   cv_mem;
  CVSpilsMem cvspils_mem;
  int retval;

  cv_mem = (CVodeMem) cvode_mem;
  cvspils_mem = (CVSpilsMem) lmem;

  retval = jtimes(v, z, tn, ycur, fcur, j_data, ytemp);
  njtimes++;
  if (retval != 0) return(retval);

  N_VLinearSum(ONE, v, -gamma, z, z);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * CVSpilsPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic Sp***Solve routine
 * (within the SPGMR, SPBCG, or SPTFQMR solver) and the
 * user's psolve routine.  It passes to psolve all required state 
 * information from cvode_mem.  Its return value is the same as that
 * returned by psolve. Note that the generic SP*** solver guarantees
 * that CVSpilsPSolve will not be called in the case in which
 * preconditioning is not done. This is the only case in which the
 * user's psolve routine is allowed to be NULL.
 * -----------------------------------------------------------------
 */

int CVSpilsPSolve(void *cvode_mem, N_Vector r, N_Vector z, int lr)
{
  CVodeMem   cv_mem;
  CVSpilsMem cvspils_mem;
  int retval;

  cv_mem = (CVodeMem) cvode_mem;
  cvspils_mem = (CVSpilsMem)lmem;

  /* This call is counted in nps within the CVSp***Solve routine */
  retval = psolve(tn, ycur, fcur, r, z, gamma, delta, lr, P_data, ytemp);

  return(retval);     
}

/*
 * -----------------------------------------------------------------
 * CVSpilsDQJtimes
 * -----------------------------------------------------------------
 * This routine generates a difference quotient approximation to
 * the Jacobian times vector f_y(t,y) * v. The approximation is 
 * Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) is
 * input, i.e. the WRMS norm of v/vnrm is 1.
 * -----------------------------------------------------------------
 */

int CVSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t, 
                    N_Vector y, N_Vector fy,
                    void *data, N_Vector work)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  realtype sig, siginv;
  int iter, retval;

  /* data is cvode_mem */
  cv_mem = (CVodeMem) data;
  cvspils_mem = (CVSpilsMem) lmem;

  /* Initialize perturbation to 1/||v|| */
  sig = ONE/N_VWrmsNorm(v, ewt);

  for (iter=0; iter<MAX_ITERS; iter++) {

    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = f(t, work, Jv, user_data); 
    nfes++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    sig *= PT25;
  }

  if (retval > 0) return(+1);

  /* Replace Jv by (Jv - fy)/sig */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return(0);
}

/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/* Readability replacements */

#define ytmp  (ca_mem->ca_ytmp)
#define yStmp (ca_mem->ca_yStmp)
#define IMget (ca_mem->ca_IMget)
#define IMinterpSensi (ca_mem->ca_IMinterpSensi)

#define pset_B     (cvspilsB_mem->s_psetB)
#define pset_BS     (cvspilsB_mem->s_psetBS)
#define psolve_B   (cvspilsB_mem->s_psolveB)
#define psolve_BS   (cvspilsB_mem->s_psolveBS)
#define jtimes_B   (cvspilsB_mem->s_jtimesB)
#define jtimes_BS   (cvspilsB_mem->s_jtimesBS)

/*
 * -----------------------------------------------------------------
 * OPTIONAL INPUT and OUTPUT FUNCTIONS
 * -----------------------------------------------------------------
 */

int CVSpilsSetPrecTypeB(void *cvode_mem, int which, int pretypeB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetPrecTypeB", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPILS", "CVSpilsSetPrecTypeB", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetPrecTypeB", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVSpilsSetPrecType(cvodeB_mem, pretypeB);

  return(flag);
}

int CVSpilsSetGSTypeB(void *cvode_mem, int which, int gstypeB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetGSTypeB", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPILS", "CVSpilsSetGSTypeB", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetGSTypeB", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVSpilsSetGSType(cvodeB_mem,gstypeB);

  return(flag);
}

int CVSpilsSetEpsLinB(void *cvode_mem, int which, realtype eplifacB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetEpsLinB", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPILS", "CVSpilsSetEpsLinB", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetEpsLinB", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVSpilsSetEpsLin(cvodeB_mem,eplifacB);

  return(flag);
}

int CVSpilsSetMaxlB(void *cvode_mem, int which, int maxlB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetMaxlB", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPILS", "CVSpilsSetMaxlB", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetMaxlB", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVSpilsSetMaxl(cvodeB_mem,maxlB);

  return(flag);
}

int CVSpilsSetPreconditionerB(void *cvode_mem, int which, 
                              CVSpilsPrecSetupFnB psetB,
                              CVSpilsPrecSolveFnB psolveB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem; 
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetPreconditionerB", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPILS", "CVSpilsSetPreconditionerB", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetPreconditionerB", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  if (cvB_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEMB_NULL, "CVSPILS", "CVSpilsSetPreconditionerB", MSGS_LMEMB_NULL);
    return(CVSPILS_LMEMB_NULL);
  }
  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  pset_B   = psetB;
  psolve_B = psolveB;

  if (psetB == NULL) {
    flag = CVSpilsSetPreconditioner(cvodeB_mem, NULL, cvSpilsPrecSolveBWrapper);
  } else {
    flag = CVSpilsSetPreconditioner(cvodeB_mem, cvSpilsPrecSetupBWrapper, cvSpilsPrecSolveBWrapper);
  }

  return(flag);
}

int CVSpilsSetPreconditionerBS(void *cvode_mem, int which, 
                               CVSpilsPrecSetupFnBS psetBS,
                               CVSpilsPrecSolveFnBS psolveBS)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem; 
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetPreconditionerBS", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPILS", "CVSpilsSetPreconditionerBS", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetPreconditionerBS", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  if (cvB_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEMB_NULL, "CVSPILS", "CVSpilsSetPreconditionerBS", MSGS_LMEMB_NULL);
    return(CVSPILS_LMEMB_NULL);
  }
  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  pset_BS   = psetBS;
  psolve_BS = psolveBS;

  if (psetBS == NULL) {
    flag = CVSpilsSetPreconditioner(cvodeB_mem, NULL, cvSpilsPrecSolveBSWrapper);
  } else {
    flag = CVSpilsSetPreconditioner(cvodeB_mem, cvSpilsPrecSetupBSWrapper, cvSpilsPrecSolveBSWrapper);
  }

  return(flag);
}

int CVSpilsSetJacTimesVecFnB(void *cvode_mem, int which, CVSpilsJacTimesVecFnB jtvB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem; 
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFnB", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPILS", "CVSpilsSetJacTimesVecFnB", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetJacTimesVecFnB", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  if (cvB_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEMB_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFnB", MSGS_LMEMB_NULL);
    return(CVSPILS_LMEMB_NULL);
  }
  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  jtimes_B = jtvB;

  if (jtvB != NULL) {
    flag = CVSpilsSetJacTimesVecFn(cvodeB_mem, cvSpilsJacTimesVecBWrapper);
  } else {
    flag = CVSpilsSetJacTimesVecFn(cvodeB_mem, NULL);
  }

  return(flag);
}

int CVSpilsSetJacTimesVecFnBS(void *cvode_mem, int which, CVSpilsJacTimesVecFnBS jtvBS)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem; 
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFnBS", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPILS", "CVSpilsSetJacTimesVecFnBS", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", "CVSpilsSetJacTimesVecFnBS", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  if (cvB_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEMB_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFnBS", MSGS_LMEMB_NULL);
    return(CVSPILS_LMEMB_NULL);
  }
  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  jtimes_BS = jtvBS;

  if (jtvBS != NULL) {
    flag = CVSpilsSetJacTimesVecFn(cvodeB_mem, cvSpilsJacTimesVecBSWrapper);
  } else {
    flag = CVSpilsSetJacTimesVecFn(cvodeB_mem, NULL);
  }

  return(flag);
}


/*
 * -----------------------------------------------------------------
 * CVSPILS private functions
 * -----------------------------------------------------------------
 */

/*
 * cvSpilsPrecSetupBWrapper
 *
 * This routine interfaces to the CVSpilsPrecSetupFnB routine 
 * provided by the user.
 */

static int cvSpilsPrecSetupBWrapper(realtype t, N_Vector yB, 
                                    N_Vector fyB, booleantype jokB, 
                                    booleantype *jcurPtrB, realtype gammaB,
                                    void *cvode_mem,
                                    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  flag = IMget(cv_mem, t, ytmp, NULL);
  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVSPILS", "cvSpilsPrecSetupBWrapper", MSGS_BAD_TINTERP);
    return(-1);
  } 

  /* Call user's adjoint precondB routine */
  retval = pset_B(t, ytmp, yB, fyB, jokB, jcurPtrB, gammaB,
                  cvB_mem->cv_user_data, tmp1B, tmp2B, tmp3B);

  return(retval);
}

/*
 * cvSpilsPrecSetupBSWrapper
 *
 * This routine interfaces to the CVSpilsPrecSetupFnBS routine 
 * provided by the user.
 */

static int cvSpilsPrecSetupBSWrapper(realtype t, N_Vector yB, 
                                     N_Vector fyB, booleantype jokB, 
                                     booleantype *jcurPtrB, realtype gammaB,
                                     void *cvode_mem,
                                     N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  if (IMinterpSensi)
    flag = IMget(cv_mem, t, ytmp, yStmp);
  else
    flag = IMget(cv_mem, t, ytmp, NULL);
  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVSPILS", "cvSpilsPrecSetupBSWrapper", MSGS_BAD_TINTERP);
    return(-1);
  } 

  /* Call user's adjoint precondB routine */
  retval = pset_BS(t, ytmp, yStmp, yB, fyB, jokB, jcurPtrB, gammaB,
                   cvB_mem->cv_user_data, tmp1B, tmp2B, tmp3B);

  return(retval);
}

/*
 * cvSpilsPrecSolveBWrapper
 *
 * This routine interfaces to the CVSpilsPrecSolveFnB routine 
 * provided by the user.
 */

static int cvSpilsPrecSolveBWrapper(realtype t, N_Vector yB, N_Vector fyB,
                                    N_Vector rB, N_Vector zB,
                                    realtype gammaB, realtype deltaB,
                                    int lrB, void *cvode_mem, N_Vector tmpB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  flag = IMget(cv_mem, t, ytmp, NULL);
  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVSPILS", "cvSpilsPrecSolveBWrapper", MSGS_BAD_TINTERP);
    return(-1);
  }

  /* Call user's adjoint psolveB routine */
  retval = psolve_B(t, ytmp, yB, fyB, rB, zB, gammaB, deltaB, 
                    lrB, cvB_mem->cv_user_data, tmpB);

  return(retval);
}

/*
 * cvSpilsPrecSolveBSWrapper
 *
 * This routine interfaces to the CVSpilsPrecSolveFnBS routine 
 * provided by the user.
 */

static int cvSpilsPrecSolveBSWrapper(realtype t, N_Vector yB, N_Vector fyB,
                                     N_Vector rB, N_Vector zB,
                                     realtype gammaB, realtype deltaB,
                                     int lrB, void *cvode_mem, N_Vector tmpB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  if (IMinterpSensi)
    flag = IMget(cv_mem, t, ytmp, yStmp);
  else
    flag = IMget(cv_mem, t, ytmp, NULL);
  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVSPILS", "cvSpilsPrecSolveBSWrapper", MSGS_BAD_TINTERP);
    return(-1);
  }

  /* Call user's adjoint psolveBS routine */
  retval = psolve_BS(t, ytmp, yStmp, yB, fyB, rB, zB, gammaB, deltaB, 
                     lrB, cvB_mem->cv_user_data, tmpB);

  return(retval);
}

/*
 * cvSpilsJacTimesVecBWrapper
 *
 * This routine interfaces to the CVSpilsJacTimesVecFnB routine 
 * provided by the user.
 */

static int cvSpilsJacTimesVecBWrapper(N_Vector vB, N_Vector JvB, realtype t, 
                                      N_Vector yB, N_Vector fyB, 
                                      void *cvode_mem, N_Vector tmpB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  flag = IMget(cv_mem, t, ytmp, NULL);
  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVSPILS", "cvSpilsJacTimesVecBWrapper", MSGS_BAD_TINTERP);
    return(-1);
  } 

  /* Call user's adjoint jtimesB routine */
  retval = jtimes_B(vB, JvB, t, ytmp, yB, fyB, cvB_mem->cv_user_data, tmpB);

  return(retval);
}

/*
 * cvSpilsJacTimesVecBSWrapper
 *
 * This routine interfaces to the CVSpilsJacTimesVecFnBS routine 
 * provided by the user.
 */

static int cvSpilsJacTimesVecBSWrapper(N_Vector vB, N_Vector JvB, realtype t, 
                                       N_Vector yB, N_Vector fyB, 
                                       void *cvode_mem, N_Vector tmpB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSpilsMemB cvspilsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  if (IMinterpSensi)
    flag = IMget(cv_mem, t, ytmp, yStmp);
  else
    flag = IMget(cv_mem, t, ytmp, NULL);
  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVSPILS", "cvSpilsJacTimesVecBSWrapper", MSGS_BAD_TINTERP);
    return(-1);
  } 

  /* Call user's adjoint jtimesBS routine */
  retval = jtimes_BS(vB, JvB, t, ytmp, yStmp, yB, fyB, cvB_mem->cv_user_data, tmpB);

  return(retval);
}
