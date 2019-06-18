/*
 * -----------------------------------------------------------------
 * $Revision: 4237 $
 * $Date: 2014-10-15 16:19:59 -0700 (Wed, 15 Oct 2014) $
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
 * This is the common implementation file for the IDAS Scaled              
 * Preconditioned Linear Solver modules.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "idas_spils_impl.h"
#include "idas_impl.h"

/* Private constants */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define PT05   RCONST(0.05)
#define ONE    RCONST(1.0)

/* Algorithmic constants */

#define MAX_ITERS  3  /* max. number of attempts to recover in DQ J*v */

/* 
 * =================================================================
 * PRIVATE FUNCTION PROTOTYPES
 * =================================================================
 */

static int IDAAspilsPrecSetupB(realtype tt, 
                               N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                               realtype c_jB, void *idaadj_mem,
                               N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

static int IDAAspilsPrecSetupBS(realtype tt, 
                                N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                                realtype c_jB, void *idaadj_mem,
                                N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

static int IDAAspilsPrecSolveB(realtype tt, 
                               N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                               N_Vector rvecB, N_Vector zvecB,
                               realtype c_jB, realtype deltaB,
                               void *idaadj_mem, N_Vector tmpB);

static int IDAAspilsPrecSolveBS(realtype tt, 
                                N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                                N_Vector rvecB, N_Vector zvecB,
                                realtype c_jB, realtype deltaB,
                                void *idaadj_mem, N_Vector tmpB);

static int IDAAspilsJacTimesVecB(realtype tt,
                                 N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                 N_Vector vB, N_Vector JvB, 
                                 realtype c_jB, void *idaadj_mem, 
                                 N_Vector tmp1B, N_Vector tmp2B);

static int IDAAspilsJacTimesVecBS(realtype tt,
                                  N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                  N_Vector vB, N_Vector JvB, 
                                  realtype c_jB, void *idaadj_mem, 
                                  N_Vector tmp1B, N_Vector tmp2B);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/* Readability Replacements */

#define lrw1      (IDA_mem->ida_lrw1)
#define liw1      (IDA_mem->ida_liw1)
#define tn        (IDA_mem->ida_tn)
#define cj        (IDA_mem->ida_cj)
#define res       (IDA_mem->ida_res)
#define user_data (IDA_mem->ida_user_data)
#define ewt       (IDA_mem->ida_ewt)
#define lmem      (IDA_mem->ida_lmem)

#define ils_type  (idaspils_mem->s_type)
#define sqrtN     (idaspils_mem->s_sqrtN)
#define epslin    (idaspils_mem->s_epslin)
#define ytemp     (idaspils_mem->s_ytemp)
#define yptemp    (idaspils_mem->s_yptemp)
#define xx        (idaspils_mem->s_xx)
#define ycur      (idaspils_mem->s_ycur)
#define ypcur     (idaspils_mem->s_ypcur)
#define rcur      (idaspils_mem->s_rcur)
#define npe       (idaspils_mem->s_npe)
#define nli       (idaspils_mem->s_nli)
#define nps       (idaspils_mem->s_nps)
#define ncfl      (idaspils_mem->s_ncfl)
#define njtimes   (idaspils_mem->s_njtimes)
#define nres      (idaspils_mem->s_nres)

#define jtimesDQ  (idaspils_mem->s_jtimesDQ)
#define jtimes    (idaspils_mem->s_jtimes)
#define jdata     (idaspils_mem->s_jdata)

#define last_flag (idaspils_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * OPTIONAL INPUT and OUTPUT
 * -----------------------------------------------------------------
 */

int IDASpilsSetGSType(void *ida_mem, int gstype)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetGSType", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsSetGSType", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  if (ils_type != SPILS_SPGMR) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetGSType", MSGS_BAD_LSTYPE);
    return(IDASPILS_ILL_INPUT);
  }

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetGSType", MSGS_BAD_GSTYPE);
    return(IDASPILS_ILL_INPUT);
  }

  idaspils_mem->s_gstype = gstype;

  return(IDASPILS_SUCCESS);
}

int IDASpilsSetMaxRestarts(void *ida_mem, int maxrs)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetMaxRestarts", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsSetMaxRestarts", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  if (ils_type != SPILS_SPGMR) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetMaxRestarts", MSGS_BAD_LSTYPE);
    return(IDASPILS_ILL_INPUT);
  }

  /* Check for legal maxrs */
  if (maxrs < 0) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetMaxRestarts", MSGS_NEG_MAXRS);
    return(IDASPILS_ILL_INPUT);
  }

  idaspils_mem->s_maxrs = maxrs;

  return(IDASPILS_SUCCESS);
}

int IDASpilsSetMaxl(void *ida_mem, int maxl)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetMaxl", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsSetMaxl", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  if (ils_type == SPILS_SPGMR) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetMaxl", MSGS_BAD_LSTYPE);
    return(IDASPILS_ILL_INPUT);
  }

  idaspils_mem->s_maxl = (maxl <= 0) ? IDA_SPILS_MAXL : maxl;

  return(IDASPILS_SUCCESS);
}

int IDASpilsSetEpsLin(void *ida_mem, realtype eplifac)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetEpsLin", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsSetEpsLin", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  /* Check for legal maxrs */
  if (eplifac < ZERO) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetEpsLin", MSGS_NEG_EPLIFAC);
    return(IDASPILS_ILL_INPUT);
  }

  if (eplifac == ZERO)
    idaspils_mem->s_eplifac = PT05;
  else
    idaspils_mem->s_eplifac = eplifac;

  return(IDASPILS_SUCCESS);
}

int IDASpilsSetIncrementFactor(void *ida_mem, realtype dqincfac)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetIncrementFactor", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsSetIncrementFactor", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  /* Check for legal maxrs */
  if (dqincfac <= ZERO) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetIncrementFactor", MSGS_NEG_DQINCFAC);
    return(IDASPILS_ILL_INPUT);
  }

  idaspils_mem->s_dqincfac = dqincfac;

  return(IDASPILS_SUCCESS);
}

int IDASpilsSetPreconditioner(void *ida_mem,
                              IDASpilsPrecSetupFn pset, IDASpilsPrecSolveFn psolve)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetPreconditioner", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  idaspils_mem->s_pset = pset;
  idaspils_mem->s_psolve = psolve;

  return(IDASPILS_SUCCESS);
}

int IDASpilsSetJacTimesVecFn(void *ida_mem, IDASpilsJacTimesVecFn jtv)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetJacTimesVecFn", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  if (jtv != NULL) {
    jtimesDQ = FALSE;
    jtimes = jtv;
  } else {
    jtimesDQ = TRUE;
  }

  return(IDASPILS_SUCCESS);
}

int IDASpilsGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  int maxl;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsGetWorkSpace", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  switch(ils_type) {
  case SPILS_SPGMR:
    maxl = idaspils_mem->s_maxl;
    *lenrwLS = lrw1*(maxl + 6) + maxl*(maxl + 4) + 1;
    *leniwLS = liw1*(maxl + 6);
    break;
  case SPILS_SPBCG:
    *lenrwLS = lrw1 * 10;
    *leniwLS = liw1 * 10;
    break;
  case SPILS_SPTFQMR:
    *lenrwLS = lrw1*13;
    *leniwLS = liw1*13;
    break;
  }

  return(IDASPILS_SUCCESS);
}

int IDASpilsGetNumPrecEvals(void *ida_mem, long int *npevals)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsGetNumPrecEvals", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  *npevals = npe;

  return(IDASPILS_SUCCESS);
}

int IDASpilsGetNumPrecSolves(void *ida_mem, long int *npsolves)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsGetNumPrecSolves", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  *npsolves = nps;

  return(IDASPILS_SUCCESS);
}

int IDASpilsGetNumLinIters(void *ida_mem, long int *nliters)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsGetNumLinIters", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  *nliters = nli;

  return(IDASPILS_SUCCESS);
}

int IDASpilsGetNumConvFails(void *ida_mem, long int *nlcfails)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsGetNumConvFails", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  *nlcfails = ncfl;

  return(IDASPILS_SUCCESS);
}

int IDASpilsGetNumJtimesEvals(void *ida_mem, long int *njvevals)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsGetNumJtimesEvals", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  *njvevals = njtimes;

  return(IDASPILS_SUCCESS);
}

int IDASpilsGetNumResEvals(void *ida_mem, long int *nrevalsLS)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsGetNumResEvals", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsGetNumResEvals", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  *nrevalsLS = nres;

  return(IDASPILS_SUCCESS);
}

int IDASpilsGetLastFlag(void *ida_mem, long int *flag)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsGetLastFlag", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASSPILS", "IDASpilsGetLastFlag", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) lmem;

  *flag = last_flag;

  return(IDASPILS_SUCCESS);
}

char *IDASpilsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case IDASPILS_SUCCESS:
    sprintf(name,"IDASPILS_SUCCESS");
    break; 
  case IDASPILS_MEM_NULL:
    sprintf(name,"IDASPILS_MEM_NULL");
    break;
  case IDASPILS_LMEM_NULL:
    sprintf(name,"IDASPILS_LMEM_NULL");
    break;
  case IDASPILS_ILL_INPUT:
    sprintf(name,"IDASPILS_ILL_INPUT");
    break;
  case IDASPILS_MEM_FAIL:
    sprintf(name,"IDASPILS_MEM_FAIL");
    break;
  case IDASPILS_PMEM_NULL:
    sprintf(name,"IDASPILS_PMEM_NULL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * IDASPILS private functions
 * -----------------------------------------------------------------
 */

#define psolve   (idaspils_mem->s_psolve)
#define pdata    (idaspils_mem->s_pdata)
#define dqincfac (idaspils_mem->s_dqincfac)

/*
 * This routine generates the matrix-vector product z = Jv, where
 * J is the system Jacobian, by calling either the user provided
 * routine or the internal DQ routine.
 */

int IDASpilsAtimes(void *ida_mem, N_Vector v, N_Vector z)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  int jtflag;

  IDA_mem = (IDAMem) ida_mem;
  idaspils_mem = (IDASpilsMem) lmem;

  jtflag = jtimes(tn, ycur, ypcur, rcur, v, z, cj, jdata, ytemp, yptemp);
  njtimes++;

  return(jtflag);
}

/*
 * This routine interfaces between the generic Solve routine and
 * the user's psolve routine.  It passes to psolve all required state 
 * information from ida_mem.  Its return value is the same as that
 * returned by psolve.  Note that the generic solver guarantees
 * that IDASilsPSolve will not be called in the case psolve = NULL.
 */

int IDASpilsPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  int retval;

  IDA_mem = (IDAMem) ida_mem;
  idaspils_mem = (IDASpilsMem) lmem;

  retval = psolve(tn, ycur, ypcur, rcur, r, z, cj, epslin, pdata, ytemp);

  /* This call is counted in nps within the IDASp**Solve routine */

  return(retval);

}

/*
 * This routine generates the matrix-vector product z = Jv, where
 * J is the system Jacobian, by using a difference quotient approximation.
 * The approximation is 
 *      Jv = [F(t,y1,yp1) - F(t,y,yp)]/sigma,  where
 *        y1 = y + sigma*v,  yp1 = yp + cj*sigma*v,
 *        sigma = sqrt(Neq)*dqincfac.
 * The return value from the call to res is saved in order to set the
 * return flag from IDASp**Solve.
 */

int IDASpilsDQJtimes(realtype tt,
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     N_Vector v, N_Vector Jv, 
                     realtype c_j, void *data, 
                     N_Vector work1, N_Vector work2)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  N_Vector y_tmp, yp_tmp;
  realtype sig=ZERO, siginv;
  int iter, retval;

  /* data is ida_mem */
  IDA_mem = (IDAMem) data;
  idaspils_mem = (IDASpilsMem) lmem;

  switch(ils_type) {
  case SPILS_SPGMR:
    sig = sqrtN*dqincfac;
    break;
  case SPILS_SPBCG:
    sig = dqincfac/N_VWrmsNorm(v, ewt);
    break;
  case SPILS_SPTFQMR:
    sig = dqincfac/N_VWrmsNorm(v, ewt);
    break;
  }

  /* Rename work1 and work2 for readibility */
  y_tmp  = work1;
  yp_tmp = work2;

  for (iter=0; iter<MAX_ITERS; iter++) {

    /* Set y_tmp = yy + sig*v, yp_tmp = yp + cj*sig*v. */
    N_VLinearSum(sig, v, ONE, yy, y_tmp);
    N_VLinearSum(c_j*sig, v, ONE, yp, yp_tmp);
    
    /* Call res for Jv = F(t, y_tmp, yp_tmp), and return if it failed. */
    retval = res(tt, y_tmp, yp_tmp, Jv, user_data); 
    nres++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    sig *= PT25;
  }

  if (retval > 0) return(+1);

  /* Set Jv to [Jv - rr]/sig and return. */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, rr, Jv);

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

#define yyTmp        (IDAADJ_mem->ia_yyTmp)
#define ypTmp        (IDAADJ_mem->ia_ypTmp)
#define noInterp     (IDAADJ_mem->ia_noInterp)
#define interpSensi  (IDAADJ_mem->ia_interpSensi)
#define yySTmp      (IDAADJ_mem->ia_yySTmp)
#define ypSTmp      (IDAADJ_mem->ia_ypSTmp)

/*
 * -----------------------------------------------------------------
 * OPTIONAL INPUT and OUTPUT FUNCTIONS
 * -----------------------------------------------------------------
 */

int IDASpilsSetGSTypeB(void *ida_mem, int which, int gstypeB)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetGSTypeB", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", "IDASpilsSetGSTypeB",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetGSTypeB", MSGS_BAD_WHICH);
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

  /* Call the corresponding Set* function for the backward problem. */
  return IDASpilsSetGSType(ida_memB, gstypeB);
}

int IDASpilsSetMaxRestartsB(void *ida_mem, int which, int maxrsB)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetMaxRestartsB", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", "IDASpilsSetMaxRestartsB",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetMaxRestartsB", MSGS_BAD_WHICH);
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

  /* Call the corresponding Set* function for the backward problem. */
  return IDASpilsSetMaxRestarts(ida_memB,maxrsB);
}

int IDASpilsSetEpsLinB(void *ida_mem, int which, realtype eplifacB)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetEpsLinB", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", "IDASpilsSetEpsLinB",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetEpsLinB", MSGS_BAD_WHICH);
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

  /* Call the corresponding Set* function for the backward problem. */
  return IDASpilsSetEpsLin(ida_memB, eplifacB);
}

int IDASpilsSetMaxlB(void *ida_mem, int which, int maxlB)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", "IDASpilsSetMaxlB", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", "IDASpilsSetMaxlB",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", "IDASpilsSetMaxlB", MSGS_BAD_WHICH);
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

  /* Call the corresponding Set* function for the backward problem. */
  return IDASpilsSetMaxl(ida_memB, maxlB);
}

int IDASpilsSetIncrementFactorB(void *ida_mem, int which, realtype dqincfacB)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", 
                    "IDASpilsSetIncrementFactorB", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", 
                    "IDASpilsSetIncrementFactorB",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", 
                    "IDASpilsSetIncrementFactorB", MSGS_BAD_WHICH);
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

  /* Call the corresponding Set* function for the backward problem. */
  return IDASpilsSetIncrementFactor(ida_memB, dqincfacB);
}

int IDASpilsSetPreconditionerB(void *ida_mem, int which, 
                               IDASpilsPrecSetupFnB psetB, IDASpilsPrecSolveFnB psolveB)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  IDASpilsMemB idaspilsB_mem;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", 
                    "IDASpilsSetPreconditionerB", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", 
                    "IDASpilsSetPreconditionerB",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", 
                    "IDASpilsSetPreconditionerB", MSGS_BAD_WHICH);
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

  if ( IDAB_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEMB_NULL, "IDASSPILS", 
                    "IDASpilsSetPreconditionerB", MSGS_LMEMB_NULL);
    return(IDASPILS_ILL_INPUT);
  }

  /* Get the IDASpilsMemB data. */
  idaspilsB_mem = (IDASpilsMemB) IDAB_mem->ida_lmem;

  /* Set preconditioners for the backward problem. */
  idaspilsB_mem->s_psetB   = psetB;
  idaspilsB_mem->s_psolveB = psolveB;

  if (psetB == NULL) {
    return IDASpilsSetPreconditioner(ida_memB, NULL, IDAAspilsPrecSolveB);
  } else {
    return IDASpilsSetPreconditioner(ida_memB, IDAAspilsPrecSetupB, IDAAspilsPrecSolveB);
  }
}

int IDASpilsSetPreconditionerBS(void *ida_mem, int which, 
                                IDASpilsPrecSetupFnBS psetBS,
                                IDASpilsPrecSolveFnBS psolveBS)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  IDASpilsMemB idaspilsB_mem;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", 
                    "IDASpilsSetPreconditionerBS", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", 
                    "IDASpilsSetPreconditionerBS",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", 
                    "IDASpilsSetPreconditionerBS", MSGS_BAD_WHICH);
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

  if ( IDAB_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEMB_NULL, "IDASSPILS", 
                    "IDASpilsSetPreconditionerBS", MSGS_LMEMB_NULL);
    return(IDASPILS_ILL_INPUT);
  }

  /* Get the IDASpilsMemB data. */
  idaspilsB_mem = (IDASpilsMemB) IDAB_mem->ida_lmem;

  /* Set preconditioners for the backward problem. */
  idaspilsB_mem->s_psetBS   = psetBS;
  idaspilsB_mem->s_psolveBS = psolveBS;

  if (psetBS == NULL) {
    return IDASpilsSetPreconditioner(ida_memB, NULL, IDAAspilsPrecSolveBS);
  } else {
    return IDASpilsSetPreconditioner(ida_memB, IDAAspilsPrecSetupBS, IDAAspilsPrecSolveBS);
  }
}


int IDASpilsSetJacTimesVecFnB(void *ida_mem, int which, IDASpilsJacTimesVecFnB jtvB)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  IDASpilsMemB idaspilsB_mem;
  void *ida_memB;
  int flag;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", 
                    "IDASpilsSetJacTimesVecFnB", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", 
                    "IDASpilsSetJacTimesVecFnB",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", 
                    "IDASpilsSetJacTimesVecFnB", MSGS_BAD_WHICH);
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

  if ( IDAB_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEMB_NULL, "IDASSPILS", 
                    "IDASpilsSetJacTimesVecFnB", MSGS_LMEMB_NULL);
    return(IDASPILS_ILL_INPUT);
  }

  idaspilsB_mem = (IDASpilsMemB) IDAB_mem->ida_lmem;

  /* Call the corresponding Set* function for the backward problem. */

  idaspilsB_mem->s_jtimesB   = jtvB;

  if (jtvB != NULL) {
    flag = IDASpilsSetJacTimesVecFn(ida_memB, IDAAspilsJacTimesVecB);
  } else {
    flag = IDASpilsSetJacTimesVecFn(ida_memB, NULL);
  }
  return(flag);
}

int IDASpilsSetJacTimesVecFnBS(void *ida_mem, int which, IDASpilsJacTimesVecFnBS jtvBS)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  IDASpilsMemB idaspilsB_mem;
  void *ida_memB;
  int flag;
  
  /* Check if ida_mem is allright. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS", 
                    "IDASpilsSetJacTimesVecFnBS", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS", 
                    "IDASpilsSetJacTimesVecFnBS",  MSGS_NO_ADJ);
    return(IDASPILS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS", 
                    "IDASpilsSetJacTimesVecFnBS", MSGS_BAD_WHICH);
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

  if ( IDAB_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEMB_NULL, "IDASSPILS", 
                    "IDASpilsSetJacTimesVecFnBS", MSGS_LMEMB_NULL);
    return(IDASPILS_ILL_INPUT);
  }

  idaspilsB_mem = (IDASpilsMemB) IDAB_mem->ida_lmem;

  /* Call the corresponding Set* function for the backward problem. */

  idaspilsB_mem->s_jtimesBS = jtvBS;

  if (jtvBS != NULL) {
    flag = IDASpilsSetJacTimesVecFn(ida_memB, IDAAspilsJacTimesVecBS);
  } else {
    flag = IDASpilsSetJacTimesVecFn(ida_memB, NULL);
  }
  return(flag);
}


/*
 * -----------------------------------------------------------------
 * IDASSPILS private functions
 *
 * These routines are interfaces to the routines provided by the
 * user for backward problems.
 *
 * -----------------------------------------------------------------
 */

static int IDAAspilsPrecSetupB(realtype tt, 
                              N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                              realtype c_jB, void *ida_mem,
                              N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDASpilsMemB idaspilsB_mem; 
  IDABMem IDAB_mem;
  int flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
  idaspilsB_mem = (IDASpilsMemB)IDAB_mem->ida_lmem;
  
  /* Get forward solution from interpolation. */
  if (noInterp==FALSE) {
    flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, NULL, NULL);
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDAB_mem->IDA_mem, -1, "IDASSPILS", "IDAAspilsPrecSetupB", MSGS_BAD_T);
      return(-1);
    }
  }

  /* Call user's adjoint precondB routine */
  flag = idaspilsB_mem->s_psetB(tt, yyTmp, ypTmp, yyB, ypB, rrB, 
                                c_jB, IDAB_mem->ida_user_data,
                                tmp1B, tmp2B, tmp3B);
  return(flag);
}

static int IDAAspilsPrecSetupBS(realtype tt, 
                               N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                               realtype c_jB, void *ida_mem,
                               N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDASpilsMemB idaspilsB_mem;
  IDABMem IDAB_mem;
  int flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
  idaspilsB_mem = (IDASpilsMemB)IDAB_mem->ida_lmem;
  
  /* Get forward solution from interpolation. */
  if( noInterp == FALSE) {
    if (interpSensi)
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, yySTmp, ypSTmp);
    else
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, NULL, NULL);
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDAB_mem->IDA_mem, -1, "IDASSPILS", "IDAAspilsPrecSetupBS", MSGS_BAD_T);
      return(-1);
    }
  }

  /* Call user's adjoint precondBS routine */
  flag = idaspilsB_mem->s_psetBS(tt, yyTmp, ypTmp, yySTmp, ypSTmp,
                                 yyB, ypB, rrB,
                                 c_jB, IDAB_mem->ida_user_data,
                                 tmp1B, tmp2B, tmp3B);
  return(flag);
}


static int IDAAspilsPrecSolveB(realtype tt, 
                               N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                               N_Vector rvecB, N_Vector zvecB,
                               realtype c_jB, realtype deltaB,
                               void *ida_mem, N_Vector tmpB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDASpilsMemB idaspilsB_mem; 
  IDABMem IDAB_mem;
  int flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
  idaspilsB_mem = (IDASpilsMemB)IDAB_mem->ida_lmem;
  
  /* Get forward solution from interpolation. */
  if (noInterp==FALSE) {
    flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, NULL, NULL);
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDAB_mem->IDA_mem, -1, "IDASSPILS", "IDAAspilsPrecSolveB", MSGS_BAD_T);
      return(-1);
    }
  }
  /* Call user's adjoint psolveB routine */
  flag = idaspilsB_mem->s_psolveB(tt, yyTmp, ypTmp, 
                                  yyB, ypB, rrB, 
                                  rvecB, zvecB, 
                                  c_jB, deltaB, 
                                  IDAB_mem->ida_user_data, tmpB);
  return(flag);
}

static int IDAAspilsPrecSolveBS(realtype tt, 
                                N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                                N_Vector rvecB, N_Vector zvecB,
                                realtype c_jB, realtype deltaB,
                                void *ida_mem, N_Vector tmpB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDASpilsMemB idaspilsB_mem;
  IDABMem IDAB_mem;
  int flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
  idaspilsB_mem = (IDASpilsMemB)IDAB_mem->ida_lmem;
  
  /* Get forward solution from interpolation. */
  if( noInterp == FALSE) {
    if (interpSensi)
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, yySTmp, ypSTmp);
    else
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, NULL, NULL);
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDAB_mem->IDA_mem, -1, "IDASSPILS", "IDAAspilsPrecSolveBS", MSGS_BAD_T);
      return(-1);
    }
  }

  /* Call user's adjoint psolveBS routine */
  flag = idaspilsB_mem->s_psolveBS(tt, yyTmp, ypTmp, yySTmp, ypSTmp,
                                   yyB, ypB, rrB,
                                   rvecB, zvecB,
                                   c_jB, deltaB,
                                   IDAB_mem->ida_user_data, tmpB);
  return(flag);
}


static int IDAAspilsJacTimesVecB(realtype tt,
                                 N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                 N_Vector vB, N_Vector JvB, 
                                 realtype c_jB, void *ida_mem, 
                                 N_Vector tmp1B, N_Vector tmp2B)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDASpilsMemB idaspilsB_mem; 
  IDABMem IDAB_mem;
  int flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
  idaspilsB_mem = (IDASpilsMemB)IDAB_mem->ida_lmem;
  
  /* Get forward solution from interpolation. */
  if (noInterp==FALSE) {
    flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, NULL, NULL);
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDAB_mem->IDA_mem, -1, "IDASSPILS", "IDAAspilsJacTimesVecB", MSGS_BAD_T);
      return(-1);
    }
  }
  /* Call user's adjoint jtimesB routine */
  flag = idaspilsB_mem->s_jtimesB(tt, yyTmp, ypTmp, 
                                  yyB, ypB, rrB, 
                                  vB, JvB, 
                                  c_jB, IDAB_mem->ida_user_data, 
                                  tmp1B, tmp2B);
  return(flag);
}

static int IDAAspilsJacTimesVecBS(realtype tt,
                                  N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                  N_Vector vB, N_Vector JvB, 
                                  realtype c_jB, void *ida_mem, 
                                  N_Vector tmp1B, N_Vector tmp2B)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDASpilsMemB idaspilsB_mem;
  IDABMem IDAB_mem;
  int flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
  idaspilsB_mem = (IDASpilsMemB)IDAB_mem->ida_lmem;
  
  /* Get forward solution from interpolation. */
  if( noInterp == FALSE) {
    if (interpSensi)
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, yySTmp, ypSTmp);
    else
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, NULL, NULL);
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDAB_mem->IDA_mem, -1, "IDASSPILS", "IDAAspilsJacTimesVecBS", MSGS_BAD_T);
      return(-1);
    }
  }

  /* Call user's adjoint jtimesBS routine */
  flag = idaspilsB_mem->s_jtimesBS(tt, yyTmp, ypTmp, yySTmp, ypSTmp,
                                   yyB, ypB, rrB,
                                   vB, JvB,
                                   c_jB, IDAB_mem->ida_user_data,
                                   tmp1B, tmp2B);
  return(flag);
}
