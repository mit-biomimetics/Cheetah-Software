/*
 * -----------------------------------------------------------------
 * $Revision: 4368 $
 * $Date: 2015-02-12 12:25:15 -0800 (Thu, 12 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
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
 * This is the implementation file for the optional input and output
 * functions for the KINSOL solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "kinsol_impl.h"
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#define ZERO      RCONST(0.0)
#define POINT1    RCONST(0.1)
#define ONETHIRD  RCONST(0.3333333333333333)
#define HALF      RCONST(0.5)
#define TWOTHIRDS RCONST(0.6666666666666667)
#define POINT9    RCONST(0.9)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)
#define TWOPT5    RCONST(2.5)

#define liw  (kin_mem->kin_liw)
#define lrw  (kin_mem->kin_lrw)
#define liw1 (kin_mem->kin_liw1)
#define lrw1 (kin_mem->kin_lrw1)

/* 
 * =================================================================
 * KINSOL optional input functions
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * KINSetErrHandlerFn
 * -----------------------------------------------------------------
 */

int KINSetErrHandlerFn(void *kinmem, KINErrHandlerFn ehfun, void *eh_data)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetErrHandlerFn", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  kin_mem->kin_ehfun = ehfun;
  kin_mem->kin_eh_data = eh_data;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetErrFile
 * -----------------------------------------------------------------
 */

int KINSetErrFile(void *kinmem, FILE *errfp)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetErrFile", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_errfp = errfp;

  return(KIN_SUCCESS);
}

#define errfp (kin_mem->kin_errfp)

/*
 * -----------------------------------------------------------------
 * Function : KINSetPrintLevel
 * -----------------------------------------------------------------
 */

int KINSetPrintLevel(void *kinmem, int printfl)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetPrintLevel", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if ((printfl < 0) || (printfl > 3)) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetPrintLevel", MSG_BAD_PRINTFL);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_printfl = printfl;

  return(KIN_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * KINSetInfoHandlerFn
 * -----------------------------------------------------------------
 */

int KINSetInfoHandlerFn(void *kinmem, KINInfoHandlerFn ihfun, void *ih_data)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetInfoHandlerFn", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  kin_mem->kin_ihfun = ihfun;
  kin_mem->kin_ih_data = ih_data;

  return(KIN_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * Function : KINSetInfoFile
 * -----------------------------------------------------------------
 */

int KINSetInfoFile(void *kinmem, FILE *infofp)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetInfoFile", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_infofp = infofp;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetUserData
 * -----------------------------------------------------------------
 */

int KINSetUserData(void *kinmem, void *user_data)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetUserData", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_user_data = user_data;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMAA
 * -----------------------------------------------------------------
 */

int KINSetMAA(void *kinmem, long int maa)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetMAA", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (maa < 0) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetMAA", MSG_BAD_MAA);
    return(KIN_ILL_INPUT);
  }

  if (maa > kin_mem->kin_mxiter) maa = kin_mem->kin_mxiter;

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_m_aa = maa;
  kin_mem->kin_aamem_aa = (maa == 0) ? FALSE : TRUE;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetAAStopCrit
 * -----------------------------------------------------------------
 */

/*  CSW: This function is currently not supported.

int KINSetAAStopCrit(void *kinmem, booleantype setstop)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetAAStopCrit", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_setstop_aa = setstop;

  return(KIN_SUCCESS);
}
*/

/*
 * -----------------------------------------------------------------
 * Function : KINSetNumMaxIters
 * -----------------------------------------------------------------
 */

int KINSetNumMaxIters(void *kinmem, long int mxiter)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetNumMaxIters", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (mxiter < 0) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetNumMaxIters", MSG_BAD_MXITER);
    return(KIN_ILL_INPUT);
  }

  if (mxiter == 0)
    kin_mem->kin_mxiter = MXITER_DEFAULT;
  else
    kin_mem->kin_mxiter = mxiter;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoInitSetup
 * -----------------------------------------------------------------
 */

int KINSetNoInitSetup(void *kinmem, booleantype noInitSetup)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetNoInitSetup", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_noInitSetup = noInitSetup;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoResMon
 * -----------------------------------------------------------------
 */

int KINSetNoResMon(void *kinmem, booleantype noResMon)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetNoResMon", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_noResMon = noResMon;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxSetupCalls
 * -----------------------------------------------------------------
 */

int KINSetMaxSetupCalls(void *kinmem, long int msbset)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetMaxSetupCalls", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (msbset < 0) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetMaxSetupCalls", MSG_BAD_MSBSET);
    return(KIN_ILL_INPUT);
  }
  
  if (msbset == 0)
    kin_mem->kin_msbset = MSBSET_DEFAULT;
  else
    kin_mem->kin_msbset = msbset;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxSubSetupCalls
 * -----------------------------------------------------------------
 */

int KINSetMaxSubSetupCalls(void *kinmem, long int msbsetsub)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetMaxSubSetupCalls", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (msbsetsub < 0) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetMaxSubSetupCalls", MSG_BAD_MSBSETSUB);
    return(KIN_ILL_INPUT);
  }
  
  if (msbsetsub == 0)
    kin_mem->kin_msbset_sub = MSBSET_SUB_DEFAULT;
  else
    kin_mem->kin_msbset_sub = msbsetsub;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaForm
 * -----------------------------------------------------------------
 */

int KINSetEtaForm(void *kinmem, int etachoice)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetEtaForm", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if ((etachoice != KIN_ETACONSTANT) && 
      (etachoice != KIN_ETACHOICE1)  && 
      (etachoice != KIN_ETACHOICE2)) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetEtaForm", MSG_BAD_ETACHOICE);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_etaflag = etachoice;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaConstValue
 * -----------------------------------------------------------------
 */

int KINSetEtaConstValue(void *kinmem, realtype eta)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetEtaConstValue", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if ((eta < ZERO) || (eta > ONE)) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetEtaConstValue", MSG_BAD_ETACONST);
    return(KIN_ILL_INPUT);
  }

  if (eta == ZERO)
    kin_mem->kin_eta = POINT1;
  else
    kin_mem->kin_eta = eta;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaParams
 * -----------------------------------------------------------------
 */

int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetEtaParams", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if ((ealpha <= ONE) || (ealpha > TWO))
    if (ealpha != ZERO) {
      KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetEtaParams", MSG_BAD_ALPHA);
      return(KIN_ILL_INPUT);
    }
  
  if (ealpha == ZERO) 
    kin_mem->kin_eta_alpha = TWO;
  else
    kin_mem->kin_eta_alpha = ealpha;

  if ((egamma <= ZERO) || (egamma > ONE))
    if (egamma != ZERO) {
      KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetEtaParams", MSG_BAD_GAMMA);
      return(KIN_ILL_INPUT);
    }

  if (egamma == ZERO)
    kin_mem->kin_eta_gamma = POINT9;
  else
    kin_mem->kin_eta_gamma = egamma;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetResMonParams
 * -----------------------------------------------------------------
 */

int KINSetResMonParams(void *kinmem, realtype omegamin, realtype omegamax)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetResMonParams", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  /* check omegamin */

  if (omegamin < ZERO) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetResMonParams", MSG_BAD_OMEGA);
    return(KIN_ILL_INPUT);
  }

  if (omegamin == ZERO) 
    kin_mem->kin_omega_min = OMEGA_MIN;
  else
    kin_mem->kin_omega_min = omegamin;

  /* check omegamax */

  if (omegamax < ZERO) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetResMonParams", MSG_BAD_OMEGA);
    return(KIN_ILL_INPUT);
  }

  if (omegamax == ZERO) {

    if (kin_mem->kin_omega_min > OMEGA_MAX) {
      KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetResMonParams", MSG_BAD_OMEGA);
      return(KIN_ILL_INPUT);
    }
    else kin_mem->kin_omega_max = OMEGA_MAX;

  } else {

    if (kin_mem->kin_omega_min > omegamax) {
      KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetResMonParams", MSG_BAD_OMEGA);
      return(KIN_ILL_INPUT);
    }
    else kin_mem->kin_omega_max = omegamax;

  }

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetResMonConstValue
 * -----------------------------------------------------------------
 */

int KINSetResMonConstValue(void *kinmem, realtype omegaconst)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetResMonConstValue", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  /* check omegaconst */

  if (omegaconst < ZERO) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetResMonConstValue", MSG_BAD_OMEGA);
    return(KIN_ILL_INPUT);
  }

  /* Load omega value. A value of 0 will force using omega_min and omega_max */
  kin_mem->kin_omega = omegaconst;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoMinEps
 * -----------------------------------------------------------------
 */

int KINSetNoMinEps(void *kinmem, booleantype noMinEps)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetNoMinEps", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_noMinEps = noMinEps;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxNewtonStep
 * -----------------------------------------------------------------
 */

int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetMaxNewtonStep", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (mxnewtstep < ZERO) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetMaxNewtonStep", MSG_BAD_MXNEWTSTEP);
    return(KIN_ILL_INPUT);
  }

  /* Note: passing a value of 0.0 will use the default
     value (computed in KINSolInit) */

  kin_mem->kin_mxnstepin = mxnewtstep;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxBetaFails
 * -----------------------------------------------------------------
 */

int KINSetMaxBetaFails(void *kinmem, long int mxnbcf)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetMaxBetaFails", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (mxnbcf < 0) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetMaxBetaFails", MSG_BAD_MXNBCF);
    return(KIN_ILL_INPUT);
  }

  if (mxnbcf == 0)
    kin_mem->kin_mxnbcf = MXNBCF_DEFAULT;
  else
    kin_mem->kin_mxnbcf = mxnbcf;

  return(KIN_SUCCESS);

}

/*
 * -----------------------------------------------------------------
 * Function : KINSetRelErrFunc
 * -----------------------------------------------------------------
 */

int KINSetRelErrFunc(void *kinmem, realtype relfunc)
{
  KINMem kin_mem;
  realtype uround;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetRelErrFunc", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (relfunc < ZERO) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetRelErrFunc", MSG_BAD_RELFUNC);
    return(KIN_ILL_INPUT);
  }

  if (relfunc == ZERO) {
    uround = kin_mem->kin_uround;
    kin_mem->kin_sqrt_relfunc = SUNRsqrt(uround);
  } else {
    kin_mem->kin_sqrt_relfunc = SUNRsqrt(relfunc);
  }

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetFuncNormTol
 * -----------------------------------------------------------------
 */

int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
{
  KINMem kin_mem;
  realtype uround;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetFuncNormTol", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (fnormtol < ZERO) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetFuncNormTol", MSG_BAD_FNORMTOL);
    return(KIN_ILL_INPUT);
  }

  if (fnormtol == ZERO) {
    uround = kin_mem->kin_uround;
    kin_mem->kin_fnormtol = SUNRpowerR(uround,ONETHIRD);
  } else {
    kin_mem->kin_fnormtol = fnormtol;
  }

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetScaledStepTol
 * -----------------------------------------------------------------
 */

int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
{
  KINMem kin_mem;
  realtype uround;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetScaledStepTol", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (scsteptol < ZERO) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetScaledStepTol", MSG_BAD_SCSTEPTOL);
    return(KIN_ILL_INPUT);
  }

  if (scsteptol == ZERO) {
    uround = kin_mem->kin_uround;
    kin_mem->kin_scsteptol = SUNRpowerR(uround,TWOTHIRDS);
  } else {
    kin_mem->kin_scsteptol = scsteptol;
  }

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetConstraints
 * -----------------------------------------------------------------
 */

int KINSetConstraints(void *kinmem, N_Vector constraints)
{
  KINMem kin_mem;
  realtype temptest;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetConstraints", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (constraints == NULL) {
    if (kin_mem->kin_constraintsSet) {
      N_VDestroy(kin_mem->kin_constraints);
      lrw -= lrw1;
      liw -= liw1;
    }
    kin_mem->kin_constraintsSet = FALSE;
    return(KIN_SUCCESS);
  }

  /* Check the constraints vector */

  temptest = N_VMaxNorm(constraints);
  if (temptest > TWOPT5){ 
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetConstraints", MSG_BAD_CONSTRAINTS);
    return(KIN_ILL_INPUT); 
  }

  if (!kin_mem->kin_constraintsSet) {
    kin_mem->kin_constraints = N_VClone(constraints);
    lrw += lrw1;
    liw += liw1;
    kin_mem->kin_constraintsSet = TRUE;
  }

  /* Load the constraint vector */

  N_VScale(ONE, constraints, kin_mem->kin_constraints);

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetSysFunc
 * -----------------------------------------------------------------
 */

int KINSetSysFunc(void *kinmem, KINSysFn func)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINSetSysFunc", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (func == NULL) {
    KINProcessError(NULL, KIN_ILL_INPUT, "KINSOL", "KINSetSysFunc", MSG_FUNC_NULL);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_func = func;

  return(KIN_SUCCESS);
}

/* 
 * =================================================================
 * Readability constants
 * =================================================================
 */

#define nni (kin_mem->kin_nni)
#define nfe (kin_mem->kin_nfe)
#define nbcf (kin_mem->kin_nbcf)  
#define nbktrk (kin_mem->kin_nbktrk)
#define stepl (kin_mem->kin_stepl)
#define fnorm (kin_mem->kin_fnorm)
#define liw (kin_mem->kin_liw)
#define lrw (kin_mem->kin_lrw)

/* 
 * =================================================================
 * KINSOL optional input functions
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : KINGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINGetWorkSpace", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  *lenrw = lrw;
  *leniw = liw;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumNonlinSolvIters
 * -----------------------------------------------------------------
 */

int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINGetNumNonlinSolvIters", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *nniters = nni;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINGetNumFuncEvals(void *kinmem, long int *nfevals)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINGetNumFuncEvals", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *nfevals = nfe;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumBetaCondFails
 * -----------------------------------------------------------------
 */

int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINGetNumBetaCondFails", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *nbcfails = nbcf;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumBacktrackOps
 * -----------------------------------------------------------------
 */

int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINGetNumBacktrackOps", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *nbacktr = nbktrk;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetFuncNorm
 * -----------------------------------------------------------------
 */

int KINGetFuncNorm(void *kinmem, realtype *funcnorm)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINGetFuncNorm", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *funcnorm = kin_mem->kin_fnorm;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetStepLength
 * -----------------------------------------------------------------
 */

int KINGetStepLength(void *kinmem, realtype *steplength)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    KINProcessError(NULL, KIN_MEM_NULL, "KINSOL", "KINGetStepLength", MSG_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *steplength = stepl;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *KINGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(24*sizeof(char));

  switch(flag) {
  case KIN_SUCCESS:
    sprintf(name, "KIN_SUCCESS");
    break;
  case KIN_INITIAL_GUESS_OK:
    sprintf(name, "KIN_INITIAL_GUESS_OK");
    break;
  case KIN_STEP_LT_STPTOL:
    sprintf(name, "KIN_STEP_LT_STPTOL");
    break;
  case KIN_WARNING:
    sprintf(name, "KIN_WARNING");
    break;
  case KIN_MEM_NULL:
    sprintf(name, "KIN_MEM_NULL");
    break;
  case KIN_ILL_INPUT:
    sprintf(name, "KIN_ILL_INPUT");
    break;
  case KIN_NO_MALLOC:
    sprintf(name, "KIN_NO_MALLOC");
    break;
  case KIN_MEM_FAIL:
    sprintf(name, "KIN_MEM_FAIL");
    break;
  case KIN_LINESEARCH_NONCONV:
    sprintf(name, "KIN_LINESEARCH_NONCONV");
    break;
  case KIN_MAXITER_REACHED:
    sprintf(name, "KIN_MAXITER_REACHED");
    break;
  case KIN_MXNEWT_5X_EXCEEDED:
    sprintf(name, "KIN_MXNEWT_5X_EXCEEDED");
    break;
  case KIN_LINESEARCH_BCFAIL:
    sprintf(name, "KIN_LINESEARCH_BCFAIL");
    break;
  case KIN_LINSOLV_NO_RECOVERY:
    sprintf(name, "KIN_LINSOLV_NO_RECOVERY");
    break;
  case KIN_LINIT_FAIL:
    sprintf(name, "KIN_LINIT_FAIL");
    break;
  case KIN_LSETUP_FAIL:
    sprintf(name, "KIN_LSETUP_FAIL");
    break;
  case KIN_LSOLVE_FAIL:
    sprintf(name, "KIN_LSOLVE_FAIL");
    break;
  default:
    sprintf(name, "NONE");
  }

  return(name);
}
