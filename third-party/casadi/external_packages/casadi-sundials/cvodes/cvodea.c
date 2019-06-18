/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
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
 * This is the implementation file for the CVODEA adjoint integrator.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_impl.h"

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

/* 
 * =================================================================
 * MACRO DEFINITIONS
 * =================================================================
 */

#define loop for(;;)

/* 
 * =================================================================
 * CVODEA PRIVATE CONSTANTS
 * =================================================================
 */

#define ZERO        RCONST(0.0)        /* real 0.0   */
#define ONE         RCONST(1.0)        /* real 1.0   */
#define TWO         RCONST(2.0)        /* real 2.0   */
#define HUNDRED     RCONST(100.0)      /* real 100.0 */
#define FUZZ_FACTOR RCONST(1000000.0)  /* fuzz factor for IMget */

/* 
 * =================================================================
 * PRIVATE FUNCTION PROTOTYPES
 * =================================================================
 */

static CkpntMem CVAckpntInit(CVodeMem cv_mem);
static CkpntMem CVAckpntNew(CVodeMem cv_mem);
static void CVAckpntDelete(CkpntMem *ck_memPtr);

static void CVAbckpbDelete(CVodeBMem *cvB_memPtr);

static int  CVAdataStore(CVodeMem cv_mem, CkpntMem ck_mem);
static int  CVAckpntGet(CVodeMem cv_mem, CkpntMem ck_mem); 

static int CVAfindIndex(CVodeMem cv_mem, realtype t, 
                        long int *indx, booleantype *newpoint);

static booleantype CVAhermiteMalloc(CVodeMem cv_mem);
static void CVAhermiteFree(CVodeMem cv_mem);
static int CVAhermiteGetY(CVodeMem cv_mem, realtype t, N_Vector y, N_Vector *yS);
static int CVAhermiteStorePnt(CVodeMem cv_mem, DtpntMem d);

static booleantype CVApolynomialMalloc(CVodeMem cv_mem);
static void CVApolynomialFree(CVodeMem cv_mem);
static int CVApolynomialGetY(CVodeMem cv_mem, realtype t, N_Vector y, N_Vector *yS);
static int CVApolynomialStorePnt(CVodeMem cv_mem, DtpntMem d);

/* Wrappers */

static int CVArhs(realtype t, N_Vector yB, 
                  N_Vector yBdot, void *cvode_mem);

static int CVArhsQ(realtype t, N_Vector yB, 
                   N_Vector qBdot, void *cvode_mem);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/*
 * CVodeAdjInit
 *
 * This routine initializes ASA and allocates space for the adjoint 
 * memory structure.
 */

int CVodeAdjInit(void *cvode_mem, long int steps, int interp)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;
  long int i, ii;

  /* ---------------
   * Check arguments
   * --------------- */

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeAdjInit", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  if (steps <= 0) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeAdjInit", MSGCV_BAD_STEPS);
    return(CV_ILL_INPUT);
  }

  if ( (interp != CV_HERMITE) && (interp != CV_POLYNOMIAL) ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeAdjInit", MSGCV_BAD_INTERP);
    return(CV_ILL_INPUT);
  } 

  /* ----------------------------
   * Allocate CVODEA memory block
   * ---------------------------- */

  ca_mem = NULL;
  ca_mem = (CVadjMem) malloc(sizeof(struct CVadjMemRec));
  if (ca_mem == NULL) {
    cvProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjInit", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  /* Attach ca_mem to CVodeMem structure */

  cv_mem->cv_adj_mem = ca_mem;

  /* ------------------------------
   * Initialization of check points
   * ------------------------------ */

  /* Set Check Points linked list to NULL */
  ca_mem->ck_mem = NULL;

  /* Initialize nckpnts to ZERO */
  ca_mem->ca_nckpnts = 0;

  /* No interpolation data is available */
  ca_mem->ca_ckpntData = NULL;

  /* ------------------------------------
   * Initialization of interpolation data
   * ------------------------------------ */

  /* Interpolation type */

  ca_mem->ca_IMtype = interp;

  /* Number of steps between check points */

  ca_mem->ca_nsteps = steps;

  /* Allocate space for the array of Data Point structures */

  ca_mem->dt_mem = NULL;
  ca_mem->dt_mem = (DtpntMem *) malloc((steps+1)*sizeof(struct DtpntMemRec *));
  if (ca_mem->dt_mem == NULL) {
    free(ca_mem); ca_mem = NULL;
    cvProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjInit", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  for (i=0; i<=steps; i++) { 
    ca_mem->dt_mem[i] = NULL;
    ca_mem->dt_mem[i] = (DtpntMem) malloc(sizeof(struct DtpntMemRec));
    if (ca_mem->dt_mem[i] == NULL) {
      for(ii=0; ii<i; ii++) {free(ca_mem->dt_mem[ii]); ca_mem->dt_mem[ii] = NULL;}
      free(ca_mem->dt_mem); ca_mem->dt_mem = NULL;
      free(ca_mem); ca_mem = NULL;
      cvProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjInit", MSGCV_MEM_FAIL);
      return(CV_MEM_FAIL);
    }
  }

  /* Attach functions for the appropriate interpolation module */
  
  switch(interp) {

  case CV_HERMITE:
    
    ca_mem->ca_IMmalloc = CVAhermiteMalloc;
    ca_mem->ca_IMfree   = CVAhermiteFree;
    ca_mem->ca_IMget    = CVAhermiteGetY;
    ca_mem->ca_IMstore  = CVAhermiteStorePnt;

    break;
    
  case CV_POLYNOMIAL:
  
    ca_mem->ca_IMmalloc = CVApolynomialMalloc;
    ca_mem->ca_IMfree   = CVApolynomialFree;
    ca_mem->ca_IMget    = CVApolynomialGetY;
    ca_mem->ca_IMstore  = CVApolynomialStorePnt;

    break;

  }

  /* The interpolation module has not been initialized yet */

  ca_mem->ca_IMmallocDone = FALSE;

  /* By default we will store but not interpolate sensitivities
   *  - IMstoreSensi will be set in CVodeF to FALSE if FSA is not enabled
   *    or if the user can force this through CVodeSetAdjNoSensi 
   *  - IMinterpSensi will be set in CVodeB to TRUE if IMstoreSensi is
   *    TRUE and if at least one backward problem requires sensitivities */

  ca_mem->ca_IMstoreSensi = TRUE;
  ca_mem->ca_IMinterpSensi = FALSE;

  /* ------------------------------------
   * Initialize list of backward problems
   * ------------------------------------ */

  ca_mem->cvB_mem = NULL;
  ca_mem->ca_bckpbCrt = NULL;
  ca_mem->ca_nbckpbs = 0;

  /* --------------------------------
   * CVodeF and CVodeB not called yet
   * -------------------------------- */

  ca_mem->ca_firstCVodeFcall = TRUE;
  ca_mem->ca_tstopCVodeFcall = FALSE;

  ca_mem->ca_firstCVodeBcall = TRUE;

  /* ---------------------------------------------
   * ASA initialized and allocated
   * --------------------------------------------- */

  cv_mem->cv_adj = TRUE;
  cv_mem->cv_adjMallocDone = TRUE;

  return(CV_SUCCESS);
} 

/* CVodeAdjReInit
 *
 * This routine reinitializes the CVODEA memory structure assuming that the
 * the number of steps between check points and the type of interpolation
 * remain unchanged.
 * The list of check points (and associated memory) is deleted.
 * The list of backward problems is kept (however, new backward problems can 
 * be added to this list by calling CVodeCreateB).
 * The CVODES memory for the forward and backward problems can be reinitialized
 * separately by calling CVodeReInit and CVodeReInitB, respectively.
 * NOTE: if a completely new list of backward problems is also needed, then
 *       simply free the adjoint memory (by calling CVodeAdjFree) and reinitialize
 *       ASA with CVodeAdjInit.
 */

int CVodeAdjReInit(void *cvode_mem)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;

  /* Check cvode_mem */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeAdjReInit", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeAdjReInit", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  /* Free current list of Check Points */

  while (ca_mem->ck_mem != NULL) CVAckpntDelete(&(ca_mem->ck_mem));

  /* Initialization of check points */
  
  ca_mem->ck_mem = NULL;
  ca_mem->ca_nckpnts = 0;
  ca_mem->ca_ckpntData = NULL;

  /* CVodeF and CVodeB not called yet */
 
  ca_mem->ca_firstCVodeFcall = TRUE;
  ca_mem->ca_tstopCVodeFcall = FALSE;
  ca_mem->ca_firstCVodeBcall = TRUE;

  return(CV_SUCCESS);
}

/*
 * CVodeAdjFree
 *
 * This routine frees the memory allocated by CVodeAdjInit.
 */

void CVodeAdjFree(void *cvode_mem)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  long int i;
  
  if (cvode_mem == NULL) return;
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_adjMallocDone) {

    ca_mem = cv_mem->cv_adj_mem;

    /* Delete check points one by one */
    while (ca_mem->ck_mem != NULL) CVAckpntDelete(&(ca_mem->ck_mem));

    /* Free vectors at all data points */
    if (ca_mem->ca_IMmallocDone) {
      ca_mem->ca_IMfree(cv_mem);
    }
    for(i=0; i<=ca_mem->ca_nsteps; i++) {
      free(ca_mem->dt_mem[i]);
      ca_mem->dt_mem[i] = NULL;
    }
    free(ca_mem->dt_mem);
    ca_mem->dt_mem = NULL;

    /* Delete backward problems one by one */
    while (ca_mem->cvB_mem != NULL) CVAbckpbDelete(&(ca_mem->cvB_mem));

    /* Free CVODEA memory */
    free(ca_mem);
    cv_mem->cv_adj_mem = NULL;

  }

}

/* 
 * -----------------------------------------------------------------
 * Readibility Constants
 * -----------------------------------------------------------------
 */

#define tinitial    (ca_mem->ca_tinitial)
#define tfinal      (ca_mem->ca_tfinal)
#define nckpnts     (ca_mem->ca_nckpnts)
#define nsteps      (ca_mem->ca_nsteps)
#define nbckpbs     (ca_mem->ca_nbckpbs)
#define ckpntData   (ca_mem->ca_ckpntData)
#define np          (ca_mem->ca_np)
#define ytmp        (ca_mem->ca_ytmp)
#define yStmp       (ca_mem->ca_yStmp)
#define Y           (ca_mem->ca_Y)
#define YS          (ca_mem->ca_YS)
#define T           (ca_mem->ca_T)

#define IMmalloc      (ca_mem->ca_IMmalloc)
#define IMfree        (ca_mem->ca_IMfree)
#define IMget         (ca_mem->ca_IMget)
#define IMstore       (ca_mem->ca_IMstore)
#define IMmallocDone  (ca_mem->ca_IMmallocDone)
#define IMstoreSensi  (ca_mem->ca_IMstoreSensi)
#define IMinterpSensi (ca_mem->ca_IMinterpSensi)
#define IMnewData     (ca_mem->ca_IMnewData)

#define uround     (cv_mem->cv_uround)
#define zn         (cv_mem->cv_zn)
#define nst        (cv_mem->cv_nst)
#define q          (cv_mem->cv_q)
#define qu         (cv_mem->cv_qu)
#define qprime     (cv_mem->cv_qprime)
#define qwait      (cv_mem->cv_qwait)
#define L          (cv_mem->cv_L)
#define gammap     (cv_mem->cv_gammap)
#define h          (cv_mem->cv_h)
#define hprime     (cv_mem->cv_hprime)
#define hscale     (cv_mem->cv_hscale)
#define eta        (cv_mem->cv_eta)
#define etamax     (cv_mem->cv_etamax)
#define tn         (cv_mem->cv_tn)
#define tretlast   (cv_mem->cv_tretlast)
#define tau        (cv_mem->cv_tau)
#define tq         (cv_mem->cv_tq)
#define l          (cv_mem->cv_l)
#define saved_tq5  (cv_mem->cv_saved_tq5)
#define forceSetup (cv_mem->cv_forceSetup)
#define f          (cv_mem->cv_f)
#define lmm        (cv_mem->cv_lmm)
#define iter       (cv_mem->cv_iter)
#define reltol     (cv_mem->cv_reltol)
#define user_data  (cv_mem->cv_user_data)
#define errfp      (cv_mem->cv_errfp)
#define h0u        (cv_mem->cv_h0u)
#define tempv      (cv_mem->cv_tempv)

#define quadr      (cv_mem->cv_quadr)
#define errconQ    (cv_mem->cv_errconQ)
#define znQ        (cv_mem->cv_znQ)
#define tempvQ     (cv_mem->cv_tempvQ)

#define sensi      (cv_mem->cv_sensi)
#define Ns         (cv_mem->cv_Ns)
#define errconS    (cv_mem->cv_errconS)
#define znS        (cv_mem->cv_znS)

#define quadr_sensi (cv_mem->cv_quadr_sensi)
#define errconQS    (cv_mem->cv_errconQS)
#define znQS        (cv_mem->cv_znQS)

#define t0_        (ck_mem->ck_t0)
#define t1_        (ck_mem->ck_t1)
#define zn_        (ck_mem->ck_zn)
#define znQ_       (ck_mem->ck_znQ)
#define znS_       (ck_mem->ck_znS)
#define znQS_      (ck_mem->ck_znQS)
#define quadr_     (ck_mem->ck_quadr)
#define sensi_     (ck_mem->ck_sensi)
#define quadr_sensi_ (ck_mem->ck_quadr_sensi)
#define Ns_        (ck_mem->ck_Ns)
#define zqm_       (ck_mem->ck_zqm)
#define nst_       (ck_mem->ck_nst)
#define tretlast_  (ck_mem->ck_tretlast)
#define q_         (ck_mem->ck_q)
#define qprime_    (ck_mem->ck_qprime)
#define qwait_     (ck_mem->ck_qwait)
#define L_         (ck_mem->ck_L)
#define gammap_    (ck_mem->ck_gammap)
#define h_         (ck_mem->ck_h)
#define hprime_    (ck_mem->ck_hprime)
#define hscale_    (ck_mem->ck_hscale)
#define eta_       (ck_mem->ck_eta)
#define etamax_    (ck_mem->ck_etamax)
#define tau_       (ck_mem->ck_tau)
#define tq_        (ck_mem->ck_tq)
#define l_         (ck_mem->ck_l)
#define saved_tq5_ (ck_mem->ck_saved_tq5)
#define next_      (ck_mem->ck_next)


/*
 * CVodeF
 *
 * This routine integrates to tout and returns solution into yout.
 * In the same time, it stores check point data every 'steps' steps. 
 * 
 * CVodeF can be called repeatedly by the user.
 *
 * ncheckPtr points to the number of check points stored so far.
 */

int CVodeF(void *cvode_mem, realtype tout, N_Vector yout, 
           realtype *tret, int itask, int *ncheckPtr)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;
  CkpntMem tmp;
  DtpntMem *dt_mem;
  int flag, i;
  booleantype iret, allocOK;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeF", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeF", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  /* Check for yout != NULL */
  if (yout == NULL) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeF", MSGCV_YOUT_NULL);
    return(CV_ILL_INPUT);
  }
  
  /* Check for tret != NULL */
  if (tret == NULL) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeF", MSGCV_TRET_NULL);
    return(CV_ILL_INPUT);
  }

  /* Check for valid itask */
  if ( (itask != CV_NORMAL) && (itask != CV_ONE_STEP) ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeF", MSGCV_BAD_ITASK);
    return(CV_ILL_INPUT);
  }

  /* All error checking done */

  dt_mem = ca_mem->dt_mem;

  /* If tstop is enabled, store some info */
  if (cv_mem->cv_tstopset) {
    ca_mem->ca_tstopCVodeFcall = TRUE;
    ca_mem->ca_tstopCVodeF = cv_mem->cv_tstop;
  }

  /* We will call CVode in CV_ONE_STEP mode, regardless
   * of what itask is, so flag if we need to return */
  if (itask == CV_ONE_STEP) iret = TRUE;
  else                      iret = FALSE;

  /* On the first step:
   *   - set tinitial
   *   - initialize list of check points
   *   - if needed, initialize the interpolation module
   *   - load dt_mem[0]
   * On subsequent steps, test if taking a new step is necessary. 
   */
  if ( ca_mem->ca_firstCVodeFcall ) {

    tinitial = tn;

    ca_mem->ck_mem = CVAckpntInit(cv_mem);
    if (ca_mem->ck_mem == NULL) {
      cvProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeF", MSGCV_MEM_FAIL);
      return(CV_MEM_FAIL);
    }

    if ( !IMmallocDone ) {

      /* Do we need to store sensitivities? */
      if (!sensi) IMstoreSensi = FALSE;

      /* Allocate space for interpolation data */
      allocOK = IMmalloc(cv_mem);
      if (!allocOK) {
        cvProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeF", MSGCV_MEM_FAIL);
        return(CV_MEM_FAIL);
      }

      /* Rename zn and, if needed, znS for use in interpolation */
      for (i=0;i<L_MAX;i++) Y[i] = zn[i];
      if (IMstoreSensi) {
        for (i=0;i<L_MAX;i++) YS[i] = znS[i];
      }

      IMmallocDone = TRUE;

    }

    dt_mem[0]->t = ca_mem->ck_mem->ck_t0;
    IMstore(cv_mem, dt_mem[0]);

    ca_mem->ca_firstCVodeFcall = FALSE;

  } else if ( (tn - tout)*h >= ZERO ) {

    /* If tout was passed, return interpolated solution. 
       No changes to ck_mem or dt_mem are needed. */
    *tret = tout;
    flag = CVodeGetDky(cv_mem, tout, 0, yout);
    *ncheckPtr = nckpnts;
    IMnewData = TRUE;
    ckpntData = ca_mem->ck_mem;
    np = nst % nsteps + 1;

    return(flag);

  }

  /* Integrate to tout (in CV_ONE_STEP mode) while loading check points */
  loop {

    /* Perform one step of the integration */

    flag = CVode(cv_mem, tout, yout, tret, CV_ONE_STEP);
    if (flag < 0) break;

    /* Test if a new check point is needed */

    if ( nst % nsteps == 0 ) {

      ca_mem->ck_mem->ck_t1 = *tret;

      /* Create a new check point, load it, and append it to the list */
      tmp = CVAckpntNew(cv_mem);
      if (tmp == NULL) {
        cvProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeF", MSGCV_MEM_FAIL);
        flag = CV_MEM_FAIL;
        break;
      }
      tmp->ck_next = ca_mem->ck_mem;
      ca_mem->ck_mem = tmp;
      nckpnts++;
      forceSetup = TRUE;
      
      /* Reset i=0 and load dt_mem[0] */
      dt_mem[0]->t = ca_mem->ck_mem->ck_t0;
      IMstore(cv_mem, dt_mem[0]);

    } else {

      /* Load next point in dt_mem */
      dt_mem[nst%nsteps]->t = *tret;
      IMstore(cv_mem, dt_mem[nst%nsteps]);

    }

    /* Set t1 field of the current ckeck point structure
       for the case in which there will be no future
       check points */
    ca_mem->ck_mem->ck_t1 = *tret;

    /* tfinal is now set to *tret */
    tfinal = *tret;

    /* Return if in CV_ONE_STEP mode */
    if (iret) break;

    /* Return if tout reached */
    if ( (*tret - tout)*h >= ZERO ) {
      *tret = tout;
      CVodeGetDky(cv_mem, tout, 0, yout);
      /* Reset tretlast in cv_mem so that CVodeGetQuad and CVodeGetSens 
       * evaluate quadratures and/or sensitivities at the proper time */
      cv_mem->cv_tretlast = tout;
      break;
    }

  } /* end of loop() */

  /* Get ncheck from ca_mem */ 
  *ncheckPtr = nckpnts;

  /* Data is available for the last interval */
  IMnewData = TRUE;
  ckpntData = ca_mem->ck_mem;
  np = nst % nsteps + 1;

  return(flag);
}



/* 
 * =================================================================
 * FUNCTIONS FOR BACKWARD PROBLEMS
 * =================================================================
 */


int CVodeCreateB(void *cvode_mem, int lmmB, int iterB, int *which)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem new_cvB_mem;
  void *cvodeB_mem;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeCreateB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeCreateB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  }
  ca_mem = cv_mem->cv_adj_mem;

  /* Allocate space for new CVodeBMem object */

  new_cvB_mem = NULL;
  new_cvB_mem = (CVodeBMem) malloc(sizeof(struct CVodeBMemRec));
  if (new_cvB_mem == NULL) {
    cvProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeCreateB", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  /* Create and set a new CVODES object for the backward problem */

  cvodeB_mem = CVodeCreate(lmmB, iterB);
  if (cvodeB_mem == NULL) {
    cvProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeCreateB", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  CVodeSetUserData(cvodeB_mem, cvode_mem);

  CVodeSetMaxHnilWarns(cvodeB_mem, -1);

  CVodeSetErrHandlerFn(cvodeB_mem, cv_mem->cv_ehfun, cv_mem->cv_eh_data);
  CVodeSetErrFile(cvodeB_mem, cv_mem->cv_errfp);

  /* Set/initialize fields in the new CVodeBMem object, new_cvB_mem */

  new_cvB_mem->cv_index   = nbckpbs;

  new_cvB_mem->cv_mem     = (CVodeMem) cvodeB_mem;

  new_cvB_mem->cv_f       = NULL;
  new_cvB_mem->cv_fs      = NULL;

  new_cvB_mem->cv_fQ      = NULL;
  new_cvB_mem->cv_fQs     = NULL;

  new_cvB_mem->cv_user_data  = NULL;

  new_cvB_mem->cv_lmem    = NULL;
  new_cvB_mem->cv_lfree   = NULL;
  new_cvB_mem->cv_pmem    = NULL;
  new_cvB_mem->cv_pfree   = NULL;

  new_cvB_mem->cv_y       = NULL;

  new_cvB_mem->cv_f_withSensi = FALSE;
  new_cvB_mem->cv_fQ_withSensi = FALSE;

  /* Attach the new object to the linked list cvB_mem */

  new_cvB_mem->cv_next = ca_mem->cvB_mem;
  ca_mem->cvB_mem = new_cvB_mem;
  
  /* Return the index of the newly created CVodeBMem object.
   * This must be passed to CVodeInitB and to other ***B 
   * functions to set optional inputs for this backward problem */

  *which = nbckpbs;

  nbckpbs++;

  return(CV_SUCCESS);
}

int CVodeInitB(void *cvode_mem, int which, 
               CVRhsFnB fB,
               realtype tB0, N_Vector yB0)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeInitB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */

  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeInitB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */

  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeInitB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */

  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);
  
  /* Allocate and set the CVODES object */

  flag = CVodeInit(cvodeB_mem, CVArhs, tB0, yB0);

  if (flag != CV_SUCCESS) return(flag);

  /* Copy fB function in cvB_mem */

  cvB_mem->cv_f_withSensi = FALSE;
  cvB_mem->cv_f = fB;

  /* Allocate space and initialize the y Nvector in cvB_mem */

  cvB_mem->cv_t0 = tB0;
  cvB_mem->cv_y = N_VClone(yB0);
  N_VScale(ONE, yB0, cvB_mem->cv_y);

  return(CV_SUCCESS);
}

int CVodeInitBS(void *cvode_mem, int which, 
                CVRhsFnBS fBs,
                realtype tB0, N_Vector yB0)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeInitBS", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */

  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeInitBS", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */

  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeInitBS", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */

  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);
  
  /* Allocate and set the CVODES object */

  flag = CVodeInit(cvodeB_mem, CVArhs, tB0, yB0);

  if (flag != CV_SUCCESS) return(flag);

  /* Copy fBs function in cvB_mem */

  cvB_mem->cv_f_withSensi = TRUE;
  cvB_mem->cv_fs = fBs;

  /* Allocate space and initialize the y Nvector in cvB_mem */

  cvB_mem->cv_t0 = tB0;
  cvB_mem->cv_y = N_VClone(yB0);
  N_VScale(ONE, yB0, cvB_mem->cv_y);

  return(CV_SUCCESS);
}


int CVodeReInitB(void *cvode_mem, int which,
                 realtype tB0, N_Vector yB0)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeReInitB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeReInitB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  }
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeReInitB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Reinitialize CVODES object */

  flag = CVodeReInit(cvodeB_mem, tB0, yB0);

  return(flag);
}


int CVodeSStolerancesB(void *cvode_mem, int which, realtype reltolB, realtype abstolB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeSStolerancesB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */

  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeSStolerancesB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */

  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeSStolerancesB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */

  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Set tolerances */

  flag = CVodeSStolerances(cvodeB_mem, reltolB, abstolB);

  return(flag);
}


int CVodeSVtolerancesB(void *cvode_mem, int which, realtype reltolB, N_Vector abstolB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeSVtolerancesB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */

  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeSVtolerancesB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */

  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeSVtolerancesB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */

  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Set tolerances */

  flag = CVodeSVtolerances(cvodeB_mem, reltolB, abstolB);

  return(flag);
}


int CVodeQuadInitB(void *cvode_mem, int which,
                     CVQuadRhsFnB fQB, N_Vector yQB0)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeQuadInitB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeQuadInitB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeQuadInitB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVodeQuadInit(cvodeB_mem, CVArhsQ, yQB0);
  if (flag != CV_SUCCESS) return(flag);

  cvB_mem->cv_fQ_withSensi = FALSE;
  cvB_mem->cv_fQ = fQB;

  return(CV_SUCCESS);
}

int CVodeQuadInitBS(void *cvode_mem, int which,
                      CVQuadRhsFnBS fQBs, N_Vector yQB0)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeQuadInitBS", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeQuadInitBS", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeQuadInitBS", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVodeQuadInit(cvodeB_mem, CVArhsQ, yQB0);
  if (flag != CV_SUCCESS) return(flag);

  cvB_mem->cv_fQ_withSensi = TRUE;
  cvB_mem->cv_fQs = fQBs;

  return(CV_SUCCESS);
}

int CVodeQuadReInitB(void *cvode_mem, int which, N_Vector yQB0)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeQuadReInitB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeQuadReInitB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeQuadReInitB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVodeQuadReInit(cvodeB_mem, yQB0);
  if (flag != CV_SUCCESS) return(flag);

  return(CV_SUCCESS);
}

int CVodeQuadSStolerancesB(void *cvode_mem, int which, realtype reltolQB, realtype abstolQB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeQuadSStolerancesB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeQuadSStolerancesB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeQuadSStolerancesB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVodeQuadSStolerances(cvodeB_mem, reltolQB, abstolQB);

  return(flag);
}

int CVodeQuadSVtolerancesB(void *cvode_mem, int which, realtype reltolQB, N_Vector abstolQB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeQuadSStolerancesB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeQuadSStolerancesB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeQuadSStolerancesB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVodeQuadSVtolerances(cvodeB_mem, reltolQB, abstolQB);

  return(flag);
}

/*
 * CVodeB
 *
 * This routine performs the backward integration towards tBout
 * of all backward problems that were defined.
 * When necessary, it performs a forward integration between two 
 * consecutive check points to update interpolation data.
 *
 * On a successful return, CVodeB returns CV_SUCCESS.
 *
 * NOTE that CVodeB DOES NOT return the solution for the backward
 * problem(s). Use CVodeGetB to extract the solution at tBret
 * for any given backward problem.
 *
 * If there are multiple backward problems and multiple check points,
 * CVodeB may not succeed in getting all problems to take one step
 * when called in ONE_STEP mode.
 */

int CVodeB(void *cvode_mem, realtype tBout, int itaskB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem, tmp_cvB_mem;
  CkpntMem ck_mem;
  int sign, flag=0;
  realtype tfuzz, tBret, tBn;
  booleantype gotCheckpoint, isActive, reachedTBout;
  
  /* Check if cvode_mem exists */

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */

  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  }
  ca_mem = cv_mem->cv_adj_mem;

  /* Check if any backward problem has been defined */

  if ( nbckpbs == 0 ) {
    cvProcessError(cv_mem, CV_NO_BCK, "CVODEA", "CVodeB", MSGCV_NO_BCK);
    return(CV_NO_BCK);
  }
  cvB_mem = ca_mem->cvB_mem;

  /* Check whether CVodeF has been called */

  if ( ca_mem->ca_firstCVodeFcall ) {
    cvProcessError(cv_mem, CV_NO_FWD, "CVODEA", "CVodeB", MSGCV_NO_FWD);
    return(CV_NO_FWD);
  }
  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* If this is the first call, loop over all backward problems and
   *   - check that tB0 is valid
   *   - check that tBout is ahead of tB0 in the backward direction
   *   - check whether we need to interpolate forward sensitivities
   */

  if ( ca_mem->ca_firstCVodeBcall ) {

    tmp_cvB_mem = cvB_mem;

    while(tmp_cvB_mem != NULL) {

      tBn = tmp_cvB_mem->cv_mem->cv_tn;

      if ( (sign*(tBn-tinitial) < ZERO) || (sign*(tfinal-tBn) < ZERO) ) {
        cvProcessError(cv_mem, CV_BAD_TB0, "CVODEA", "CVodeB", MSGCV_BAD_TB0,
                       tmp_cvB_mem->cv_index);
        return(CV_BAD_TB0);
      }

      if (sign*(tBn-tBout) <= ZERO) {
        cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeB", MSGCV_BAD_TBOUT,
                       tmp_cvB_mem->cv_index);
        return(CV_ILL_INPUT);
      }

      if ( tmp_cvB_mem->cv_f_withSensi || tmp_cvB_mem->cv_fQ_withSensi )
          IMinterpSensi = TRUE;

      tmp_cvB_mem = tmp_cvB_mem->cv_next;

    }

    if ( IMinterpSensi && !IMstoreSensi) {
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeB", MSGCV_BAD_SENSI);
      return(CV_ILL_INPUT);
    }

    ca_mem->ca_firstCVodeBcall = FALSE;
  }

  /* Check if itaskB is legal */

  if ( (itaskB != CV_NORMAL) && (itaskB != CV_ONE_STEP) ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeB", MSGCV_BAD_ITASKB);
    return(CV_ILL_INPUT);
  }

  /* Check if tBout is legal */

  if ( (sign*(tBout-tinitial) < ZERO) || (sign*(tfinal-tBout) < ZERO) ) {
    tfuzz = HUNDRED*uround*(SUNRabs(tinitial) + SUNRabs(tfinal));
    if ( (sign*(tBout-tinitial) < ZERO) && (SUNRabs(tBout-tinitial) < tfuzz) ) {
      tBout = tinitial;
    } else {
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeB", MSGCV_BAD_TBOUT);
      return(CV_ILL_INPUT);
    }
  }

  /* Loop through the check points and stop as soon as a backward
   * problem has its tn value behind the current check point's t0_
   * value (in the backward direction) */

  ck_mem = ca_mem->ck_mem;

  gotCheckpoint = FALSE;

  loop {

    tmp_cvB_mem = cvB_mem;
    while(tmp_cvB_mem != NULL) {
      tBn = tmp_cvB_mem->cv_mem->cv_tn;

      if ( sign*(tBn-t0_) > ZERO ) {
        gotCheckpoint = TRUE;
        break;
      }

      if ( (itaskB==CV_NORMAL) && (tBn == t0_) && (sign*(tBout-t0_) >= ZERO) ) {
        gotCheckpoint = TRUE;
        break;
      }

      tmp_cvB_mem = tmp_cvB_mem->cv_next;
    }

    if (gotCheckpoint) break;

    if (next_ == NULL) break;

    ck_mem = next_;
  }

  /* Starting with the current check point from above, loop over check points
     while propagating backward problems */

  loop {

    /* Store interpolation data if not available.
       This is the 2nd forward integration pass */

    if (ck_mem != ckpntData) {
      flag = CVAdataStore(cv_mem, ck_mem);
      if (flag != CV_SUCCESS) break;
    }

    /* Loop through all backward problems and, if needed,
     * propagate their solution towards tBout */

    tmp_cvB_mem = cvB_mem;
    while (tmp_cvB_mem != NULL) {

      /* Decide if current backward problem is "active" in this check point */

      isActive = TRUE;

      tBn = tmp_cvB_mem->cv_mem->cv_tn;

      if ( (tBn == t0_) && (sign*(tBout-t0_) < ZERO ) ) isActive = FALSE;
      if ( (tBn == t0_) && (itaskB==CV_ONE_STEP) ) isActive = FALSE;

      if ( sign * (tBn - t0_) < ZERO ) isActive = FALSE;

      if ( isActive ) {

        /* Store the address of current backward problem memory 
         * in ca_mem to be used in the wrapper functions */
        ca_mem->ca_bckpbCrt = tmp_cvB_mem;

        /* Integrate current backward problem */
        CVodeSetStopTime(tmp_cvB_mem->cv_mem, t0_);
        flag = CVode(tmp_cvB_mem->cv_mem, tBout, tmp_cvB_mem->cv_y, &tBret, itaskB);

        /* Set the time at which we will report solution and/or quadratures */
        tmp_cvB_mem->cv_tout = tBret;

        /* If an error occurred, exit while loop */
        if (flag < 0) break;

      } else {
        flag = CV_SUCCESS;
        tmp_cvB_mem->cv_tout = tBn;
      }

      /* Move to next backward problem */

      tmp_cvB_mem = tmp_cvB_mem->cv_next;
    }

    /* If an error occurred, return now */

    if (flag <0) {
      cvProcessError(cv_mem, flag, "CVODEA", "CVodeB", MSGCV_BACK_ERROR,
                     tmp_cvB_mem->cv_index);
      return(flag);
    }

    /* If in CV_ONE_STEP mode, return now (flag = CV_SUCCESS) */

    if (itaskB == CV_ONE_STEP) break;

    /* If all backward problems have succesfully reached tBout, return now */

    reachedTBout = TRUE;

    tmp_cvB_mem = cvB_mem;
    while(tmp_cvB_mem != NULL) {
      if ( sign*(tmp_cvB_mem->cv_tout - tBout) > ZERO ) {
        reachedTBout = FALSE;
        break;
      }
      tmp_cvB_mem = tmp_cvB_mem->cv_next;
    }

    if ( reachedTBout ) break;

    /* Move check point in linked list to next one */

    ck_mem = next_;

  } 

  return(flag);
}


int CVodeGetB(void *cvode_mem, int which, realtype *tret, N_Vector yB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeGetB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeGetB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeGetB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  } 

  N_VScale(ONE, cvB_mem->cv_y, yB);
  *tret = cvB_mem->cv_tout;

  return(CV_SUCCESS);
}


/*
 * CVodeGetQuadB
 */

int CVodeGetQuadB(void *cvode_mem, int which, realtype *tret, N_Vector qB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  long int nstB;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeGetQuadB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeGetQuadB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeGetQuadB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  } 

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* If the integration for this backward problem has not started yet,
   * simply return the current value of qB (i.e. the final conditions) */

  flag = CVodeGetNumSteps(cvodeB_mem, &nstB);
  
  if (nstB == 0) {
    N_VScale(ONE, cvB_mem->cv_mem->cv_znQ[0], qB);
    *tret = cvB_mem->cv_tout;
  } else {
    flag = CVodeGetQuad(cvodeB_mem, tret, qB);
  }

  return(flag);
}


/* 
 * =================================================================
 * PRIVATE FUNCTIONS FOR CHECK POINTS
 * =================================================================
 */

/*
 * CVAckpntInit
 *
 * This routine initializes the check point linked list with 
 * information from the initial time.
 */

static CkpntMem CVAckpntInit(CVodeMem cv_mem)
{
  CkpntMem ck_mem;
  int is;

  /* Allocate space for ckdata */
  ck_mem = NULL;
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  if (ck_mem == NULL) return(NULL);

  zn_[0] = N_VClone(tempv);
  if (zn_[0] == NULL) {
    free(ck_mem); ck_mem = NULL;
    return(NULL);
  }
  
  zn_[1] = N_VClone(tempv);
  if (zn_[1] == NULL) {
    N_VDestroy(zn_[0]);
    free(ck_mem); ck_mem = NULL;
    return(NULL);
  }

  /* zn_[qmax] was not allocated */
  zqm_ = 0;

  /* Load ckdata from cv_mem */
  N_VScale(ONE, zn[0], zn_[0]);
  t0_    = tn;
  nst_   = 0;
  q_     = 1;
  h_     = 0.0;
  
  /* Do we need to carry quadratures */
  quadr_ = quadr && errconQ;

  if (quadr_) {

    znQ_[0] = N_VClone(tempvQ);
    if (znQ_[0] == NULL) {
      N_VDestroy(zn_[0]);
      N_VDestroy(zn_[1]);
      free(ck_mem); ck_mem = NULL;
      return(NULL);
    }

    N_VScale(ONE, znQ[0], znQ_[0]);

  }

  /* Do we need to carry sensitivities? */
  sensi_ = sensi;

  if (sensi_) {

    Ns_ = Ns;

    znS_[0] = N_VCloneVectorArray(Ns, tempv);
    if (znS_[0] == NULL) {
      N_VDestroy(zn_[0]);
      N_VDestroy(zn_[1]);
      if (quadr_) N_VDestroy(znQ_[0]);
      free(ck_mem); ck_mem = NULL;
      return(NULL);
    }

    for (is=0; is<Ns; is++)
      N_VScale(ONE, znS[0][is], znS_[0][is]);

  }

  /* Do we need to carry quadrature sensitivities? */
  quadr_sensi_ = quadr_sensi && errconQS;

  if (quadr_sensi_) {
    znQS_[0] = N_VCloneVectorArray(Ns, tempvQ);
    if (znQS_[0] == NULL) {
      N_VDestroy(zn_[0]);
      N_VDestroy(zn_[1]);
      if (quadr_) N_VDestroy(znQ_[0]);
      N_VDestroyVectorArray(znS_[0], Ns);
      free(ck_mem); ck_mem = NULL;
      return(NULL);
    }
    
    for (is=0; is<Ns; is++)
      N_VScale(ONE, znQS[0][is], znQS_[0][is]);

  }

  /* Next in list */
  next_  = NULL;

  return(ck_mem);
}

/*
 * CVAckpntNew
 *
 * This routine allocates space for a new check point and sets 
 * its data from current values in cv_mem.
 */

static CkpntMem CVAckpntNew(CVodeMem cv_mem)
{
  CkpntMem ck_mem;
  int j, jj, is, qmax; 

  /* Allocate space for ckdata */
  ck_mem = NULL;
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  if (ck_mem == NULL) return(NULL);

  /* Set cv_next to NULL */
  ck_mem->ck_next = NULL;

  /* Test if we need to allocate space for the last zn.
   * NOTE: zn(qmax) may be needed for a hot restart, if an order
   * increase is deemed necessary at the first step after a check point */
  qmax = cv_mem->cv_qmax;
  zqm_ = (q < qmax) ? qmax : 0;

  for (j=0; j<=q; j++) {
    zn_[j] = N_VClone(tempv);
    if (zn_[j] == NULL) {
      for (jj=0; jj<j; jj++) N_VDestroy(zn_[jj]);
      free(ck_mem); ck_mem = NULL;
      return(NULL);
    }
  }

  if (q < qmax) {
    zn_[qmax] = N_VClone(tempv);
    if (zn_[qmax] == NULL) {
      for (jj=0; jj<=q; jj++) N_VDestroy(zn_[jj]);
      free(ck_mem); ck_mem = NULL;
      return(NULL);
    }
  }

  /* Test if we need to carry quadratures */
  quadr_ = quadr && errconQ;

  if (quadr_) {

    for (j=0; j<=q; j++) {
      znQ_[j] = N_VClone(tempvQ);
      if(znQ_[j] == NULL) {
        for (jj=0; jj<j; jj++) N_VDestroy(znQ_[jj]);
        if (q < qmax) N_VDestroy(zn_[qmax]);
        for (jj=0; jj<=q; j++) N_VDestroy(zn_[jj]);
        free(ck_mem); ck_mem = NULL;
        return(NULL);
      }
    }

    if (q < qmax) {
      znQ_[qmax] = N_VClone(tempvQ);
      if (znQ_[qmax] == NULL) {
        for (jj=0; jj<=q; jj++) N_VDestroy(znQ_[jj]);
        N_VDestroy(zn_[qmax]);
        for (jj=0; jj<=q; jj++) N_VDestroy(zn_[jj]);
        free(ck_mem); ck_mem = NULL;
        return(NULL);
      }
    }

  }

  /* Test if we need to carry sensitivities */
  sensi_ = sensi;

  if (sensi_) {

    Ns_ = Ns;

    for (j=0; j<=q; j++) {
      znS_[j] = N_VCloneVectorArray(Ns, tempv);
      if (znS_[j] == NULL) {
        for (jj=0; jj<j; jj++) N_VDestroyVectorArray(znS_[jj], Ns);
        if (quadr_) {
          if (q < qmax) N_VDestroy(znQ_[qmax]);
          for (jj=0; jj<=q; jj++) N_VDestroy(znQ_[jj]);
        }
        if (q < qmax) N_VDestroy(zn_[qmax]);
        for (jj=0; jj<=q; jj++) N_VDestroy(zn_[jj]);
        free(ck_mem); ck_mem = NULL;
        return(NULL);
      }
    }

    if ( q < qmax) {
      znS_[qmax] = N_VCloneVectorArray(Ns, tempv);
      if (znS_[qmax] == NULL) {
        for (jj=0; jj<=q; jj++) N_VDestroyVectorArray(znS_[jj], Ns);
        if (quadr_) {
          N_VDestroy(znQ_[qmax]);
          for (jj=0; jj<=q; jj++) N_VDestroy(znQ_[jj]);
        }
        N_VDestroy(zn_[qmax]);
        for (jj=0; jj<=q; jj++) N_VDestroy(zn_[jj]);
        free(ck_mem); ck_mem = NULL;
        return(NULL);
      }
    }

  }

  /* Test if we need to carry quadrature sensitivities */
  quadr_sensi_ = quadr_sensi && errconQS;

  if (quadr_sensi_) {

    for (j=0; j<=q; j++) {
      znQS_[j] = N_VCloneVectorArray(Ns, tempvQ);
      if (znQS_[j] == NULL) {
        for (jj=0; jj<j; jj++) N_VDestroyVectorArray(znQS_[jj], Ns);
        if (q < qmax) N_VDestroyVectorArray(znS_[qmax], Ns);
        for (jj=0; jj<=q; jj++) N_VDestroyVectorArray(znS_[jj], Ns);
        if (quadr_) {
          if (q < qmax) N_VDestroy(znQ_[qmax]);
          for (jj=0; jj<=q; jj++) N_VDestroy(znQ_[jj]);
        }
        if (q < qmax) N_VDestroy(zn_[qmax]);
        for (jj=0; jj<=q; jj++) N_VDestroy(zn_[jj]);
        free(ck_mem); ck_mem = NULL;
        return(NULL);
      }
    }

    if ( q < qmax) {
      znQS_[qmax] = N_VCloneVectorArray(Ns, tempvQ);
      if (znQS_[qmax] == NULL) {
        for (jj=0; jj<=q; jj++) N_VDestroyVectorArray(znQS_[jj], Ns);
        N_VDestroyVectorArray(znS_[qmax], Ns);
        for (jj=0; jj<=q; jj++) N_VDestroyVectorArray(znS_[jj], Ns);
        if (quadr_) {
          N_VDestroy(znQ_[qmax]);
          for (jj=0; jj<=q; jj++) N_VDestroy(zn_[jj]);
        }
        N_VDestroy(zn_[qmax]);
        for (jj=0; jj<=q; jj++) N_VDestroy(zn_[jj]);
        free(ck_mem); ck_mem = NULL;
        return(NULL);
      }
    }

  }

  /* Load check point data from cv_mem */

  for (j=0; j<=q; j++) N_VScale(ONE, zn[j], zn_[j]);
  if ( q < qmax ) N_VScale(ONE, zn[qmax], zn_[qmax]);

  if (quadr_) {
    for (j=0; j<=q; j++) N_VScale(ONE, znQ[j], znQ_[j]);
    if ( q < qmax ) N_VScale(ONE, znQ[qmax], znQ_[qmax]);
  }

  if (sensi_) {
    for (is=0; is<Ns; is++) {
      for (j=0; j<=q; j++) N_VScale(ONE, znS[j][is], znS_[j][is]);
      if ( q < qmax ) N_VScale(ONE, znS[qmax][is], znS_[qmax][is]);
    }
  }

  if (quadr_sensi_) {
    for (is=0; is<Ns; is++) {
      for (j=0; j<=q; j++) N_VScale(ONE, znQS[j][is], znQS_[j][is]);
      if ( q < qmax ) N_VScale(ONE, znQS[qmax][is], znQS_[qmax][is]);
    }
  }

  for (j=0; j<=L_MAX; j++)     tau_[j] = tau[j];
  for (j=0; j<=NUM_TESTS; j++) tq_[j] = tq[j];
  for (j=0; j<=q; j++)         l_[j] = l[j];
  nst_       = nst;
  tretlast_  = tretlast;
  q_         = q;
  qprime_    = qprime;
  qwait_     = qwait;
  L_         = L;
  gammap_    = gammap;
  h_         = h;
  hprime_    = hprime;
  hscale_    = hscale;
  eta_       = eta;
  etamax_    = etamax;
  t0_        = tn;
  saved_tq5_ = saved_tq5;

  return(ck_mem);
}

/*
 * CVAckpntDelete
 *
 * This routine deletes the first check point in list and returns
 * the new list head
 */

static void CVAckpntDelete(CkpntMem *ck_memPtr)
{
  CkpntMem tmp;
  int j;

  if (*ck_memPtr == NULL) return;

  /* store head of list */
  tmp = *ck_memPtr;

  /* move head of list */
  *ck_memPtr = (*ck_memPtr)->ck_next;

  /* free N_Vectors in tmp */
  for (j=0;j<=tmp->ck_q;j++) N_VDestroy(tmp->ck_zn[j]);
  if (tmp->ck_zqm != 0) N_VDestroy(tmp->ck_zn[tmp->ck_zqm]);

  /* free N_Vectors for quadratures in tmp 
   * Note that at the check point at t_initial, only znQ_[0] 
   * was allocated */
  if (tmp->ck_quadr) {

    if (tmp->ck_next != NULL) {
      for (j=0;j<=tmp->ck_q;j++) N_VDestroy(tmp->ck_znQ[j]);
      if (tmp->ck_zqm != 0) N_VDestroy(tmp->ck_znQ[tmp->ck_zqm]);
    } else {
      N_VDestroy(tmp->ck_znQ[0]);
    }
    
  }

  /* free N_Vectors for sensitivities in tmp
   * Note that at the check point at t_initial, only znS_[0] 
   * was allocated */
  if (tmp->ck_sensi) {
    
    if (tmp->ck_next != NULL) {
      for (j=0;j<=tmp->ck_q;j++) N_VDestroyVectorArray(tmp->ck_znS[j], tmp->ck_Ns);
      if (tmp->ck_zqm != 0) N_VDestroyVectorArray(tmp->ck_znS[tmp->ck_zqm], tmp->ck_Ns);
    } else {
      N_VDestroyVectorArray(tmp->ck_znS[0], tmp->ck_Ns);
    }
    
  }

  /* free N_Vectors for quadrature sensitivities in tmp
   * Note that at the check point at t_initial, only znQS_[0] 
   * was allocated */
  if (tmp->ck_quadr_sensi) {
    
    if (tmp->ck_next != NULL) {
      for (j=0;j<=tmp->ck_q;j++) N_VDestroyVectorArray(tmp->ck_znQS[j], tmp->ck_Ns);
      if (tmp->ck_zqm != 0) N_VDestroyVectorArray(tmp->ck_znQS[tmp->ck_zqm], tmp->ck_Ns);
    } else {
      N_VDestroyVectorArray(tmp->ck_znQS[0], tmp->ck_Ns);
    }
    
  }

  free(tmp); tmp = NULL;

}

/* 
 * =================================================================
 * PRIVATE FUNCTIONS FOR BACKWARD PROBLEMS
 * =================================================================
 */

static void CVAbckpbDelete(CVodeBMem *cvB_memPtr)
{
  CVodeBMem tmp;
  void *cvode_mem;

  if (*cvB_memPtr != NULL) {

    /* Save head of the list */
    tmp = *cvB_memPtr;

    /* Move head of the list */
    *cvB_memPtr = (*cvB_memPtr)->cv_next;

    /* Free CVODES memory in tmp */
    cvode_mem = (void *)(tmp->cv_mem);
    CVodeFree(&cvode_mem);

    /* Free linear solver memory */
    if (tmp->cv_lfree != NULL) tmp->cv_lfree(tmp);

    /* Free preconditioner memory */
    if (tmp->cv_pfree != NULL) tmp->cv_pfree(tmp);

    /* Free workspace Nvector */
    N_VDestroy(tmp->cv_y);

    free(tmp); tmp = NULL;

  }

}

/* 
 * =================================================================
 * PRIVATE FUNCTIONS FOR INTERPOLATION
 * =================================================================
 */

/*
 * CVAdataStore
 *
 * This routine integrates the forward model starting at the check
 * point ck_mem and stores y and yprime at all intermediate steps.
 *
 * Return values:
 * CV_SUCCESS
 * CV_REIFWD_FAIL
 * CV_FWD_FAIL
 */

static int CVAdataStore(CVodeMem cv_mem, CkpntMem ck_mem)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  realtype t;
  long int i;
  int flag, sign;

  ca_mem = cv_mem->cv_adj_mem;
  dt_mem = ca_mem->dt_mem;

  /* Initialize cv_mem with data from ck_mem */
  flag = CVAckpntGet(cv_mem, ck_mem);
  if (flag != CV_SUCCESS)
    return(CV_REIFWD_FAIL);

  /* Set first structure in dt_mem[0] */
  dt_mem[0]->t = t0_;
  IMstore(cv_mem, dt_mem[0]);

  /* Decide whether TSTOP must be activated */
  if (ca_mem->ca_tstopCVodeFcall) {
    CVodeSetStopTime(cv_mem, ca_mem->ca_tstopCVodeF);
  }

  sign = (tfinal - tinitial > ZERO) ? 1 : -1;


  /* Run CVode to set following structures in dt_mem[i] */
  i = 1;
  do {

    flag = CVode(cv_mem, t1_, ytmp, &t, CV_ONE_STEP);
    if (flag < 0) return(CV_FWD_FAIL);

    dt_mem[i]->t = t;
    IMstore(cv_mem, dt_mem[i]);
    i++;

  } while ( sign*(t1_ - t) > ZERO );


  IMnewData = TRUE;     /* New data is now available    */
  ckpntData = ck_mem;   /* starting at this check point */
  np = i;               /* and we have this many points */

  return(CV_SUCCESS);
}

/*
 * CVAckpntGet
 *
 * This routine prepares CVODES for a hot restart from
 * the check point ck_mem
 */

static int CVAckpntGet(CVodeMem cv_mem, CkpntMem ck_mem) 
{
  int flag, j, is, qmax;

  if (next_ == NULL) {

    /* In this case, we just call the reinitialization routine,
     * but make sure we use the same initial stepsize as on 
     * the first run. */

    CVodeSetInitStep(cv_mem, h0u);

    flag = CVodeReInit(cv_mem, t0_, zn_[0]);
    if (flag != CV_SUCCESS) return(flag);

    if (quadr_) {
      flag = CVodeQuadReInit(cv_mem, znQ_[0]);
      if (flag != CV_SUCCESS) return(flag);
    }

    if (sensi_) {
      flag = CVodeSensReInit(cv_mem, cv_mem->cv_ism, znS_[0]);
      if (flag != CV_SUCCESS) return(flag);
    }

    if (quadr_sensi_) {
      flag = CVodeQuadSensReInit(cv_mem, znQS_[0]);
      if (flag != CV_SUCCESS) return(flag);
    }

  } else {
    
    qmax = cv_mem->cv_qmax;

    /* Copy parameters from check point data structure */

    nst       = nst_;
    tretlast  = tretlast_;
    q         = q_;
    qprime    = qprime_;
    qwait     = qwait_;
    L         = L_;
    gammap    = gammap_;
    h         = h_;
    hprime    = hprime_;
    hscale    = hscale_;
    eta       = eta_;
    etamax    = etamax_;
    tn        = t0_;
    saved_tq5 = saved_tq5_;
    
    /* Copy the arrays from check point data structure */

    for (j=0; j<=q; j++) N_VScale(ONE, zn_[j], zn[j]);
    if ( q < qmax ) N_VScale(ONE, zn_[qmax], zn[qmax]);

    if (quadr_) {
      for (j=0; j<=q; j++) N_VScale(ONE, znQ_[j], znQ[j]);
      if ( q < qmax ) N_VScale(ONE, znQ_[qmax], znQ[qmax]);
    }

    if (sensi_) {
      for (is=0; is<Ns; is++) {
        for (j=0; j<=q; j++) N_VScale(ONE, znS_[j][is], znS[j][is]);
        if ( q < qmax ) N_VScale(ONE, znS_[qmax][is], znS[qmax][is]);
      }
    }

    if (quadr_sensi_) {
      for (is=0; is<Ns; is++) {
        for (j=0; j<=q; j++) N_VScale(ONE, znQS_[j][is], znQS[j][is]);
        if ( q < qmax ) N_VScale(ONE, znQS_[qmax][is], znQS[qmax][is]);
      }
    }

    for (j=0; j<=L_MAX; j++)     tau[j] = tau_[j];
    for (j=0; j<=NUM_TESTS; j++) tq[j] = tq_[j];
    for (j=0; j<=q; j++)         l[j] = l_[j];
    
    /* Force a call to setup */

    forceSetup = TRUE;

  }

  return(CV_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Functions for interpolation
 * -----------------------------------------------------------------
 */

/*
 * CVAfindIndex
 *
 * Finds the index in the array of data point strctures such that
 *     dt_mem[indx-1].t <= t < dt_mem[indx].t
 * If indx is changed from the previous invocation, then newpoint = TRUE
 *
 * If t is beyond the leftmost limit, but close enough, indx=0.
 *
 * Returns CV_SUCCESS if successful and CV_GETY_BADT if unable to
 * find indx (t is too far beyond limits).
 */

static int CVAfindIndex(CVodeMem cv_mem, realtype t, 
                        long int *indx, booleantype *newpoint)
{
  CVadjMem ca_mem;
  static long int ilast;
  DtpntMem *dt_mem;
  int sign;
  booleantype to_left, to_right;

  ca_mem = cv_mem->cv_adj_mem;
  dt_mem = ca_mem->dt_mem;

  *newpoint = FALSE;

  /* Find the direction of integration */
  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* If this is the first time we use new data */
  if (IMnewData) {
    ilast     = np-1;
    *newpoint = TRUE;
    IMnewData   = FALSE;
  }

  /* Search for indx starting from ilast */
  to_left  = ( sign*(t - dt_mem[ilast-1]->t) < ZERO);
  to_right = ( sign*(t - dt_mem[ilast]->t)   > ZERO);

  if ( to_left ) {
    /* look for a new indx to the left */

    *newpoint = TRUE;
    
    *indx = ilast;
    loop {
      if ( *indx == 0 ) break;
      if ( sign*(t - dt_mem[*indx-1]->t) <= ZERO ) (*indx)--;
      else                                         break;
    }

    if ( *indx == 0 )
      ilast = 1;
    else
      ilast = *indx;

    if ( *indx == 0 ) {
      /* t is beyond leftmost limit. Is it too far? */  
      if ( SUNRabs(t - dt_mem[0]->t) > FUZZ_FACTOR * uround ) {
        return(CV_GETY_BADT);
      }
    }

  } else if ( to_right ) {
    /* look for a new indx to the right */

    *newpoint = TRUE;

    *indx = ilast;
    loop {
      if ( sign*(t - dt_mem[*indx]->t) > ZERO) (*indx)++;
      else                                     break;
    }

    ilast = *indx;


  } else {
    /* ilast is still OK */

    *indx = ilast;

  }

  return(CV_SUCCESS);


}

/*
 * CVodeGetAdjY
 *
 * This routine returns the interpolated forward solution at time t.
 * The user must allocate space for y.
 */

int CVodeGetAdjY(void *cvode_mem, realtype t, N_Vector y)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  int flag;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeGetAdjY", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  flag = IMget(cv_mem, t, y, NULL);

  return(flag);
}

/* 
 * -----------------------------------------------------------------
 * Functions specific to cubic Hermite interpolation
 * -----------------------------------------------------------------
 */

/*
 * CVAhermiteMalloc
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points. 
 * This data is then used to interpolate the forward solution 
 * at any other time.
 */

static booleantype CVAhermiteMalloc(CVodeMem cv_mem)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  HermiteDataMem content;
  long int i, ii=0;
  booleantype allocOK;

  allocOK = TRUE;

  ca_mem = cv_mem->cv_adj_mem;

  /* Allocate space for the vectors ytmp and yStmp */

  ytmp = N_VClone(tempv);
  if (ytmp == NULL) {
    return(FALSE);
  }

  if (IMstoreSensi) {
    yStmp = N_VCloneVectorArray(Ns, tempv);
    if (yStmp == NULL) {
      N_VDestroy(ytmp);
      return(FALSE);
    }
  }

  /* Allocate space for the content field of the dt structures */

  dt_mem = ca_mem->dt_mem;

  for (i=0; i<=nsteps; i++) {

    content = NULL;
    content = (HermiteDataMem) malloc(sizeof(struct HermiteDataMemRec));
    if (content == NULL) {
      ii = i;
      allocOK = FALSE;
      break;
    }

    content->y = N_VClone(tempv);
    if (content->y == NULL) {
      free(content); content = NULL;
      ii = i;
      allocOK = FALSE;
      break;
    }

    content->yd = N_VClone(tempv);
    if (content->yd == NULL) {
      N_VDestroy(content->y);
      free(content); content = NULL;
      ii = i;
      allocOK = FALSE;
      break;
    }

    if (IMstoreSensi) {

      content->yS = N_VCloneVectorArray(Ns, tempv);
      if (content->yS == NULL) {
        N_VDestroy(content->y);
        N_VDestroy(content->yd);
        free(content); content = NULL;
        ii = i;
        allocOK = FALSE;
        break;
      }

      content->ySd = N_VCloneVectorArray(Ns, tempv);
      if (content->ySd == NULL) {
        N_VDestroy(content->y);
        N_VDestroy(content->yd);
        N_VDestroyVectorArray(content->yS, Ns);
        free(content); content = NULL;
        ii = i;
        allocOK = FALSE;
        break;
      }
      
    }
    
    dt_mem[i]->content = content;

  } 

  /* If an error occurred, deallocate and return */

  if (!allocOK) {

    N_VDestroy(ytmp);

    if (IMstoreSensi) {
      N_VDestroyVectorArray(yStmp, Ns);
    }

    for (i=0; i<ii; i++) {
      content = (HermiteDataMem) (dt_mem[i]->content);
      N_VDestroy(content->y);
      N_VDestroy(content->yd);
      if (IMstoreSensi) {
        N_VDestroyVectorArray(content->yS, Ns);
        N_VDestroyVectorArray(content->ySd, Ns);
      }
      free(dt_mem[i]->content); dt_mem[i]->content = NULL;
    }

  }

  return(allocOK);
}

/*
 * CVAhermiteFree
 *
 * This routine frees the memory allocated for data storage.
 */

static void CVAhermiteFree(CVodeMem cv_mem)
{  
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  HermiteDataMem content;
  long int i;

  ca_mem = cv_mem->cv_adj_mem;

  N_VDestroy(ytmp);

  if (IMstoreSensi) {
    N_VDestroyVectorArray(yStmp, Ns);
  }

  dt_mem = ca_mem->dt_mem;

  for (i=0; i<=nsteps; i++) {
    content = (HermiteDataMem) (dt_mem[i]->content);
    N_VDestroy(content->y);
    N_VDestroy(content->yd);
    if (IMstoreSensi) {
      N_VDestroyVectorArray(content->yS, Ns);
      N_VDestroyVectorArray(content->ySd, Ns);
    }
    free(dt_mem[i]->content); dt_mem[i]->content = NULL;
  }
}

/*
 * CVAhermiteStorePnt ( -> IMstore )
 *
 * This routine stores a new point (y,yd) in the structure d for use
 * in the cubic Hermite interpolation.
 * Note that the time is already stored.
 */

static int CVAhermiteStorePnt(CVodeMem cv_mem, DtpntMem d)
{
  CVadjMem ca_mem;
  HermiteDataMem content;
  int is, retval;

  ca_mem = cv_mem->cv_adj_mem;

  content = (HermiteDataMem) d->content;

  /* Load solution */

  N_VScale(ONE, zn[0], content->y);
  
  if (IMstoreSensi) {
    for (is=0; is<Ns; is++) 
      N_VScale(ONE, znS[0][is], content->yS[is]);
  }

  /* Load derivative */

  if (nst == 0) {

    retval = f(tn, content->y, content->yd, user_data);

    if (IMstoreSensi) {
      retval = cvSensRhsWrapper(cv_mem, tn, content->y, content->yd,
                                content->yS, content->ySd,
                                cv_mem->cv_tempv, cv_mem->cv_ftemp);
    }

  } else {

    N_VScale(ONE/h, zn[1], content->yd);

    if (IMstoreSensi) {
      for (is=0; is<Ns; is++) 
        N_VScale(ONE/h, znS[1][is], content->ySd[is]);
    }

  }

  return(0);
}

/*
 * CVAhermiteGetY ( -> IMget )
 *
 * This routine uses cubic piece-wise Hermite interpolation for 
 * the forward solution vector. 
 * It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB) but
 * can be directly called by the user through CVodeGetAdjY
 */

static int CVAhermiteGetY(CVodeMem cv_mem, realtype t,
                          N_Vector y, N_Vector *yS)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  HermiteDataMem content0, content1;

  realtype t0, t1, delta;
  realtype factor1, factor2, factor3;

  N_Vector y0, yd0, y1, yd1;
  N_Vector *yS0=NULL, *ySd0=NULL, *yS1, *ySd1;

  int flag, is, NS;
  long int indx;
  booleantype newpoint;

 
  ca_mem = cv_mem->cv_adj_mem;
  dt_mem = ca_mem->dt_mem;
 
  /* Local value of Ns */
 
  NS = IMinterpSensi ? Ns : 0;

  /* Get the index in dt_mem */

  flag = CVAfindIndex(cv_mem, t, &indx, &newpoint);
  if (flag != CV_SUCCESS) return(flag);

  /* If we are beyond the left limit but close enough,
     then return y at the left limit. */

  if (indx == 0) {
    content0 = (HermiteDataMem) (dt_mem[0]->content);
    N_VScale(ONE, content0->y, y);
    for (is=0; is<NS; is++) N_VScale(ONE, content0->yS[is], yS[is]);
    return(CV_SUCCESS);
  }

  /* Extract stuff from the appropriate data points */

  t0 = dt_mem[indx-1]->t;
  t1 = dt_mem[indx]->t;
  delta = t1 - t0;

  content0 = (HermiteDataMem) (dt_mem[indx-1]->content);
  y0  = content0->y;
  yd0 = content0->yd;
  if (IMinterpSensi) {
    yS0  = content0->yS;
    ySd0 = content0->ySd;
  }

  if (newpoint) {
    
    /* Recompute Y0 and Y1 */

    content1 = (HermiteDataMem) (dt_mem[indx]->content);

    y1  = content1->y;
    yd1 = content1->yd;

    N_VLinearSum(ONE, y1, -ONE, y0, Y[0]);
    N_VLinearSum(ONE, yd1,  ONE, yd0, Y[1]);
    N_VLinearSum(delta, Y[1], -TWO, Y[0], Y[1]);
    N_VLinearSum(ONE, Y[0], -delta, yd0, Y[0]);


    yS1  = content1->yS;
    ySd1 = content1->ySd;
      
    for (is=0; is<NS; is++) {
      N_VLinearSum(ONE, yS1[is], -ONE, yS0[is], YS[0][is]);
      N_VLinearSum(ONE, ySd1[is],  ONE, ySd0[is], YS[1][is]);
      N_VLinearSum(delta, YS[1][is], -TWO, YS[0][is], YS[1][is]);
      N_VLinearSum(ONE, YS[0][is], -delta, ySd0[is], YS[0][is]);
    }

  }

  /* Perform the actual interpolation. */

  factor1 = t - t0;

  factor2 = factor1/delta;
  factor2 = factor2*factor2;

  factor3 = factor2*(t-t1)/delta;

  N_VLinearSum(ONE, y0, factor1, yd0, y);
  N_VLinearSum(ONE, y, factor2, Y[0], y);
  N_VLinearSum(ONE, y, factor3, Y[1], y);

  for (is=0; is<NS; is++) {
    N_VLinearSum(ONE, yS0[is], factor1, ySd0[is], yS[is]);
    N_VLinearSum(ONE, yS[is], factor2, YS[0][is], yS[is]);
    N_VLinearSum(ONE, yS[is], factor3, YS[1][is], yS[is]);
  }


  return(CV_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Functions specific to Polynomial interpolation
 * -----------------------------------------------------------------
 */

/*
 * CVApolynomialMalloc
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points. 
 * This data is then used to interpolate the forward solution 
 * at any other time.
 */

static booleantype CVApolynomialMalloc(CVodeMem cv_mem)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  PolynomialDataMem content;
  long int i, ii=0;
  booleantype allocOK;

  allocOK = TRUE;

  ca_mem = cv_mem->cv_adj_mem;

  /* Allocate space for the vectors ytmp and yStmp */

  ytmp = N_VClone(tempv);
  if (ytmp == NULL) {
    return(FALSE);
  }

  if (IMstoreSensi) {
    yStmp = N_VCloneVectorArray(Ns, tempv);
    if (yStmp == NULL) {
      N_VDestroy(ytmp);
      return(FALSE);
    }
  }

  /* Allocate space for the content field of the dt structures */

  dt_mem = ca_mem->dt_mem;

  for (i=0; i<=nsteps; i++) {

    content = NULL;
    content = (PolynomialDataMem) malloc(sizeof(struct PolynomialDataMemRec));
    if (content == NULL) {
      ii = i;
      allocOK = FALSE;
      break;
    }

    content->y = N_VClone(tempv);
    if (content->y == NULL) {
      free(content); content = NULL;
      ii = i;
      allocOK = FALSE;
      break;
    }

    if (IMstoreSensi) {

      content->yS = N_VCloneVectorArray(Ns, tempv);
      if (content->yS == NULL) {
        N_VDestroy(content->y);
        free(content); content = NULL;
        ii = i;
        allocOK = FALSE;
        break;
      }

    }

    dt_mem[i]->content = content;

  } 

  /* If an error occurred, deallocate and return */

  if (!allocOK) {

    N_VDestroy(ytmp);

    if (IMstoreSensi) {
      N_VDestroyVectorArray(yStmp, Ns);
    }

    for (i=0; i<ii; i++) {
      content = (PolynomialDataMem) (dt_mem[i]->content);
      N_VDestroy(content->y);
      if (IMstoreSensi) {
        N_VDestroyVectorArray(content->yS, Ns);
      }
      free(dt_mem[i]->content); dt_mem[i]->content = NULL;
    }

  }

  return(allocOK);

}

/*
 * CVApolynomialFree
 *
 * This routine frees the memeory allocated for data storage.
 */

static void CVApolynomialFree(CVodeMem cv_mem)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  PolynomialDataMem content;
  long int i;

  ca_mem = cv_mem->cv_adj_mem;

  N_VDestroy(ytmp);

  if (IMstoreSensi) {
    N_VDestroyVectorArray(yStmp, Ns);
  }

  dt_mem = ca_mem->dt_mem;

  for (i=0; i<=nsteps; i++) {
    content = (PolynomialDataMem) (dt_mem[i]->content);
    N_VDestroy(content->y);
    if (IMstoreSensi) {
      N_VDestroyVectorArray(content->yS, Ns);
    }
    free(dt_mem[i]->content); dt_mem[i]->content = NULL;
  }
}

/*
 * CVApolynomialStorePnt ( -> IMstore )
 *
 * This routine stores a new point y in the structure d for use
 * in the Polynomial interpolation.
 * Note that the time is already stored.
 */

static int CVApolynomialStorePnt(CVodeMem cv_mem, DtpntMem d)
{
  CVadjMem ca_mem;
  PolynomialDataMem content;
  int is;

  ca_mem = cv_mem->cv_adj_mem;

  content = (PolynomialDataMem) d->content;

  N_VScale(ONE, zn[0], content->y);

  if (IMstoreSensi) {
    for (is=0; is<Ns; is++) 
      N_VScale(ONE, znS[0][is], content->yS[is]);
  }

  content->order = qu;

  return(0);
}

/*
 * CVApolynomialGetY ( -> IMget )
 *
 * This routine uses polynomial interpolation for the forward solution vector. 
 * It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB)) but
 * can be directly called by the user through CVodeGetAdjY.
 */

static int CVApolynomialGetY(CVodeMem cv_mem, realtype t,
                             N_Vector y, N_Vector *yS)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  PolynomialDataMem content;

  int flag, dir, order, i, j, is, NS;
  long int indx, base;
  booleantype newpoint;
  realtype dt, factor;

  ca_mem = cv_mem->cv_adj_mem;
  dt_mem = ca_mem->dt_mem;
  
  /* Local value of Ns */
 
  NS = IMinterpSensi ? Ns : 0;

  /* Get the index in dt_mem */

  flag = CVAfindIndex(cv_mem, t, &indx, &newpoint);
  if (flag != CV_SUCCESS) return(flag);

  /* If we are beyond the left limit but close enough,
     then return y at the left limit. */

  if (indx == 0) {
    content = (PolynomialDataMem) (dt_mem[0]->content);
    N_VScale(ONE, content->y, y);
    for (is=0; is<NS; is++) N_VScale(ONE, content->yS[is], yS[is]);
    return(CV_SUCCESS);
  }

  /* Scaling factor */

  dt = SUNRabs(dt_mem[indx]->t - dt_mem[indx-1]->t);

  /* Find the direction of the forward integration */

  dir = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* Establish the base point depending on the integration direction.
     Modify the base if there are not enough points for the current order */

  if (dir == 1) {
    base = indx;
    content = (PolynomialDataMem) (dt_mem[base]->content);
    order = content->order;
    if(indx < order) base += order-indx;
  } else {
    base = indx-1;
    content = (PolynomialDataMem) (dt_mem[base]->content);
    order = content->order;
    if (np-indx > order) base -= indx+order-np;
  }

  /* Recompute Y (divided differences for Newton polynomial) if needed */

  if (newpoint) {

    /* Store 0-th order DD */
    if (dir == 1) {
      for(j=0;j<=order;j++) {
        T[j] = dt_mem[base-j]->t;
        content = (PolynomialDataMem) (dt_mem[base-j]->content);
        N_VScale(ONE, content->y, Y[j]);
        for (is=0; is<NS; is++) N_VScale(ONE, content->yS[is], YS[j][is]);
      }
    } else {
      for(j=0;j<=order;j++) {
        T[j] = dt_mem[base-1+j]->t;
        content = (PolynomialDataMem) (dt_mem[base-1+j]->content);
        N_VScale(ONE, content->y, Y[j]);
        for (is=0; is<NS; is++) N_VScale(ONE, content->yS[is], YS[j][is]);
      }
    }

    /* Compute higher-order DD */
    for(i=1;i<=order;i++) {
      for(j=order;j>=i;j--) {
        factor = dt/(T[j]-T[j-i]);
        N_VLinearSum(factor, Y[j], -factor, Y[j-1], Y[j]);
        for (is=0; is<NS; is++) N_VLinearSum(factor, YS[j][is], -factor, YS[j-1][is], YS[j][is]);
      }
    }
  }

  /* Perform the actual interpolation using nested multiplications */

  N_VScale(ONE, Y[order], y);
  for (is=0; is<NS; is++) N_VScale(ONE, YS[order][is], yS[is]);
  for (i=order-1; i>=0; i--) {
    factor = (t-T[i])/dt;
    N_VLinearSum(factor, y, ONE, Y[i], y);
    for (is=0; is<NS; is++) N_VLinearSum(factor, yS[is], ONE, YS[i][is], yS[is]);
  }

  return(CV_SUCCESS);

}

/* 
 * =================================================================
 * WRAPPERS FOR ADJOINT SYSTEM
 * =================================================================
 */
/*
 * CVArhs
 *
 * This routine interfaces to the CVRhsFnB (or CVRhsFnBS) routine 
 * provided by the user.
 */

static int CVArhs(realtype t, N_Vector yB, 
                  N_Vector yBdot, void *cvode_mem)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  int flag, retval;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  /* Get forward solution from interpolation */

  if (IMinterpSensi)
    flag = IMget(cv_mem, t, ytmp, yStmp);
  else 
    flag = IMget(cv_mem, t, ytmp, NULL);

  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVODEA", "CVArhs", MSGCV_BAD_TINTERP, t);
    return(-1);
  }

  /* Call the user's RHS function */

  if (cvB_mem->cv_f_withSensi)
    retval = (cvB_mem->cv_fs)(t, ytmp, yStmp, yB, yBdot, cvB_mem->cv_user_data);
  else
    retval = (cvB_mem->cv_f)(t, ytmp, yB, yBdot, cvB_mem->cv_user_data);

  return(retval);
}

/*
 * CVArhsQ
 *
 * This routine interfaces to the CVQuadRhsFnB (or CVQuadRhsFnBS) routine
 * provided by the user.
 */

static int CVArhsQ(realtype t, N_Vector yB, 
                   N_Vector qBdot, void *cvode_mem)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  int flag, retval;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  /* Get forward solution from interpolation */

  if (IMinterpSensi)
    flag = IMget(cv_mem, t, ytmp, yStmp);
  else 
    flag = IMget(cv_mem, t, ytmp, NULL);

  /* Call the user's RHS function */

  if (cvB_mem->cv_fQ_withSensi)
    retval = (cvB_mem->cv_fQs)(t, ytmp, yStmp, yB, qBdot, cvB_mem->cv_user_data);
  else
    retval = (cvB_mem->cv_fQ)(t, ytmp, yB, qBdot, cvB_mem->cv_user_data);

  return(retval);
}
