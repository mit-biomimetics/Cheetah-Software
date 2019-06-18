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
 * This is the implementation file for the IDAA adjoint integrator.
 * -----------------------------------------------------------------
 */

/*=================================================================*/
/*                  Import Header Files                            */
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>

#include "idas_impl.h"
#include <sundials/sundials_math.h>

/*=================================================================*/
/*                  Macros                                         */
/*=================================================================*/

#define loop for(;;)

/*=================================================================*/
/*                 IDAA Private Constants                          */
/*=================================================================*/

#define ZERO        RCONST(0.0)    /* real   0.0 */
#define ONE         RCONST(1.0)    /* real   1.0 */
#define TWO         RCONST(2.0)    /* real   2.0 */
#define HUNDRED     RCONST(100.0)  /* real 100.0 */
#define FUZZ_FACTOR RCONST(1000000.0)  /* fuzz factor for IDAAgetY */


/*=================================================================*/
/*               Private Functions Prototypes                      */
/*=================================================================*/

static CkpntMem IDAAckpntInit(IDAMem IDA_mem);
static CkpntMem IDAAckpntNew(IDAMem IDA_mem);
static void IDAAckpntCopyVectors(IDAMem IDA_mem, CkpntMem ck_mem);
static booleantype IDAAckpntAllocVectors(IDAMem IDA_mem, CkpntMem ck_mem);
static void IDAAckpntDelete(CkpntMem *ck_memPtr);

static void IDAAbckpbDelete(IDABMem *IDAB_memPtr);

static booleantype IDAAdataMalloc(IDAMem IDA_mem);
static void IDAAdataFree(IDAMem IDA_mem);
static int  IDAAdataStore(IDAMem IDA_mem, CkpntMem ck_mem);

static int  IDAAckpntGet(IDAMem IDA_mem, CkpntMem ck_mem); 

static booleantype IDAAhermiteMalloc(IDAMem IDA_mem);
static void        IDAAhermiteFree(IDAMem IDA_mem);
static int         IDAAhermiteStorePnt(IDAMem IDA_mem, DtpntMem d);
static int         IDAAhermiteGetY(IDAMem IDA_mem, realtype t, 
                                   N_Vector yy, N_Vector yp,
                                   N_Vector *yyS, N_Vector *ypS);

static booleantype IDAApolynomialMalloc(IDAMem IDA_mem);
static void        IDAApolynomialFree(IDAMem IDA_mem);
static int         IDAApolynomialStorePnt(IDAMem IDA_mem, DtpntMem d);
static int         IDAApolynomialGetY(IDAMem IDA_mem, realtype t, 
                                      N_Vector yy, N_Vector yp,
                                      N_Vector *yyS, N_Vector *ypS);

static int IDAAfindIndex(IDAMem ida_mem, realtype t, 
                         long int *indx, booleantype *newpoint);                         

static int IDAAres(realtype tt, 
                   N_Vector yyB, N_Vector ypB, 
                   N_Vector resvalB,  void *ida_mem);

static int IDAArhsQ(realtype tt, 
                     N_Vector yyB, N_Vector ypB,
                     N_Vector rrQB, void *ida_mem);

static int IDAAGettnSolutionYp(IDAMem IDA_mem, N_Vector yp);
static int IDAAGettnSolutionYpS(IDAMem IDA_mem, N_Vector *ypS);

extern int IDAGetSolution(void *ida_mem, realtype t, N_Vector yret, N_Vector ypret);

/*=================================================================*/
/*             Readibility Constants                               */
/*=================================================================*/

/* IDAADJ memory block */
#define tinitial     (IDAADJ_mem->ia_tinitial)
#define tfinal       (IDAADJ_mem->ia_tfinal)
#define nckpnts      (IDAADJ_mem->ia_nckpnts)
#define nbckpbs      (IDAADJ_mem->ia_nbckpbs)
#define nsteps       (IDAADJ_mem->ia_nsteps)
#define ckpntData    (IDAADJ_mem->ia_ckpntData)
#define newData      (IDAADJ_mem->ia_newData)
#define np           (IDAADJ_mem->ia_np)
#define dt           (IDAADJ_mem->ia_dt)
#define yyTmp        (IDAADJ_mem->ia_yyTmp)
#define ypTmp        (IDAADJ_mem->ia_ypTmp)
#define yySTmp       (IDAADJ_mem->ia_yySTmp)
#define ypSTmp       (IDAADJ_mem->ia_ypSTmp)
#define res_B        (IDAADJ_mem->ia_resB)
#define djac_B       (IDAADJ_mem->ia_djacB)
#define bjac_B       (IDAADJ_mem->ia_bjacB)
#define pset_B       (IDAADJ_mem->ia_psetB)
#define psolve_B     (IDAADJ_mem->ia_psolveB)
#define jtimes_B     (IDAADJ_mem->ia_jtimesB)
#define jdata_B      (IDAADJ_mem->ia_jdataB)
#define pdata_B      (IDAADJ_mem->ia_pdataB)
#define rhsQ_B       (IDAADJ_mem->ia_rhsQB)

#define Y            (IDAADJ_mem->ia_Y)
#define YS           (IDAADJ_mem->ia_YS)
#define T            (IDAADJ_mem->ia_T)
#define mallocDone   (IDAADJ_mem->ia_mallocDone)

#define interpSensi  (IDAADJ_mem->ia_interpSensi)
#define storeSensi   (IDAADJ_mem->ia_storeSensi)
#define noInterp     (IDAADJ_mem->ia_noInterp)

/* Forward IDAS memory block */
#define uround     (IDA_mem->ida_uround)
#define res        (IDA_mem->ida_res)
#define itol       (IDA_mem->ida_itol)
#define reltol     (IDA_mem->ida_reltol)
#define abstol     (IDA_mem->ida_abstol)
#define user_data  (IDA_mem->ida_user_data)

#define forceSetup (IDA_mem->ida_forceSetup)
#define h0u        (IDA_mem->ida_h0u)

#define phi        (IDA_mem->ida_phi)
#define psi        (IDA_mem->ida_psi)
#define alpha      (IDA_mem->ida_alpha)
#define beta       (IDA_mem->ida_beta)
#define sigma      (IDA_mem->ida_sigma)
#define gamma      (IDA_mem->ida_gamma)
#define tn         (IDA_mem->ida_tn)
#define kk         (IDA_mem->ida_kk)
#define nst        (IDA_mem->ida_nst)
#define tretlast   (IDA_mem->ida_tretlast)
#define kk         (IDA_mem->ida_kk)
#define kused      (IDA_mem->ida_kused)
#define knew       (IDA_mem->ida_knew)
#define maxord     (IDA_mem->ida_maxord)
#define phase      (IDA_mem->ida_phase)
#define ns         (IDA_mem->ida_ns)
#define hh         (IDA_mem->ida_hh)
#define hused      (IDA_mem->ida_hused)
#define rr         (IDA_mem->ida_rr)
#define cj         (IDA_mem->ida_cj)
#define cjlast     (IDA_mem->ida_cjlast)
#define cjold      (IDA_mem->ida_cjold)
#define cjratio    (IDA_mem->ida_cjratio) 
#define ss         (IDA_mem->ida_ss)
#define ssS        (IDA_mem->ida_ssS)

#define tempv      (IDA_mem->ida_tempv1)

#define sensi      (IDA_mem->ida_sensi)
#define Ns         (IDA_mem->ida_Ns)
#define phiS       (IDA_mem->ida_phiS)

#define quadr      (IDA_mem->ida_quadr)
#define errconQ    (IDA_mem->ida_errconQ)
#define phiQ       (IDA_mem->ida_phiQ)
#define rhsQ       (IDA_mem->ida_rhsQ)

#define quadr_sensi (IDA_mem->ida_quadr_sensi)
#define errconQS   (IDA_mem->ida_errconQS)
#define phiQS      (IDA_mem->ida_phiQS)

#define tempvQ     (IDA_mem->ida_eeQ)

/* Checkpoint memory block */

#define t0_        (ck_mem->ck_t0)
#define t1_        (ck_mem->ck_t1)
#define phi_       (ck_mem->ck_phi)
#define phiQ_      (ck_mem->ck_phiQ)
#define psi_       (ck_mem->ck_psi)
#define alpha_     (ck_mem->ck_alpha)
#define beta_      (ck_mem->ck_beta)
#define sigma_     (ck_mem->ck_sigma)
#define gamma_     (ck_mem->ck_gamma)
#define nst_       (ck_mem->ck_nst)
#define tretlast_  (ck_mem->ck_tretlast)
#define kk_        (ck_mem->ck_kk)
#define kused_     (ck_mem->ck_kused)
#define knew_      (ck_mem->ck_knew)
#define phase_     (ck_mem->ck_phase)
#define ns_        (ck_mem->ck_ns)
#define hh_        (ck_mem->ck_hh)
#define hused_     (ck_mem->ck_hused)
#define rr_        (ck_mem->ck_rr)
#define cj_        (ck_mem->ck_cj)
#define cjlast_    (ck_mem->ck_cjlast)
#define cjold_     (ck_mem->ck_cjold)
#define cjratio_   (ck_mem->ck_cjratio)
#define ss_        (ck_mem->ck_ss)
#define ssS_       (ck_mem->ck_ssS)
#define next_      (ck_mem->ck_next)
#define phi_alloc_ (ck_mem->ck_phi_alloc)

#define sensi_     (ck_mem->ck_sensi)
#define Ns_        (ck_mem->ck_Ns)
#define phiS_      (ck_mem->ck_phiS)

#define quadr_     (ck_mem->ck_quadr)
#define phiQS_     (ck_mem->ck_phiQS)

#define quadr_sensi_ (ck_mem->ck_quadr_sensi)

/*=================================================================*/
/*                  Exported Functions                             */
/*=================================================================*/

/*
 * IDAAdjInit 
 *
 * This routine allocates space for the global IDAA memory
 * structure.
 */


int IDAAdjInit(void *ida_mem, long int steps, int interp)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;

  /* Check arguments */

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAAdjInit", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  if (steps <= 0) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAAdjInit", MSGAM_BAD_STEPS);
    return(IDA_ILL_INPUT);
  }

  if ( (interp != IDA_HERMITE) && (interp != IDA_POLYNOMIAL) ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAAdjInit", MSGAM_BAD_INTERP);
    return(IDA_ILL_INPUT);
  } 

  /* Allocate memory block for IDAadjMem. */
  IDAADJ_mem = (IDAadjMem) malloc(sizeof(struct IDAadjMemRec));
  if (IDAADJ_mem == NULL) {
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDAA", "IDAAdjInit", MSGAM_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  /* Attach IDAS memory for forward runs */
  IDA_mem->ida_adj_mem = IDAADJ_mem;

  /* Initialization of check points. */
  IDAADJ_mem->ck_mem = NULL;
  IDAADJ_mem->ia_nckpnts = 0;
  IDAADJ_mem->ia_ckpntData = NULL;

  /* Initialize wrapper function workspace to NULL for safe deletion if unused */
  IDAADJ_mem->ia_yyTmp = NULL;
  IDAADJ_mem->ia_ypTmp = NULL;
  IDAADJ_mem->ia_yySTmp = NULL;
  IDAADJ_mem->ia_ypSTmp = NULL;

  /* Initialization of interpolation data. */
  IDAADJ_mem->ia_interpType = interp;
  IDAADJ_mem->ia_nsteps = steps;

  /* Allocate space for the array of Data Point structures. */
  if (IDAAdataMalloc(IDA_mem) == FALSE) {
    free(IDAADJ_mem); IDAADJ_mem = NULL;
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDAA", "IDAAdjInit", MSGAM_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  /* Attach functions for the appropriate interpolation module */
  switch(interp) {

  case IDA_HERMITE:
    IDAADJ_mem->ia_malloc    = IDAAhermiteMalloc;
    IDAADJ_mem->ia_free      = IDAAhermiteFree;
    IDAADJ_mem->ia_getY      = IDAAhermiteGetY;
    IDAADJ_mem->ia_storePnt  = IDAAhermiteStorePnt;
    break;
    
    case IDA_POLYNOMIAL:
    
    IDAADJ_mem->ia_malloc    = IDAApolynomialMalloc;
    IDAADJ_mem->ia_free      = IDAApolynomialFree;
    IDAADJ_mem->ia_getY      = IDAApolynomialGetY;
    IDAADJ_mem->ia_storePnt  = IDAApolynomialStorePnt;
    break;
  }

 /* The interpolation module has not been initialized yet */
  IDAADJ_mem->ia_mallocDone = FALSE;

  /* By default we will store but not interpolate sensitivities
   *  - storeSensi will be set in IDASolveF to FALSE if FSA is not enabled
   *    or if the user forced this through IDAAdjSetNoSensi 
   *  - interpSensi will be set in IDASolveB to TRUE if storeSensi is TRUE 
   *    and if at least one backward problem requires sensitivities 
   *  - noInterp will be set in IDACalcICB to TRUE before the call to
   *    IDACalcIC and FALSE after.*/

  IDAADJ_mem->ia_storeSensi  = TRUE;
  IDAADJ_mem->ia_interpSensi = FALSE;
  IDAADJ_mem->ia_noInterp    = FALSE;

  /* Initialize backward problems. */
  IDAADJ_mem->IDAB_mem = NULL;
  IDAADJ_mem->ia_bckpbCrt = NULL;
  IDAADJ_mem->ia_nbckpbs = 0;

  /* Flags for tracking the first calls to IDASolveF and IDASolveF. */
  IDAADJ_mem->ia_firstIDAFcall = TRUE;
  IDAADJ_mem->ia_tstopIDAFcall = FALSE;
  IDAADJ_mem->ia_firstIDABcall = TRUE;

  /* Adjoint module initialized and allocated. */
  IDA_mem->ida_adj = TRUE;
  IDA_mem->ida_adjMallocDone = TRUE;

  return(IDA_SUCCESS);
} 

/*
 * IDAAdjReInit
 *
 * IDAAdjReInit reinitializes the IDAS memory structure for ASA
 */

int IDAAdjReInit(void *ida_mem)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;

  /* Check arguments */

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAAdjReInit", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Was ASA previously initialized? */
  if(IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAAdjReInit",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }

  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Free all stored  checkpoints. */
  while (IDAADJ_mem->ck_mem != NULL) 
      IDAAckpntDelete(&(IDAADJ_mem->ck_mem));

  IDAADJ_mem->ck_mem = NULL;
  IDAADJ_mem->ia_nckpnts = 0;
  IDAADJ_mem->ia_ckpntData = NULL;

  /* Flags for tracking the first calls to IDASolveF and IDASolveF. */
  IDAADJ_mem->ia_firstIDAFcall = TRUE;
  IDAADJ_mem->ia_tstopIDAFcall = FALSE;
  IDAADJ_mem->ia_firstIDABcall = TRUE;

  return(IDA_SUCCESS);
} 

/*
 * IDAAdjFree
 *
 * IDAAdjFree routine frees the memory allocated by IDAAdjInit.
*/


void IDAAdjFree(void *ida_mem)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;

  if (ida_mem == NULL) return;
  IDA_mem = (IDAMem) ida_mem;

  if(IDA_mem->ida_adjMallocDone) {

    /* Data for adjoint. */
    IDAADJ_mem = IDA_mem->ida_adj_mem;
    
    /* Delete check points one by one */
    while (IDAADJ_mem->ck_mem != NULL) {
      IDAAckpntDelete(&(IDAADJ_mem->ck_mem));
    }

    IDAAdataFree(IDA_mem);

    /* Free all backward problems. */
    while (IDAADJ_mem->IDAB_mem != NULL)
      IDAAbckpbDelete( &(IDAADJ_mem->IDAB_mem) );

    /* Free IDAA memory. */
    free(IDAADJ_mem);

    IDA_mem->ida_adj_mem = NULL;
  }
}

/* 
 * =================================================================
 * PRIVATE FUNCTIONS FOR BACKWARD PROBLEMS
 * =================================================================
 */

static void IDAAbckpbDelete(IDABMem *IDAB_memPtr)
{
  IDABMem IDAB_mem = (*IDAB_memPtr);
  void * ida_mem;

  if (IDAB_mem == NULL) return;

  /* Move head to the next element in list. */
  *IDAB_memPtr = IDAB_mem->ida_next;

  /* IDAB_mem is going to be deallocated. */

  /* Free IDAS memory for this backward problem. */
  ida_mem = (void *)IDAB_mem->IDA_mem;
  IDAFree(&ida_mem);

  /* Free linear solver memory. */
  if (IDAB_mem->ida_lfree != NULL) IDAB_mem->ida_lfree(IDAB_mem);

  /* Free preconditioner memory. */
  if (IDAB_mem->ida_pfree != NULL) IDAB_mem->ida_pfree(IDAB_mem);

  /* Free any workspace vectors. */
  N_VDestroy(IDAB_mem->ida_yy);
  N_VDestroy(IDAB_mem->ida_yp);

  /* Free the node itself. */
  free(IDAB_mem);
  IDAB_mem = NULL;
}

/*=================================================================*/
/*                    Wrappers for IDAA                            */
/*=================================================================*/

/*
 *                      IDASolveF 
 *
 * This routine integrates to tout and returns solution into yout.
 * In the same time, it stores check point data every 'steps' steps. 
 *  
 * IDASolveF can be called repeatedly by the user. The last tout
 *  will be used as the starting time for the backward integration.
 * 
 *  ncheckPtr points to the number of check points stored so far.
*/

int IDASolveF(void *ida_mem, realtype tout, realtype *tret,
              N_Vector yret, N_Vector ypret, int itask, int *ncheckPtr)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  CkpntMem tmp;
  DtpntMem *dt_mem;
  int flag, i;
  booleantype iret, allocOK;

  /* Is the mem OK? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDASolveF", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized ? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDASolveF",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check for yret != NULL */
  if (yret == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASolveF", MSG_YRET_NULL);
    return(IDA_ILL_INPUT);
  }
  
  /* Check for ypret != NULL */
  if (ypret == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASolveF", MSG_YPRET_NULL);
    return(IDA_ILL_INPUT);
  }
  /* Check for tret != NULL */
  if (tret == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASolveF", MSG_TRET_NULL);
    return(IDA_ILL_INPUT);
  }
  
  /* Check for valid itask */
  if ( (itask != IDA_NORMAL) && (itask != IDA_ONE_STEP) ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASolveF", MSG_BAD_ITASK);
    return(IDA_ILL_INPUT);
  }
  
  /* All memory checks done, proceed ... */
  
  dt_mem = IDAADJ_mem->dt_mem;

  /* If tstop is enabled, store some info */
  if (IDA_mem->ida_tstopset) {
    IDAADJ_mem->ia_tstopIDAFcall = TRUE;
    IDAADJ_mem->ia_tstopIDAF = IDA_mem->ida_tstop;
  }
  
  /* We will call IDASolve in IDA_ONE_STEP mode, regardless
     of what itask is, so flag if we need to return */
  if (itask == IDA_ONE_STEP) iret = TRUE;
  else                       iret = FALSE;

  /* On the first step:
   *   - set tinitial
   *   - initialize list of check points
   *   - if needed, initialize the interpolation module
   *   - load dt_mem[0]
   * On subsequent steps, test if taking a new step is necessary. 
   */
  if ( IDAADJ_mem->ia_firstIDAFcall ) {
    
    tinitial = tn;
    IDAADJ_mem->ck_mem = IDAAckpntInit(IDA_mem);
    if (IDAADJ_mem->ck_mem == NULL) {
      IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDAA", "IDASolveF", MSG_MEM_FAIL);
      return(IDA_MEM_FAIL);
    }

    if (!mallocDone) {
      /* Do we need to store sensitivities? */
      if (!sensi) storeSensi = FALSE;

      /* Allocate space for interpolation data */
      allocOK = IDAADJ_mem->ia_malloc(IDA_mem);
      if (!allocOK) {
        IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDAA", "IDASolveF", MSG_MEM_FAIL);
        return(IDA_MEM_FAIL);
      }

      /* Rename phi and, if needed, phiS for use in interpolation */
      for (i=0;i<MXORDP1;i++) Y[i] = phi[i];
      if (storeSensi) {
        for (i=0;i<MXORDP1;i++) YS[i] = phiS[i];
      }

      mallocDone = TRUE;
    }

    dt_mem[0]->t = IDAADJ_mem->ck_mem->ck_t0;
    IDAADJ_mem->ia_storePnt(IDA_mem, dt_mem[0]);

    IDAADJ_mem->ia_firstIDAFcall = FALSE;

  } else if ( (tn-tout)*hh >= ZERO ) {

    /* If tout was passed, return interpolated solution. 
       No changes to ck_mem or dt_mem are needed. */
    *tret = tout;
    flag = IDAGetSolution(IDA_mem, tout, yret, ypret);
    *ncheckPtr = nckpnts;
    newData = TRUE;
    ckpntData = IDAADJ_mem->ck_mem;
    np = nst % nsteps + 1;

    return(flag);
  }
  /* Integrate to tout while loading check points */
  loop {

    /* Perform one step of the integration */

    flag = IDASolve(IDA_mem, tout, tret, yret, ypret, IDA_ONE_STEP);

    if (flag < 0) break;

    /* Test if a new check point is needed */

    if ( nst % nsteps == 0 ) {

      IDAADJ_mem->ck_mem->ck_t1 = *tret;

      /* Create a new check point, load it, and append it to the list */
      tmp = IDAAckpntNew(IDA_mem);
      if (tmp == NULL) {
        flag = IDA_MEM_FAIL;
        break;
      }

      tmp->ck_next = IDAADJ_mem->ck_mem;
      IDAADJ_mem->ck_mem = tmp;
      nckpnts++;
      
      forceSetup = TRUE;
      
      /* Reset i=0 and load dt_mem[0] */
      dt_mem[0]->t = IDAADJ_mem->ck_mem->ck_t0;
      IDAADJ_mem->ia_storePnt(IDA_mem, dt_mem[0]);

    } else {
      
      /* Load next point in dt_mem */
      dt_mem[nst%nsteps]->t = *tret;
      IDAADJ_mem->ia_storePnt(IDA_mem, dt_mem[nst%nsteps]);
    }

    /* Set t1 field of the current ckeck point structure
       for the case in which there will be no future
       check points */
    IDAADJ_mem->ck_mem->ck_t1 = *tret;

    /* tfinal is now set to *t */
    tfinal = *tret;

    /* In IDA_ONE_STEP mode break from loop */
    if (itask == IDA_ONE_STEP) break;

    /* Return if tout reached */
    if ( (*tret - tout)*hh >= ZERO ) {
      *tret = tout;
      IDAGetSolution(IDA_mem, tout, yret, ypret);
      /* Reset tretlast in IDA_mem so that IDAGetQuad and IDAGetSens 
       * evaluate quadratures and/or sensitivities at the proper time */
      IDA_mem->ida_tretlast = tout;
      break;
    }    
  }

  /* Get ncheck from IDAADJ_mem */ 
  *ncheckPtr = nckpnts;

  /* Data is available for the last interval */
  newData = TRUE;
  ckpntData = IDAADJ_mem->ck_mem;
  np = nst % nsteps + 1;

  return(flag);
}




/* 
 * =================================================================
 * FUNCTIONS FOR BACKWARD PROBLEMS
 * =================================================================
 */

int IDACreateB(void *ida_mem, int *which)
{
  IDAMem IDA_mem;
  void* ida_memB;
  IDABMem new_IDAB_mem;
  IDAadjMem IDAADJ_mem;
  
  /* Is the mem OK? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDACreateB", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized ? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDACreateB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Allocate a new IDABMem struct. */
  new_IDAB_mem = (IDABMem) malloc( sizeof( struct IDABMemRec ) );
  if (new_IDAB_mem == NULL) {
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDAA", "IDACreateB",  MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }
  
  /* Allocate the IDAMem struct needed by this backward problem. */
  ida_memB = IDACreate();
  if (ida_memB == NULL) {
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDAA", "IDACreateB",  MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  /* Save ida_mem in ida_memB as user data. */
  IDASetUserData(ida_memB, ida_mem);
  
  /* Set same error output and handler for ida_memB. */
  IDASetErrHandlerFn(ida_memB, IDA_mem->ida_ehfun, IDA_mem->ida_eh_data);
  IDASetErrFile(ida_memB, IDA_mem->ida_errfp);

  /* Initialize fields in the IDABMem struct. */
  new_IDAB_mem->ida_index   = IDAADJ_mem->ia_nbckpbs;
  new_IDAB_mem->IDA_mem     = (IDAMem) ida_memB;

  new_IDAB_mem->ida_res      = NULL;
  new_IDAB_mem->ida_resS     = NULL;
  new_IDAB_mem->ida_rhsQ     = NULL;
  new_IDAB_mem->ida_rhsQS    = NULL;


  new_IDAB_mem->ida_user_data = NULL;

  new_IDAB_mem->ida_lmem     = NULL;
  new_IDAB_mem->ida_lfree    = NULL;
  new_IDAB_mem->ida_pmem     = NULL;
  new_IDAB_mem->ida_pfree    = NULL;

  new_IDAB_mem->ida_yy       = NULL;
  new_IDAB_mem->ida_yp       = NULL;

  new_IDAB_mem->ida_res_withSensi = FALSE;
  new_IDAB_mem->ida_rhsQ_withSensi = FALSE;
  
  /* Attach the new object to the beginning of the linked list IDAADJ_mem->IDAB_mem. */
  new_IDAB_mem->ida_next = IDAADJ_mem->IDAB_mem;
  IDAADJ_mem->IDAB_mem = new_IDAB_mem;

  /* Return the assigned index. This id is used as identificator and has to be passed 
     to IDAInitB and other ***B functions that set the optional inputs for  this 
     backward problem. */
  *which = IDAADJ_mem->ia_nbckpbs;

  /*Increase the counter of the backward problems stored. */
  IDAADJ_mem->ia_nbckpbs++;

  return(IDA_SUCCESS);

}

int IDAInitB(void *ida_mem, int which, IDAResFnB resB,
             realtype tB0, N_Vector yyB0, N_Vector ypB0)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void * ida_memB;
  int flag;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAInitB", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized ? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAInitB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the initial time for this backward problem against the adjoint data. */
  if ( (tB0 < tinitial) || (tB0 > tfinal) ) {
    IDAProcessError(IDA_mem, IDA_BAD_TB0, "IDAA", "IDAInitB", MSGAM_BAD_TB0);
    return(IDA_BAD_TB0);
  }

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAInitB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Get the IDAMem corresponding to this backward problem. */
  ida_memB = (void*) IDAB_mem->IDA_mem;

  /* Call the IDAInit for this backward problem. */
  flag = IDAInit(ida_memB, IDAAres, tB0, yyB0, ypB0);
  if (IDA_SUCCESS != flag) return(flag);

  /* Copy residual function in IDAB_mem. */
  IDAB_mem->ida_res = resB;
  IDAB_mem->ida_res_withSensi = FALSE;

  /* Initialized the initial time field. */
  IDAB_mem->ida_t0 = tB0;

  /* Allocate and initialize space workspace vectors. */
  IDAB_mem->ida_yy = N_VClone(yyB0);
  IDAB_mem->ida_yp = N_VClone(yyB0);
  N_VScale(ONE, yyB0, IDAB_mem->ida_yy);
  N_VScale(ONE, ypB0, IDAB_mem->ida_yp);

  return(flag);

}

int IDAInitBS(void *ida_mem, int which, IDAResFnBS resS,
                realtype tB0, N_Vector yyB0, N_Vector ypB0)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void * ida_memB;
  int flag;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAInitBS", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized ? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAInitBS",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the initial time for this backward problem against the adjoint data. */
  if ( (tB0 < tinitial) || (tB0 > tfinal) ) {
    IDAProcessError(IDA_mem, IDA_BAD_TB0, "IDAA", "IDAInitBS", MSGAM_BAD_TB0);
    return(IDA_BAD_TB0);
  }

  /* Were sensitivities active during the forward integration? */
  if (!storeSensi) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAInitBS", MSGAM_BAD_SENSI);
    return(IDA_ILL_INPUT);
  }

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAInitBS", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Get the IDAMem corresponding to this backward problem. */
  ida_memB = (void*) IDAB_mem->IDA_mem;
  
  /* Allocate and set the IDAS object */
  flag = IDAInit(ida_memB, IDAAres, tB0, yyB0, ypB0);

  if (flag != IDA_SUCCESS) return(flag);

  /* Copy residual function pointer in IDAB_mem. */
  IDAB_mem->ida_res_withSensi = TRUE;
  IDAB_mem->ida_resS = resS;

  /* Allocate space and initialize the yy and yp vectors. */
  IDAB_mem->ida_t0 = tB0;
  IDAB_mem->ida_yy = N_VClone(yyB0);
  IDAB_mem->ida_yp = N_VClone(ypB0);
  N_VScale(ONE, yyB0, IDAB_mem->ida_yy);
  N_VScale(ONE, ypB0, IDAB_mem->ida_yp);

  return(IDA_SUCCESS);
}


int IDAReInitB(void *ida_mem, int which,
               realtype tB0, N_Vector yyB0, N_Vector ypB0)
{

  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void * ida_memB;
  int flag;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAReInitB", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized ? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAReInitB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the initial time for this backward problem against the adjoint data. */
  if ( (tB0 < tinitial) || (tB0 > tfinal) ) {
    IDAProcessError(IDA_mem, IDA_BAD_TB0, "IDAA", "IDAReInitB", MSGAM_BAD_TB0);
    return(IDA_BAD_TB0);
  }

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAReInitB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Get the IDAMem corresponding to this backward problem. */
  ida_memB = (void*) IDAB_mem->IDA_mem;


  /* Call the IDAReInit for this backward problem. */
  flag = IDAReInit(ida_memB, tB0, yyB0, ypB0);
  return(flag);
}

int IDASStolerancesB(void *ida_mem, int which, 
                     realtype relTolB, realtype absTolB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDASStolerancesB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDASStolerancesB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASStolerancesB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Get the IDAMem corresponding to this backward problem. */
  ida_memB = (void*) IDAB_mem->IDA_mem;

  /* Set tolerances and return. */
  return IDASStolerances(ida_memB, relTolB, absTolB);
  
}
int IDASVtolerancesB(void *ida_mem, int which, 
                     realtype relTolB, N_Vector absTolB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDASVtolerancesB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDASVtolerancesB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASVtolerancesB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Get the IDAMem corresponding to this backward problem. */
  ida_memB = (void*) IDAB_mem->IDA_mem;

  /* Set tolerances and return. */
  return IDASVtolerances(ida_memB, relTolB, absTolB);
}

int IDAQuadSStolerancesB(void *ida_mem, int which,
                         realtype reltolQB, realtype abstolQB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAQuadSStolerancesB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAQuadSStolerancesB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAQuadSStolerancesB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void *) IDAB_mem->IDA_mem;
  
  return IDAQuadSStolerances(ida_memB, reltolQB, abstolQB);
}


int IDAQuadSVtolerancesB(void *ida_mem, int which,
                         realtype reltolQB, N_Vector abstolQB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAQuadSVtolerancesB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAQuadSVtolerancesB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAQuadSVtolerancesB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void *) IDAB_mem->IDA_mem;
  
  return IDAQuadSVtolerances(ida_memB, reltolQB, abstolQB);
}


int IDAQuadInitB(void *ida_mem, int which, IDAQuadRhsFnB rhsQB, N_Vector yQB0)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  int flag;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAQuadInitB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAQuadInitB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAQuadInitB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void *) IDAB_mem->IDA_mem;

  flag = IDAQuadInit(ida_memB, IDAArhsQ, yQB0);
  if (IDA_SUCCESS != flag) return flag;

  IDAB_mem->ida_rhsQ_withSensi = FALSE;
  IDAB_mem->ida_rhsQ = rhsQB;

  return(flag);
}


int IDAQuadInitBS(void *ida_mem, int which, 
                  IDAQuadRhsFnBS rhsQS, N_Vector yQB0)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  void * ida_memB;
  int flag;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAQuadInitBS", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized ? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAQuadInitBS",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAQuadInitBS", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Get the IDAMem corresponding to this backward problem. */
  ida_memB = (void*) IDAB_mem->IDA_mem;
  
  /* Allocate and set the IDAS object */
  flag = IDAQuadInit(ida_memB, IDAArhsQ, yQB0);

  if (flag != IDA_SUCCESS) return(flag);

  /* Copy RHS function pointer in IDAB_mem and enable quad sensitivities. */
  IDAB_mem->ida_rhsQ_withSensi = TRUE;
  IDAB_mem->ida_rhsQS = rhsQS;

  return(IDA_SUCCESS);
}


int IDAQuadReInitB(void *ida_mem, int which, N_Vector yQB0)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAQuadInitB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAQuadInitB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAQuadInitB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void *) IDAB_mem->IDA_mem;

  return IDAQuadReInit(ida_mem, yQB0);
}


/*
 * ----------------------------------------------------------------
 * Function : IDACalcICB                                         
 * ----------------------------------------------------------------
 * IDACalcIC calculates corrected initial conditions for a DAE  
 * backward system (index-one in semi-implicit form).
 * It uses Newton iteration combined with a Linesearch algorithm. 
 * Calling IDACalcICB is optional. It is only necessary when the   
 * initial conditions do not solve the given system.  I.e., if    
 * yB0 and ypB0 are known to satisfy the backward problem, then       
 * a call to IDACalcIC is NOT necessary (for index-one problems). 
*/

int IDACalcICB(void *ida_mem, int which, realtype tout1, 
               N_Vector yy0, N_Vector yp0)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  int flag;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDACalcICB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDACalcICB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDACalcICB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void *) IDAB_mem->IDA_mem;

  /* The wrapper for user supplied res function requires ia_bckpbCrt from
     IDAAdjMem to be set to curent problem. */
  IDAADJ_mem->ia_bckpbCrt = IDAB_mem;

  /* Save (y, y') in yyTmp and ypTmp for use in the res wrapper.*/
  /* yyTmp and ypTmp workspaces are safe to use if IDAADataStore is not called.*/
  N_VScale(ONE, yy0, yyTmp);
  N_VScale(ONE, yp0, ypTmp);
  
  /* Set noInterp flag to true, so IDAARes will use user provided values for
     y and y' and will not call the interpolation routine(s). */
  noInterp = TRUE;
  
  flag = IDACalcIC(ida_memB, IDA_YA_YDP_INIT, tout1);

  /* Set interpolation on in IDAARes. */
  noInterp = FALSE;

  return(flag);
}

/*
 * ----------------------------------------------------------------
 * Function : IDACalcICBS                                        
 * ----------------------------------------------------------------
 * IDACalcIC calculates corrected initial conditions for a DAE  
 * backward system (index-one in semi-implicit form) that also 
 * dependes on the sensivities.
 *
 * It calls IDACalcIC for the 'which' backward problem.
*/

int IDACalcICBS(void *ida_mem, int which, realtype tout1, 
               N_Vector yy0, N_Vector yp0, 
               N_Vector *yyS0, N_Vector *ypS0)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  int flag, is;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDACalcICBS", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDACalcICBS",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Were sensitivities active during the forward integration? */
  if (!storeSensi) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDACalcICBS", MSGAM_BAD_SENSI);
    return(IDA_ILL_INPUT);
  }

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDACalcICBS", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void *) IDAB_mem->IDA_mem;

  /* Was InitBS called for this problem? */
  if (!IDAB_mem->ida_res_withSensi) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDACalcICBS", MSGAM_NO_INITBS);
    return(IDA_ILL_INPUT);    
  }

  /* The wrapper for user supplied res function requires ia_bckpbCrt from
     IDAAdjMem to be set to curent problem. */
  IDAADJ_mem->ia_bckpbCrt = IDAB_mem;

  /* Save (y, y') and (y_p, y'_p) in yyTmp, ypTmp and yySTmp, ypSTmp.The wrapper 
     for residual will use these values instead of calling interpolation routine.*/

  /* The four workspaces variables are safe to use if IDAADataStore is not called.*/
  N_VScale(ONE, yy0, yyTmp);
  N_VScale(ONE, yp0, ypTmp);

  for (is=0; is<Ns; is++) {
    N_VScale(ONE, yyS0[is], yySTmp[is]);
    N_VScale(ONE, ypS0[is], ypSTmp[is]);
  }
  
  /* Set noInterp flag to true, so IDAARes will use user provided values for
     y and y' and will not call the interpolation routine(s). */
  noInterp = TRUE;
  
  flag = IDACalcIC(ida_memB, IDA_YA_YDP_INIT, tout1);

  /* Set interpolation on in IDAARes. */
  noInterp = FALSE;

  return(flag);
}


/*
 * IDASolveB
 *
 * This routine performs the backward integration from tB0 
 * to tinitial through a sequence of forward-backward runs in
 * between consecutive check points. It returns the values of
 * the adjoint variables and any existing quadrature variables
 * at tinitial.
 *
 * On a successful return, IDASolveB returns IDA_SUCCESS.
 *
 * NOTE that IDASolveB DOES NOT return the solution for the 
 * backward problem(s). Use IDAGetB to extract the solution 
 * for any given backward problem.
 *
 * If there are multiple backward problems and multiple check points,
 * IDASolveB may not succeed in getting all problems to take one step
 * when called in ONE_STEP mode.
 */

int IDASolveB(void *ida_mem, realtype tBout, int itaskB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  CkpntMem ck_mem;
  IDABMem IDAB_mem, tmp_IDAB_mem;
  int flag=0, sign;
  realtype tfuzz, tBret, tBn;
  booleantype gotCkpnt, reachedTBout, isActive;

  /* Is the mem OK? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDASolveB", MSGAM_NULL_IDAMEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized ? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDASolveB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  if ( nbckpbs == 0 ) {
    IDAProcessError(IDA_mem, IDA_NO_BCK, "IDAA", "IDASolveB", MSGAM_NO_BCK);
    return(IDA_NO_BCK);
  }
  IDAB_mem = IDAADJ_mem->IDAB_mem;

  /* Check whether IDASolveF has been called */
  if ( IDAADJ_mem->ia_firstIDAFcall ) {
    IDAProcessError(IDA_mem, IDA_NO_FWD, "IDAA", "IDASolveB", MSGAM_NO_FWD);
    return(IDA_NO_FWD);
  }
  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* If this is the first call, loop over all backward problems and
   *   - check that tB0 is valid
   *   - check that tBout is ahead of tB0 in the backward direction
   *   - check whether we need to interpolate forward sensitivities
   */
  if (IDAADJ_mem->ia_firstIDABcall) {

    /* First IDABMem struct. */
    tmp_IDAB_mem = IDAB_mem;
    
    while (tmp_IDAB_mem != NULL) {

      tBn = tmp_IDAB_mem->IDA_mem->ida_tn;

      if ( (sign*(tBn-tinitial) < ZERO) || (sign*(tfinal-tBn) < ZERO) ) {
        IDAProcessError(IDA_mem, IDA_BAD_TB0, "IDAA", "IDASolveB", 
                        MSGAM_BAD_TB0, tmp_IDAB_mem->ida_index);
        return(IDA_BAD_TB0);
      }

      if (sign*(tBn-tBout) <= ZERO) {
        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASolveB", MSGAM_BAD_TBOUT,
                       tmp_IDAB_mem->ida_index);
        return(IDA_ILL_INPUT);
      }

      if ( tmp_IDAB_mem->ida_res_withSensi || 
           tmp_IDAB_mem->ida_rhsQ_withSensi ) interpSensi = TRUE;

      /* Advance in list. */
      tmp_IDAB_mem = tmp_IDAB_mem->ida_next;      
    }

    if ( interpSensi && !storeSensi) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASolveB", MSGAM_BAD_SENSI);
      return(IDA_ILL_INPUT);
    }

    IDAADJ_mem->ia_firstIDABcall = FALSE;
  }

  /* Check for valid itask */
  if ( (itaskB != IDA_NORMAL) && (itaskB != IDA_ONE_STEP) ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASolveB", MSG_BAD_ITASK);
    return(IDA_ILL_INPUT);
  }

  /* Check if tBout is legal */
  if ( (sign*(tBout-tinitial) < ZERO) || (sign*(tfinal-tBout) < ZERO) ) {
    tfuzz = HUNDRED*uround*(SUNRabs(tinitial) + SUNRabs(tfinal));
    if ( (sign*(tBout-tinitial) < ZERO) && (SUNRabs(tBout-tinitial) < tfuzz) ) {
      tBout = tinitial;
    } else {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDASolveB", MSGAM_BAD_TBOUT);
      return(IDA_ILL_INPUT);
    }
  }

  /* Loop through the check points and stop as soon as a backward
   * problem has its tn value behind the current check point's t0_
   * value (in the backward direction) */

  ck_mem = IDAADJ_mem->ck_mem;

  gotCkpnt = FALSE;

  loop {
    tmp_IDAB_mem = IDAB_mem;
    while(tmp_IDAB_mem != NULL) {
      tBn = tmp_IDAB_mem->IDA_mem->ida_tn;

      if ( sign*(tBn-t0_) > ZERO ) {
        gotCkpnt = TRUE;
        break;
      }

      if ( (itaskB == IDA_NORMAL) && (tBn == t0_) && (sign*(tBout-t0_) >= ZERO) ) {
        gotCkpnt = TRUE;
        break;
      }

      tmp_IDAB_mem = tmp_IDAB_mem->ida_next;
    }

    if (gotCkpnt) break;

    if (ck_mem->ck_next == NULL) break;

    ck_mem = ck_mem->ck_next;
  }

  /* Loop while propagating backward problems */
  loop {

    /* Store interpolation data if not available.
       This is the 2nd forward integration pass */
    if (ck_mem != ckpntData) {

      flag = IDAAdataStore(IDA_mem, ck_mem);
      if (flag != IDA_SUCCESS) break;
    }

    /* Starting with the current check point from above, loop over check points
       while propagating backward problems */

    tmp_IDAB_mem = IDAB_mem;
    while (tmp_IDAB_mem != NULL) {

      /* Decide if current backward problem is "active" in this check point */
      isActive = TRUE;

      tBn = tmp_IDAB_mem->IDA_mem->ida_tn;

      if ( (tBn == t0_) && (sign*(tBout-t0_) < ZERO ) ) isActive = FALSE;
      if ( (tBn == t0_) && (itaskB == IDA_ONE_STEP) ) isActive = FALSE;
      if ( sign*(tBn - t0_) < ZERO ) isActive = FALSE;

      if ( isActive ) {
        /* Store the address of current backward problem memory 
         * in IDAADJ_mem to be used in the wrapper functions */
        IDAADJ_mem->ia_bckpbCrt = tmp_IDAB_mem;

        /* Integrate current backward problem */
        IDASetStopTime(tmp_IDAB_mem->IDA_mem, t0_);
        flag = IDASolve(tmp_IDAB_mem->IDA_mem, tBout, &tBret, 
                        tmp_IDAB_mem->ida_yy, tmp_IDAB_mem->ida_yp, 
                        itaskB);

        /* Set the time at which we will report solution and/or quadratures */
        tmp_IDAB_mem->ida_tout = tBret;

        /* If an error occurred, exit while loop */
        if (flag < 0) break;

      } else {

        flag = IDA_SUCCESS;
        tmp_IDAB_mem->ida_tout = tBn;
      }

      /* Move to next backward problem */
      tmp_IDAB_mem = tmp_IDAB_mem->ida_next;
    } /* End of while: iteration through backward problems. */
    
    /* If an error occurred, return now */
    if (flag <0) {
      IDAProcessError(IDA_mem, flag, "IDAA", "IDASolveB",
                      MSGAM_BACK_ERROR, tmp_IDAB_mem->ida_index);
      return(flag);
    }

    /* If in IDA_ONE_STEP mode, return now (flag = IDA_SUCCESS) */
    if (itaskB == IDA_ONE_STEP) break;

    /* If all backward problems have succesfully reached tBout, return now */
    reachedTBout = TRUE;

    tmp_IDAB_mem = IDAB_mem;
    while(tmp_IDAB_mem != NULL) {
      if ( sign*(tmp_IDAB_mem->ida_tout - tBout) > ZERO ) {
        reachedTBout = FALSE;
        break;
      }
      tmp_IDAB_mem = tmp_IDAB_mem->ida_next;
    }

    if ( reachedTBout ) break;

    /* Move check point in linked list to next one */
    ck_mem = ck_mem->ck_next;

  } /* End of loop. */

  return(flag);
}


/*
 * IDAGetB
 *
 * IDAGetB returns the state variables at the same time (also returned 
 * in tret) as that at which IDASolveBreturned the solution.
 */

SUNDIALS_EXPORT int IDAGetB(void* ida_mem, int which, realtype *tret,
                            N_Vector yy, N_Vector yp)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAGetB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAGetB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAGetB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  N_VScale(ONE, IDAB_mem->ida_yy, yy);
  N_VScale(ONE, IDAB_mem->ida_yp, yp);
  *tret = IDAB_mem->ida_tout;

  return(IDA_SUCCESS);
}



/*
 * IDAGetQuadB
 *
 * IDAGetQuadB returns the quadrature variables at the same 
 * time (also returned in tret) as that at which IDASolveB 
 * returned the solution.
 */

int IDAGetQuadB(void *ida_mem, int which, realtype *tret, N_Vector qB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void *ida_memB;
  int flag;
  long int nstB;
  
  /* Is ida_mem valid? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAGetQuadB", MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, "IDAA", "IDAGetQuadB",  MSGAM_NO_ADJ);
    return(IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDAA", "IDAGetQuadB", MSGAM_BAD_WHICH);
    return(IDA_ILL_INPUT);
  }
  
  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void *) IDAB_mem->IDA_mem;

  /* If the integration for this backward problem has not started yet,
   * simply return the current value of qB (i.e. the final conditions) */

  flag = IDAGetNumSteps(ida_memB, &nstB);
  if (IDA_SUCCESS != flag) return(flag);

  if (nstB == 0) {
    N_VScale(ONE, IDAB_mem->IDA_mem->ida_phiQ[0], qB);
    *tret = IDAB_mem->ida_tout;
  } else {
    flag = IDAGetQuad(ida_memB, tret, qB);
  }
  return(flag);
}

/*=================================================================*/
/*                Private Functions Implementation                 */
/*=================================================================*/

/*
 * IDAAckpntInit
 *
 * This routine initializes the check point linked list with 
 * information from the initial time.
*/

static CkpntMem IDAAckpntInit(IDAMem IDA_mem)
{
  CkpntMem ck_mem;

  /* Allocate space for ckdata */
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  if (NULL==ck_mem) return(NULL);

  t0_    = tn;
  nst_   = 0;
  kk_    = 1;
  hh_    = ZERO;

  /* Test if we need to carry quadratures */
  quadr_ = quadr && errconQ;

  /* Test if we need to carry sensitivities */
  sensi_ = sensi;
  if(sensi_) Ns_    = Ns;

  /* Test if we need to carry quadrature sensitivities */
  quadr_sensi_ = quadr_sensi && errconQS;

  /* Alloc 3: current order, i.e. 1,  +   2. */
  phi_alloc_ = 3;
  
  if (!IDAAckpntAllocVectors(IDA_mem, ck_mem)) {
    free(ck_mem); ck_mem = NULL;
    return(NULL);
  }
  /* Save phi* vectors from IDA_mem to ck_mem. */
  IDAAckpntCopyVectors(IDA_mem, ck_mem);

  /* Next in list */
  next_  = NULL;

  return(ck_mem);
}

/*
 * IDAAckpntNew
 *
 * This routine allocates space for a new check point and sets 
 * its data from current values in IDA_mem.
*/

static CkpntMem IDAAckpntNew(IDAMem IDA_mem)
{
  CkpntMem ck_mem;
  int j;

  /* Allocate space for ckdata */
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  if (ck_mem == NULL) return(NULL);

  nst_       = nst;
  tretlast_  = tretlast;
  kk_        = kk;
  kused_     = kused;
  knew_      = knew;
  phase_     = phase;
  ns_        = ns;
  hh_        = hh;
  hused_     = hused;
  rr_        = rr;
  cj_        = cj;
  cjlast_    = cjlast;
  cjold_     = cjold;
  cjratio_   = cjratio;
  ss_        = ss;
  ssS_       = ssS;
  t0_        = tn;

  for (j=0; j<MXORDP1; j++) {
    psi_[j]   = psi[j];
    alpha_[j] = alpha[j];
    beta_[j]  = beta[j];
    sigma_[j] = sigma[j];
    gamma_[j] = gamma[j];
  }

  /* Test if we need to carry quadratures */
  quadr_ = quadr && errconQ;

  /* Test if we need to carry sensitivities */
  sensi_ = sensi;
  if(sensi_) Ns_    = Ns;

  /* Test if we need to carry quadrature sensitivities */
  quadr_sensi_ = quadr_sensi && errconQS;

  phi_alloc_ =  kk+2 < MXORDP1 ? kk+2 : MXORDP1;

  if (!IDAAckpntAllocVectors(IDA_mem, ck_mem)) {
    free(ck_mem); ck_mem = NULL;
    return(NULL);
  }

  /* Save phi* vectors from IDA_mem to ck_mem. */
  IDAAckpntCopyVectors(IDA_mem, ck_mem);

  return(ck_mem);
}

/* IDAAckpntDelete 
 *
 * This routine deletes the first check point in list.
*/

static void IDAAckpntDelete(CkpntMem *ck_memPtr)
{
  CkpntMem tmp;
  int j;

  if (*ck_memPtr != NULL) {
    /* store head of list */
    tmp = *ck_memPtr;
    /* move head of list */
    *ck_memPtr = (*ck_memPtr)->ck_next;

    /* free N_Vectors in tmp */
    for (j=0; j<tmp->ck_phi_alloc; j++) 
      N_VDestroy(tmp->ck_phi[j]);

    /* free N_Vectors for quadratures in tmp */
    if (tmp->ck_quadr) {
      for (j=0; j<tmp->ck_phi_alloc; j++) 
        N_VDestroy(tmp->ck_phiQ[j]);
    }

    /* Free sensitivity related data. */
    if (tmp->ck_sensi) {
      for (j=0; j<tmp->ck_phi_alloc; j++) 
        N_VDestroyVectorArray(tmp->ck_phiS[j], tmp->ck_Ns);
    }
    
    if (tmp->ck_quadr_sensi) {
      for (j=0; j<tmp->ck_phi_alloc; j++) 
        N_VDestroyVectorArray(tmp->ck_phiQS[j], tmp->ck_Ns);
    }

    free(tmp); tmp=NULL;
  }
}

/* 
 * IDAAckpntAllocVectors
 *
 * Allocate checkpoint's phi, phiQ, phiS, phiQS vectors needed to save 
 * current state of IDAMem.
 *
 */
static booleantype IDAAckpntAllocVectors(IDAMem IDA_mem, CkpntMem ck_mem)
{
  int j, jj;

  for (j=0; j<phi_alloc_; j++) {
    phi_[j] = N_VClone(tempv);
    if(phi_[j] == NULL) {    
      for(jj=0; jj<j; j++) N_VDestroy(phi_[jj]);
      return(FALSE);
    }
  }

  /* Do we need to carry quadratures? */
  if(quadr_) {
    for (j=0; j<phi_alloc_; j++) {
      phiQ_[j] = N_VClone(tempvQ);
      if(phiQ_[j] == NULL)  {        
        for (jj=0; jj<j; j++) N_VDestroy(phiQ_[jj]);

        for(jj=0; jj<phi_alloc_; j++) N_VDestroy(phi_[jj]);

        return(FALSE);
      }
    }
  }

  /* Do we need to carry sensitivities? */
  if(sensi_) {

    for (j=0; j<phi_alloc_; j++) {
      phiS_[j] = N_VCloneVectorArray(Ns, tempv);
      if (phiS_[j] == NULL) {
        for (jj=0; jj<j; jj++) N_VDestroyVectorArray(phiS_[jj], Ns);

        if (quadr_)
          for (jj=0; jj<phi_alloc_; jj++) N_VDestroy(phiQ_[jj]);

        for (jj=0; jj<phi_alloc_; jj++) N_VDestroy(phi_[jj]);

        return(FALSE);
      }
    }
  }

  /* Do we need to carry quadrature sensitivities? */
  if (quadr_sensi_) {

    for (j=0; j<phi_alloc_; j++) {
      phiQS_[j] = N_VCloneVectorArray(Ns, tempvQ);
      if (phiQS_[j] == NULL) {

        for (jj=0; jj<j; jj++) N_VDestroyVectorArray(phiQS_[jj], Ns);

        for (jj=0; jj<phi_alloc_; jj++) N_VDestroyVectorArray(phiS_[jj], Ns);

        if (quadr_) 
          for (jj=0; jj<phi_alloc_; jj++) N_VDestroy(phiQ_[jj]);

        for (jj=0; jj<phi_alloc_; jj++) N_VDestroy(phi_[jj]);

        return(FALSE);
      }
    }
  }
  return(TRUE);
}

/* 
 * IDAAckpntCopyVectors
 *
 * Copy phi* vectors from IDAMem in the corresponding vectors from checkpoint
 *
 */
static void IDAAckpntCopyVectors(IDAMem IDA_mem, CkpntMem ck_mem)
{
  int j, is;

  /* Save phi* arrays from IDA_mem */

  for (j=0; j<phi_alloc_; j++) N_VScale(ONE, phi[j], phi_[j]);

  if (quadr_) {
    for (j=0; j<phi_alloc_; j++) N_VScale(ONE, phiQ[j], phiQ_[j]);
  }

  if (sensi_) {
    for (is=0; is<Ns; is++)
      for (j=0; j<phi_alloc_; j++) 
        N_VScale(ONE, phiS[j][is], phiS_[j][is]);
  }

  if(quadr_sensi_) {
    for (is=0; is<Ns; is++)
      for (j=0; j<phi_alloc_; j++) 
        N_VScale(ONE, phiQS[j][is], phiQS_[j][is]);

  }
}

/*
 * IDAAdataMalloc
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points. 
 * This data is then used to interpolate the forward solution 
 * at any other time.
*/

static booleantype IDAAdataMalloc(IDAMem IDA_mem)
{
  IDAadjMem IDAADJ_mem;
  DtpntMem *dt_mem;
  long int i, j;

  IDAADJ_mem = IDA_mem->ida_adj_mem;
  IDAADJ_mem->dt_mem = NULL;

  dt_mem = (DtpntMem *)malloc((nsteps+1)*sizeof(struct DtpntMemRec *));
  if (dt_mem==NULL) return(FALSE);

  for (i=0; i<=nsteps; i++) {
    
    dt_mem[i] = (DtpntMem)malloc(sizeof(struct DtpntMemRec));
    
    /* On failure, free any allocated memory and return NULL. */
    if (dt_mem[i] == NULL) {

      for(j=0; j<i; j++) 
        free(dt_mem[j]);

      free(dt_mem);
      return(FALSE);
    }
    dt_mem[i]->content = NULL;
  }
  /* Attach the allocated dt_mem to IDAADJ_mem. */
  IDAADJ_mem->dt_mem = dt_mem;
  return(TRUE);
}

/*
 * IDAAdataFree
 *
 * This routine frees the memory allocated for data storage.
 */

static void IDAAdataFree(IDAMem IDA_mem)
{
  IDAadjMem IDAADJ_mem;
  long int i;

  IDAADJ_mem = IDA_mem->ida_adj_mem;

  if (IDAADJ_mem == NULL) return;

  /* Destroy data points by calling the interpolation's 'free' routine. */
  IDAADJ_mem->ia_free(IDA_mem);

  for (i=0; i<=nsteps; i++) {
     free(IDAADJ_mem->dt_mem[i]);
     IDAADJ_mem->dt_mem[i] = NULL;
  }

  free(IDAADJ_mem->dt_mem);
  IDAADJ_mem->dt_mem = NULL;
}


/*
 * IDAAdataStore 
 *
 * This routine integrates the forward model starting at the check
 * point ck_mem and stores y and yprime at all intermediate 
 * steps. 
 *
 * Return values: 
 *   - the flag that IDASolve may return on error
 *   - IDA_REIFWD_FAIL if no check point is available for this hot start
 *   - IDA_SUCCESS
 */

static int IDAAdataStore(IDAMem IDA_mem, CkpntMem ck_mem)
{
  IDAadjMem IDAADJ_mem;
  DtpntMem *dt_mem;
  realtype t;
  long int i;
  int flag, sign;

  IDAADJ_mem = IDA_mem->ida_adj_mem;
  dt_mem = IDAADJ_mem->dt_mem;

  /* Initialize IDA_mem with data from ck_mem. */
  flag = IDAAckpntGet(IDA_mem, ck_mem);
  if (flag != IDA_SUCCESS)
    return(IDA_REIFWD_FAIL);

  /* Set first structure in dt_mem[0] */
  dt_mem[0]->t = t0_;
  IDAADJ_mem->ia_storePnt(IDA_mem, dt_mem[0]);

  /* Decide whether TSTOP must be activated */
  if (IDAADJ_mem->ia_tstopIDAFcall) {
    IDASetStopTime(IDA_mem, IDAADJ_mem->ia_tstopIDAF);
  }

  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* Run IDASolve in IDA_ONE_STEP mode to set following structures in dt_mem[i]. */
  i = 1;
  do {

    flag = IDASolve(IDA_mem, t1_, &t, yyTmp, ypTmp, IDA_ONE_STEP);
    if (flag < 0) return(IDA_FWD_FAIL);

    dt_mem[i]->t = t;
    IDAADJ_mem->ia_storePnt(IDA_mem, dt_mem[i]);

    i++;
  } while ( sign*(t1_ - t) > ZERO );

  /* New data is now available. */
  ckpntData = ck_mem;
  newData = TRUE;
  np  = i;

  return(IDA_SUCCESS);
}

/*
 * CVAckpntGet
 *
 * This routine prepares IDAS for a hot restart from
 * the check point ck_mem
 */

static int IDAAckpntGet(IDAMem IDA_mem, CkpntMem ck_mem) 
{
  int flag, j, is;

  if (next_ == NULL) {

    /* In this case, we just call the reinitialization routine,
     * but make sure we use the same initial stepsize as on 
     * the first run. */

    IDASetInitStep(IDA_mem, h0u);

    flag = IDAReInit(IDA_mem, t0_, phi_[0], phi_[1]);
    if (flag != IDA_SUCCESS) return(flag);

    if (quadr_) {
      flag = IDAQuadReInit(IDA_mem, phiQ_[0]);
      if (flag != IDA_SUCCESS) return(flag);
    }

    if (sensi_) {
      flag = IDASensReInit(IDA_mem, IDA_mem->ida_ism, phiS_[0], phiS_[1]);
      if (flag != IDA_SUCCESS) return(flag);
    }

    if (quadr_sensi_) {
      flag = IDAQuadSensReInit(IDA_mem, phiQS_[0]);
      if (flag != IDA_SUCCESS) return(flag);
    }

  } else {

    /* Copy parameters from check point data structure */
    nst       = nst_;
    tretlast  = tretlast_;
    kk        = kk_;
    kused     = kused_;
    knew      = knew_;
    phase     = phase_;
    ns        = ns_;
    hh        = hh_;
    hused     = hused_;
    rr        = rr_;
    cj        = cj_;
    cjlast    = cjlast_;
    cjold     = cjold_;
    cjratio   = cjratio_;
    tn        = t0_;
    ss        = ss_;
    ssS       = ssS_;

    
    /* Copy the arrays from check point data structure */
    for (j=0; j<phi_alloc_; j++) N_VScale(ONE, phi_[j], phi[j]);

    if(quadr_) {
      for (j=0; j<phi_alloc_; j++) N_VScale(ONE, phiQ_[j], phiQ[j]);
    }

    if (sensi_) {
      for (is=0; is<Ns; is++) {
        for (j=0; j<phi_alloc_; j++) N_VScale(ONE, phiS_[j][is], phiS[j][is]);
      }
    }

    if (quadr_sensi_) {
      for (is=0; is<Ns; is++) {
        for (j=0; j<phi_alloc_; j++) N_VScale(ONE, phiQS_[j][is], phiQS[j][is]);
      }
    }

    for (j=0; j<MXORDP1; j++) {
      psi[j]   = psi_[j];
      alpha[j] = alpha_[j];
      beta[j]  = beta_[j];
      sigma[j] = sigma_[j];
      gamma[j] = gamma_[j];
    }

    /* Force a call to setup */
    forceSetup = TRUE;
  }

  return(IDA_SUCCESS);
}


/* 
 * -----------------------------------------------------------------
 * Functions specific to cubic Hermite interpolation
 * -----------------------------------------------------------------
 */

/*
 * IDAAhermiteMalloc
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points. 
 * This data is then used to interpolate the forward solution 
 * at any other time.
 */

static booleantype IDAAhermiteMalloc(IDAMem IDA_mem)
{
  IDAadjMem IDAADJ_mem;
  DtpntMem *dt_mem;
  HermiteDataMem content;
  long int i, ii=0;
  booleantype allocOK;

  allocOK = TRUE;

  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Allocate space for the vectors yyTmp and ypTmp. */
  yyTmp = N_VClone(tempv);
  if (yyTmp == NULL) {
    return(FALSE);
  }
  ypTmp = N_VClone(tempv);
  if (ypTmp == NULL) {
    return(FALSE);
  }

  /* Allocate space for sensitivities temporary vectors. */
  if (storeSensi) {
    
    yySTmp = N_VCloneVectorArray(Ns, tempv);
    if (yySTmp == NULL) {
      N_VDestroy(yyTmp);
      N_VDestroy(ypTmp);
      return(FALSE);
    }

    ypSTmp = N_VCloneVectorArray(Ns, tempv);
    if (ypSTmp == NULL) {
      N_VDestroy(yyTmp);
      N_VDestroy(ypTmp);
      N_VDestroyVectorArray(yySTmp, Ns);
      return(FALSE);

    }
  }

  /* Allocate space for the content field of the dt structures */

  dt_mem = IDAADJ_mem->dt_mem;

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

    if (storeSensi) {
      
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

    N_VDestroy(yyTmp);
    N_VDestroy(ypTmp);  

    if (storeSensi) {     
      N_VDestroyVectorArray(yySTmp, Ns);
      N_VDestroyVectorArray(ypSTmp, Ns);
    }

    for (i=0; i<ii; i++) {
      content = (HermiteDataMem) (dt_mem[i]->content);
      N_VDestroy(content->y);
      N_VDestroy(content->yd);

      if (storeSensi) {
        N_VDestroyVectorArray(content->yS, Ns);
        N_VDestroyVectorArray(content->ySd, Ns);        
      }

      free(dt_mem[i]->content); dt_mem[i]->content = NULL;
    }

  }

  return(allocOK);
}

/*
 * IDAAhermiteFree
 *
 * This routine frees the memory allocated for data storage.
 */

static void IDAAhermiteFree(IDAMem IDA_mem)
{  
  IDAadjMem IDAADJ_mem;
  DtpntMem *dt_mem;
  HermiteDataMem content;
  long int i;

  IDAADJ_mem = IDA_mem->ida_adj_mem;

  N_VDestroy(yyTmp);
  N_VDestroy(ypTmp);

  if (storeSensi) {    
    N_VDestroyVectorArray(yySTmp, Ns);
    N_VDestroyVectorArray(ypSTmp, Ns);
  }

  dt_mem = IDAADJ_mem->dt_mem;

  for (i=0; i<=nsteps; i++) {

    content = (HermiteDataMem) (dt_mem[i]->content);
    /* content might be NULL, if IDAAdjInit was called but IDASolveF was not. */
    if(content) {

      N_VDestroy(content->y);
      N_VDestroy(content->yd);

      if (storeSensi) {
        N_VDestroyVectorArray(content->yS, Ns);
        N_VDestroyVectorArray(content->ySd, Ns);      
      }
      free(dt_mem[i]->content); 
      dt_mem[i]->content = NULL;
    }
  }
}

/*
 * IDAAhermiteStorePnt
 *
 * This routine stores a new point (y,yd) in the structure d for use
 * in the cubic Hermite interpolation.
 * Note that the time is already stored.
 */

static int IDAAhermiteStorePnt(IDAMem IDA_mem, DtpntMem d)
{
  IDAadjMem IDAADJ_mem;
  HermiteDataMem content;
  int is;

  IDAADJ_mem = IDA_mem->ida_adj_mem;

  content = (HermiteDataMem) d->content;

  /* Load solution(s) */
  N_VScale(ONE, phi[0], content->y);
  
  if (storeSensi) {
    for (is=0; is<Ns; is++) 
      N_VScale(ONE, phiS[0][is], content->yS[is]);
  }

  /* Load derivative(s). */
  IDAAGettnSolutionYp(IDA_mem, content->yd);

  if (storeSensi) {
    IDAAGettnSolutionYpS(IDA_mem, content->ySd);
  }

  return(0);
}


/*
 * IDAAhermiteGetY
 *
 * This routine uses cubic piece-wise Hermite interpolation for 
 * the forward solution vector. 
 * It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB) but
 * can be directly called by the user through IDAGetAdjY
 */

static int IDAAhermiteGetY(IDAMem IDA_mem, realtype t,
                           N_Vector yy, N_Vector yp,
                          N_Vector *yyS, N_Vector *ypS)
{
  IDAadjMem IDAADJ_mem;
  DtpntMem *dt_mem;
  HermiteDataMem content0, content1;

  realtype t0, t1, delta;
  realtype factor1, factor2, factor3;

  N_Vector y0, yd0, y1, yd1;
  N_Vector *yS0=NULL, *ySd0=NULL, *yS1, *ySd1;

  int flag, is, NS;
  long int indx;
  booleantype newpoint;

 
  IDAADJ_mem = IDA_mem->ida_adj_mem;
  dt_mem = IDAADJ_mem->dt_mem;
 
  /* Local value of Ns */
  NS = interpSensi ? Ns : 0;

  /* Get the index in dt_mem */
  flag = IDAAfindIndex(IDA_mem, t, &indx, &newpoint);
  if (flag != IDA_SUCCESS) return(flag);

  /* If we are beyond the left limit but close enough,
     then return y at the left limit. */

  if (indx == 0) {
    content0 = (HermiteDataMem) (dt_mem[0]->content);
    N_VScale(ONE, content0->y,  yy);
    N_VScale(ONE, content0->yd, yp);

    for (is=0; is<NS; is++) {
      N_VScale(ONE, content0->yS[is], yyS[is]);
      N_VScale(ONE, content0->ySd[is],ypS[is]);
    }
    return(IDA_SUCCESS);
  }

  /* Extract stuff from the appropriate data points */
  t0 = dt_mem[indx-1]->t;
  t1 = dt_mem[indx]->t;
  delta = t1 - t0;

  content0 = (HermiteDataMem) (dt_mem[indx-1]->content);
  y0  = content0->y;
  yd0 = content0->yd;
  if (interpSensi) {
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

  /* For y. */
  factor1 = t - t0;

  factor2 = factor1/delta;
  factor2 = factor2*factor2;

  factor3 = factor2*(t-t1)/delta;

  N_VLinearSum(ONE, y0, factor1, yd0, yy);
  N_VLinearSum(ONE, yy, factor2, Y[0], yy);
  N_VLinearSum(ONE, yy, factor3, Y[1], yy);

  /* Sensi Interpolation. */
  for (is=0; is<NS; is++) {
    N_VLinearSum(ONE, yS0[is], factor1, ySd0[is], yyS[is]);
    N_VLinearSum(ONE, yyS[is], factor2, YS[0][is], yyS[is]);
    N_VLinearSum(ONE, yyS[is], factor3, YS[1][is], yyS[is]);
  }

  /*For y'. */
  factor1 = factor1/delta/delta; /* factor1 = 2(t-t0)/(t1-t0)^2 */
  factor2 = factor1*((3*t-2*t1-t0)/delta); /* factor2 = (t-t0)(3*t-2*t1-t0)/(t1-t0)^3 */
  factor1 *= 2;

  N_VLinearSum(ONE, yd0, factor1, Y[0], yp);
  N_VLinearSum(ONE, yp,  factor2, Y[1], yp);
                                            
  /* Sensi interpolation for 1st derivative. */
  for (is=0; is<NS; is++) {
    N_VLinearSum(ONE, ySd0[is], factor1, YS[0][is], ypS[is]);
    N_VLinearSum(ONE, ypS[is],  factor2, YS[1][is], ypS[is]);    
  }

  return(IDA_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Functions specific to Polynomial interpolation
 * -----------------------------------------------------------------
 */

/*
 * IDAApolynomialMalloc 
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points. 
 * This data is then used to interpolate the forward solution 
 * at any other time.
 *
 * Information about the first derivative is stored only for the first
 * data point.
 */

static booleantype IDAApolynomialMalloc(IDAMem IDA_mem)
{
  IDAadjMem IDAADJ_mem;
  DtpntMem *dt_mem;
  PolynomialDataMem content;
  long int i, ii=0;
  booleantype allocOK;

  allocOK = TRUE;

  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Allocate space for the vectors yyTmp and ypTmp */
  yyTmp = N_VClone(tempv);
  if (yyTmp == NULL) {
    return(FALSE);
  }
  ypTmp = N_VClone(tempv);
  if (ypTmp == NULL) {
    return(FALSE);
  }

  if (storeSensi) {
    
    yySTmp = N_VCloneVectorArray(Ns, tempv);
    if (yySTmp == NULL) {
      N_VDestroy(yyTmp);
      N_VDestroy(ypTmp);
      return(FALSE);
    }

    ypSTmp = N_VCloneVectorArray(Ns, tempv);
    if (ypSTmp == NULL) {
      N_VDestroy(yyTmp);
      N_VDestroy(ypTmp);
      N_VDestroyVectorArray(yySTmp, Ns);
      return(FALSE);

    }
  }

  /* Allocate space for the content field of the dt structures */
  dt_mem = IDAADJ_mem->dt_mem;

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

    /* Allocate space for yp also. Needed for the most left point interpolation. */
    if (i == 0) {
      content->yd = N_VClone(tempv);
      
      /* Memory allocation failure ? */
      if (content->yd == NULL) {
        N_VDestroy(content->y);
        free(content); content = NULL;
        ii = i;
        allocOK = FALSE;
      }
    } else {
      /* Not the first data point. */
      content->yd = NULL;
    }

    if (storeSensi) {
      
      content->yS = N_VCloneVectorArray(Ns, tempv);
      if (content->yS == NULL) {
        N_VDestroy(content->y);
        if (content->yd) N_VDestroy(content->yd);
        free(content); content = NULL;
        ii = i;
        allocOK = FALSE;
        break;
      }
      
      if (i==0) {
        content->ySd = N_VCloneVectorArray(Ns, tempv);
        if (content->ySd == NULL) {
          N_VDestroy(content->y);
          if (content->yd) N_VDestroy(content->yd);
          N_VDestroyVectorArray(content->yS, Ns);
          free(content); content = NULL;
          ii = i;
          allocOK = FALSE;
        }
      } else {
        content->ySd = NULL;
      }
    }

    dt_mem[i]->content = content;
  } 

  /* If an error occurred, deallocate and return */
  if (!allocOK) {

    N_VDestroy(yyTmp);
    N_VDestroy(ypTmp);
    if (storeSensi) {

        N_VDestroyVectorArray(yySTmp, Ns);
        N_VDestroyVectorArray(ypSTmp, Ns);      
    }

    for (i=0; i<ii; i++) {
      content = (PolynomialDataMem) (dt_mem[i]->content);
      N_VDestroy(content->y);

      if (content->yd) N_VDestroy(content->yd);

      if (storeSensi) {
        
          N_VDestroyVectorArray(content->yS, Ns);
        
          if (content->ySd)
            N_VDestroyVectorArray(content->ySd, Ns);
      }
      free(dt_mem[i]->content); dt_mem[i]->content = NULL;
    }

  }
  return(allocOK);
}

/*
 * IDAApolynomialFree
 *
 * This routine frees the memory allocated for data storage.
 */

static void IDAApolynomialFree(IDAMem IDA_mem)
{
  IDAadjMem IDAADJ_mem;
  DtpntMem *dt_mem;
  PolynomialDataMem content;
  long int i;

  IDAADJ_mem = IDA_mem->ida_adj_mem;

  N_VDestroy(yyTmp);
  N_VDestroy(ypTmp);

  if (storeSensi) {
    N_VDestroyVectorArray(yySTmp, Ns);
    N_VDestroyVectorArray(ypSTmp, Ns);
  }

  dt_mem = IDAADJ_mem->dt_mem;

  for (i=0; i<=nsteps; i++) {

    content = (PolynomialDataMem) (dt_mem[i]->content);

    /* content might be NULL, if IDAAdjInit was called but IDASolveF was not. */
    if(content) {
      N_VDestroy(content->y);

      if (content->yd) N_VDestroy(content->yd);

      if (storeSensi) {
        
        N_VDestroyVectorArray(content->yS, Ns);
        
        if (content->ySd)
          N_VDestroyVectorArray(content->ySd, Ns);
      }
      free(dt_mem[i]->content); dt_mem[i]->content = NULL;
    }
  }
}

/*
 * IDAApolynomialStorePnt
 *
 * This routine stores a new point y in the structure d for use
 * in the Polynomial interpolation.
 *
 * Note that the time is already stored. Information about the 
 * first derivative is available only for the first data point, 
 * in which case content->yp is non-null.
 */

static int IDAApolynomialStorePnt(IDAMem IDA_mem, DtpntMem d)
{
  IDAadjMem IDAADJ_mem;
  PolynomialDataMem content;
  int is;

  IDAADJ_mem = IDA_mem->ida_adj_mem;
  content = (PolynomialDataMem) d->content;

  N_VScale(ONE, phi[0], content->y);

  /* copy also the derivative for the first data point (in this case
     content->yp is non-null). */
  if (content->yd)
    IDAAGettnSolutionYp(IDA_mem, content->yd);

  if (storeSensi) {
    
    for (is=0; is<Ns; is++) 
      N_VScale(ONE, phiS[0][is], content->yS[is]);
    
    /* store the derivative if it is the first data point. */
    if(content->ySd)
      IDAAGettnSolutionYpS(IDA_mem, content->ySd);
  }

  content->order = kused;

  return(0);
}

/*
 * IDAApolynomialGetY
 *
 * This routine uses polynomial interpolation for the forward solution vector. 
 * It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB)) but
 * can be directly called by the user through CVodeGetAdjY.
 */

static int IDAApolynomialGetY(IDAMem IDA_mem, realtype t,
                              N_Vector yy, N_Vector yp,
                              N_Vector *yyS, N_Vector *ypS)
{
  IDAadjMem IDAADJ_mem;
  DtpntMem *dt_mem;
  PolynomialDataMem content;

  int flag, dir, order, i, j, is, NS;
  long int indx, base;
  booleantype newpoint;
  realtype delt, factor, Psi, Psiprime;

  IDAADJ_mem = IDA_mem->ida_adj_mem;
  dt_mem = IDAADJ_mem->dt_mem;
 
  /* Local value of Ns */
 
  NS = interpSensi ? Ns : 0;

  /* Get the index in dt_mem */

  flag = IDAAfindIndex(IDA_mem, t, &indx, &newpoint);
  if (flag != IDA_SUCCESS) return(flag);

  /* If we are beyond the left limit but close enough,
     then return y at the left limit. */

  if (indx == 0) {
    content = (PolynomialDataMem) (dt_mem[0]->content);
    N_VScale(ONE, content->y,  yy);
    N_VScale(ONE, content->yd, yp);

    
    for (is=0; is<NS; is++) {
      N_VScale(ONE, content->yS[is], yyS[is]);
      N_VScale(ONE, content->ySd[is], ypS[is]);
    }
    
    return(IDA_SUCCESS);
  }

  /* Scaling factor */
  delt = SUNRabs(dt_mem[indx]->t - dt_mem[indx-1]->t);

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
        
        for (is=0; is<NS; is++) 
          N_VScale(ONE, content->yS[is], YS[j][is]);
       
      }
    } else {
      for(j=0;j<=order;j++) {
        T[j] = dt_mem[base-1+j]->t;
        content = (PolynomialDataMem) (dt_mem[base-1+j]->content);
        N_VScale(ONE, content->y, Y[j]);
        
        for (is=0; is<NS; is++) 
          N_VScale(ONE, content->yS[is], YS[j][is]);
        
      }
    }

    /* Compute higher-order DD */
    for(i=1;i<=order;i++) {
      for(j=order;j>=i;j--) {
        factor = delt/(T[j]-T[j-i]);
        N_VLinearSum(factor, Y[j], -factor, Y[j-1], Y[j]);
        
        for (is=0; is<NS; is++) 
          N_VLinearSum(factor, YS[j][is], -factor, YS[j-1][is], YS[j][is]);
        
      }
    }
  }

  /* Perform the actual interpolation for yy using nested multiplications */
  N_VScale(ONE, Y[order], yy);
  
  for (is=0; is<NS; is++) 
    N_VScale(ONE, YS[order][is], yyS[is]);
  
  for (i=order-1; i>=0; i--) {
    factor = (t-T[i])/delt;
    N_VLinearSum(factor, yy, ONE, Y[i], yy);
    
    for (is=0; is<NS; is++) 
      N_VLinearSum(factor, yyS[is], ONE, YS[i][is], yyS[is]);
    
  }
  
  /* Perform the actual interpolation for yp.

     Writing p(t) = y0 + (t-t0)*f[t0,t1] + ... + (t-t0)(t-t1)...(t-tn)*f[t0,t1,...tn],
     denote psi_k(t) = (t-t0)(t-t1)...(t-tk).

     The formula used for p'(t) is: 
       - p'(t) = f[t0,t1] + psi_1'(t)*f[t0,t1,t2] + ... + psi_n'(t)*f[t0,t1,...,tn]
     
     We reccursively compute psi_k'(t) from:
       - psi_k'(t) = (t-tk)*psi_{k-1}'(t) + psi_{k-1}

     psi_k is rescaled with 1/delt each time is computed, because the Newton DDs from Y were
     scaled with delt.
  */

  Psi = ONE; Psiprime = ZERO; 
  N_VConst(ZERO, yp);

  for (is=0; is<NS; is++)
    N_VConst(ZERO, ypS[is]);

  for(i=1; i<=order; i++) {
    factor = (t-T[i-1])/delt;

    Psiprime = Psi/delt +  factor * Psiprime;
    Psi = Psi * factor;

    N_VLinearSum(ONE, yp, Psiprime, Y[i], yp);
    
    for (is=0; is<NS; is++)
      N_VLinearSum(ONE, ypS[is], Psiprime, YS[i][is], ypS[is]);
  }

  return(IDA_SUCCESS);
}

/* 
 * IDAAGetSolutionYp
 *
 * Evaluates the first derivative of the solution at the last time returned by
 * IDASolve (tretlast).
 * 
 * The function implements the same algorithm as in IDAGetSolution but in the 
 * particular case when  t=tn (i.e. delta=0).
 *
 * This function was implemented to avoid calls to IDAGetSolution which computes 
 * y by doing a loop that is not necessary for this particular situation.
 */

static int IDAAGettnSolutionYp(IDAMem IDA_mem, N_Vector yp)
{
  int j, kord;
  realtype C, D, gam;

  if (nst==0) {

    /* If no integration was done, return the yp supplied by user.*/
      N_VScale(ONE, phi[1], yp);

    return(0);
  }

  /* Compute yp as in IDAGetSolution for this particular case when t=tn. */
  N_VConst(ZERO, yp);
  
  kord = kused;
  if(kused==0) kord=1;
  
  C = ONE; D = ZERO;
  gam = ZERO;
  for (j=1; j <= kord; j++) {
    D = D*gam + C/psi[j-1];
    C = C*gam;
    gam = psi[j-1]/psi[j];
    N_VLinearSum(ONE, yp, D, phi[j], yp);
  }   

  return(0);
}


/* 
 * IDAAGettnSolutionYpS
 *
 * Same as IDAAGettnSolutionYp, but for first derivative of the sensitivities.
 *
 */

static int IDAAGettnSolutionYpS(IDAMem IDA_mem, N_Vector *ypS)
{
  int j, kord, is;
  realtype C, D, gam;

  if (nst==0) {

    /* If no integration was done, return the ypS supplied by user.*/
    for (is=0; is<Ns; is++) 
      N_VScale(ONE, phiS[1][is], ypS[is]);

    return(0);
  }

  for (is=0; is<Ns; is++) 
    N_VConst(ZERO, ypS[is]);
  
  kord = kused;
  if(kused==0) kord=1;
  
  C = ONE; D = ZERO;
  gam = ZERO;
  for (j=1; j <= kord; j++) {
    D = D*gam + C/psi[j-1];
    C = C*gam;
    gam = psi[j-1]/psi[j];
  
    for (is=0; is<Ns; is++)
      N_VLinearSum(ONE, ypS[is], D, phiS[j][is], ypS[is]);
  }   

  return(0);
}



/*
 * IDAAfindIndex
 *
 * Finds the index in the array of data point strctures such that
 *     dt_mem[indx-1].t <= t < dt_mem[indx].t
 * If indx is changed from the previous invocation, then newpoint = TRUE
 *
 * If t is beyond the leftmost limit, but close enough, indx=0.
 *
 * Returns IDA_SUCCESS if successful and IDA_GETY_BADT if unable to
 * find indx (t is too far beyond limits).
 */

static int IDAAfindIndex(IDAMem ida_mem, realtype t, 
                        long int *indx, booleantype *newpoint)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  static long int ilast;
  DtpntMem *dt_mem;
  int sign;
  booleantype to_left, to_right;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;
  dt_mem = IDAADJ_mem->dt_mem;

  *newpoint = FALSE;

  /* Find the direction of integration */
  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* If this is the first time we use new data */
  if (newData) {
    ilast     = np-1;
    *newpoint = TRUE;
    newData   = FALSE;
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
        return(IDA_GETY_BADT);
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
  return(IDA_SUCCESS);
}


/*
 * IDAGetAdjY
 *
 * This routine returns the interpolated forward solution at time t.
 * The user must allocate space for y.
 */

int IDAGetAdjY(void *ida_mem, realtype t, N_Vector yy, N_Vector yp)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  int flag;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDAA", "IDAGetAdjY", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;                              
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  flag = IDAADJ_mem->ia_getY(IDA_mem, t, yy, yp, NULL, NULL);

  return(flag);
}

/*=================================================================*/
/*             Wrappers for adjoint system                         */
/*=================================================================*/

/*
 * IDAAres
 *
 * This routine interfaces to the RhsFnB routine provided by
 * the user.
*/

static int IDAAres(realtype tt, 
                   N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                   void *ida_mem)
{
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  IDAMem IDA_mem;
  int flag, retval;

  IDA_mem = (IDAMem) ida_mem;

  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Get the current backward problem. */
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;

  /* Get forward solution from interpolation. */
  if( noInterp == FALSE) {
    if (interpSensi)
      flag = IDAADJ_mem->ia_getY(ida_mem, tt, yyTmp, ypTmp, yySTmp, ypSTmp);
    else
      flag = IDAADJ_mem->ia_getY(ida_mem, tt, yyTmp, ypTmp, NULL, NULL);
  
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDA_mem, -1, "IDAA", "IDAAres", MSGAM_BAD_TINTERP, tt);
      return(-1);
    }
  }

  /* Call the user supplied residual. */
  if(IDAB_mem->ida_res_withSensi) {
    retval = IDAB_mem->ida_resS(tt, yyTmp, ypTmp, 
                                yySTmp, ypSTmp,
                                yyB, ypB, 
                                rrB, IDAB_mem->ida_user_data);
  }else {
    retval = IDAB_mem->ida_res(tt, yyTmp, ypTmp, yyB, ypB, rrB, IDAB_mem->ida_user_data);
  }
  return(retval);
}

/*
 *IDAArhsQ
 *
 * This routine interfaces to the IDAQuadRhsFnB routine provided by
 * the user.
 *
 * It is passed to IDAQuadInit calls for backward problem, so it must
 * be of IDAQuadRhsFn type.
*/

static int IDAArhsQ(realtype tt, 
                    N_Vector yyB, N_Vector ypB,
                    N_Vector resvalQB, void *ida_mem)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  int retval, flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Get current backward problem. */
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;

  retval = IDA_SUCCESS;

  /* Get forward solution from interpolation. */
  if (noInterp == FALSE) {
    if (interpSensi) {
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, yySTmp, ypSTmp);
    } else {
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, NULL, NULL);
    }
    
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDA_mem, -1, "IDAA", "IDAArhsQ", MSGAM_BAD_TINTERP, tt);
      return(-1);
    }  
  }

  /* Call user's adjoint quadrature RHS routine */
  if (IDAB_mem->ida_rhsQ_withSensi) {
    retval = IDAB_mem->ida_rhsQS(tt, yyTmp, ypTmp, yySTmp, ypSTmp, 
                                 yyB, ypB, 
                                 resvalQB, IDAB_mem->ida_user_data);
  } else {
    retval = IDAB_mem->ida_rhsQ(tt,
                                yyTmp, ypTmp,
                                yyB, ypB,
                                resvalQB, IDAB_mem->ida_user_data);
  }
  return(retval);
}


