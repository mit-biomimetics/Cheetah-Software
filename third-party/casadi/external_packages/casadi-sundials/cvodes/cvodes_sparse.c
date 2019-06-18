/*
 * -----------------------------------------------------------------
 * $Revision: 4086 $
 * $Date: 2014-04-30 16:04:45 -0700 (Wed, 30 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
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
 * This is the implementation file for an CVSLS linear solver.
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
#include "cvodes_sparse_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PRIVATE FUNCTION PROTOTYPES
 * =================================================================
 */

static int cvSlsSparseJacBWrapper(realtype t,
				  N_Vector yB, N_Vector fyB, 
				  SlsMat JB, void *cvode_mem,
				  N_Vector tmp1B, N_Vector tmp2B, 
				  N_Vector tmp3B);

static int cvSlsSparseJacBSWrapper(realtype t,
				   N_Vector yB, N_Vector fyB, 
				   SlsMat JB, void *cvode_mem,
				   N_Vector tmp1B, N_Vector tmp2B, 
				   N_Vector tmp3B);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * CVSlsSetSparseJacFn specifies the sparse Jacobian function.
 */
int CVSlsSetSparseJacFn(void *cvode_mem, CVSlsSparseJacFn jac)
{
  CVodeMem cv_mem;
  CVSlsMem cvsls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSLS", "CVSlsSetSparseJacFn", 
		    MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSLS_LMEM_NULL, "CVSLS", 
		    "CVSlsSetSparseJacFn", MSGSP_LMEM_NULL);
    return(CVSLS_LMEM_NULL);
  }
  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;

  cvsls_mem->s_jaceval = jac;

  return(CVSLS_SUCCESS);
}

/*
 * CVSlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int CVSlsGetNumJacEvals(void *cvode_mem, long int *njevals)
{
  CVodeMem cv_mem;
  CVSlsMem cvsls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSLS", "CVSlsGetNumJacEvals", MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSLS_LMEM_NULL, "CVSLS", 
		    "CVSlsGetNumJacEvals", MSGSP_LMEM_NULL);
    return(CVSLS_LMEM_NULL);
  }
  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;

  *njevals = cvsls_mem->s_nje;

  return(CVSLS_SUCCESS);
}

/*
 * CVSlsGetReturnFlagName returns the name associated with a CVSLS
 * return value.
 */
char *CVSlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVSLS_SUCCESS:
    sprintf(name,"CVSLS_SUCCESS");
    break;   
  case CVSLS_MEM_NULL:
    sprintf(name,"CVSLS_MEM_NULL");
    break;
  case CVSLS_LMEM_NULL:
    sprintf(name,"CVSLS_LMEM_NULL");
    break;
  case CVSLS_ILL_INPUT:
    sprintf(name,"CVSLS_ILL_INPUT");
    break;
  case CVSLS_MEM_FAIL:
    sprintf(name,"CVSLS_MEM_FAIL");
    break;
  case CVSLS_JAC_NOSET:
    sprintf(name,"CVSLS_JAC_NOSET");
    break;
  case CVSLS_JACFUNC_UNRECVR:
    sprintf(name,"CVSLS_JACFUNC_UNRECVR");
    break;
  case CVSLS_JACFUNC_RECVR:
    sprintf(name,"CVSLS_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * CVSlsGetLastFlag returns the last flag set in a CVSLS function.
 */
int CVSlsGetLastFlag(void *cvode_mem, long int *flag)
{
  CVodeMem cv_mem;
  CVSlsMem cvsls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSLS", "CVSlsGetLastFlag", 
		    MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSLS_LMEM_NULL, "CVSLS", 
		    "CVSlsGetLastFlag", MSGSP_LMEM_NULL);
    return(CVSLS_LMEM_NULL);
  }
  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;

  *flag = cvsls_mem->s_last_flag;

  return(CVSLS_SUCCESS);
}

/* 
 * =================================================================
 * BACKWARD INTEGRATION SUPPORT
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * EXPORTED FUNCTIONS
 * -----------------------------------------------------------------
 */

int CVSlsSetSparseJacFnB(void *cvode_mem, int which, CVSlsSparseJacFnB jacB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSlsMemB cvslsB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSSLS", "CVSlsSetSparseJacFnB", MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSLS_NO_ADJ, "CVSSLS", "CVSlsSetSparseJacFnB", MSGSP_NO_ADJ);
    return(CVSLS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSLS_ILL_INPUT, "CVSSLS", "CVSlsSetSparseJacFnB", MSGSP_BAD_WHICH);
    return(CVSLS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  if (cvB_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSLS_LMEMB_NULL, "CVSSLS", "CVSlsSetSparseJacFnB", MSGSP_LMEMB_NULL);
    return(CVSLS_LMEMB_NULL);
  }
  cvslsB_mem = (CVSlsMemB) (cvB_mem->cv_lmem);

  cvslsB_mem->s_djacB = jacB;

  if (jacB != NULL) {
    flag = CVSlsSetSparseJacFn(cvodeB_mem, cvSlsSparseJacBWrapper);
  } else {
    flag = CVSlsSetSparseJacFn(cvodeB_mem, NULL);
  }

  return(flag);
}

int CVSlsSetSparseJacFnBS(void *cvode_mem, int which, CVSlsSparseJacFnBS jacBS)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSlsMemB cvslsB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSSLS", "CVSlsSetSparseJacFnBS", MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSLS_NO_ADJ, "CVSSLS", "CVSlsSetSparseJacFnBS", MSGSP_NO_ADJ);
    return(CVSLS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSLS_ILL_INPUT, "CVSSLS", "CVSlsSetSparseJacFnBS", MSGSP_BAD_WHICH);
    return(CVSLS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  if (cvB_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSLS_LMEMB_NULL, "CVSSLS", "CVSlsSetSparseJacFnBS", MSGSP_LMEMB_NULL);
    return(CVSLS_LMEMB_NULL);
  }
  cvslsB_mem = (CVSlsMemB) (cvB_mem->cv_lmem);

  cvslsB_mem->s_djacBS = jacBS;

  if (jacBS != NULL) {
    flag = CVSlsSetSparseJacFn(cvodeB_mem, cvSlsSparseJacBSWrapper);
  } else {
    flag = CVSlsSetSparseJacFn(cvodeB_mem, NULL);
  }

  return(flag);
}

/*
 * -----------------------------------------------------------------
 * PRIVATE INTERFACE FUNCTIONS
 * -----------------------------------------------------------------
 */

/*
 * cvSlsSparseJacBWrapper
 *
 * This routine interfaces to the CVSlsSparseJacFnB routine provided 
 * by the user. cvSlsSparseJacBWrapper is of type CVSlsSparseJacFn.
 * NOTE: data here contains cvode_mem
 */

static int cvSlsSparseJacBWrapper(realtype t,
				  N_Vector yB, N_Vector fyB, 
				  SlsMat JB, void *cvode_mem,
				  N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSlsMemB cvslsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvslsB_mem = (CVSlsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, NULL);
  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVSSLS", "cvSlsSparseJacBWrapper", MSGSP_BAD_TINTERP);
    return(-1);
  }

  /* Call user's adjoint dense djacB routine (of type CVSlsSparseJacFnB) */
  retval = cvslsB_mem->s_djacB(t, ca_mem->ca_ytmp, yB, fyB, JB, 
			       cvB_mem->cv_user_data, 
			       tmp1B, tmp2B, tmp3B);

  return(retval);
}

/*
 * cvSlsSparseJacBSWrapper
 *
 * This routine interfaces to the CVSlsSparseJacFnBS routine provided 
 * by the user. cvSlsSparseJacBSWrapper is of type CVSlsSparseJacFn.
 * NOTE: data here contains cvode_mem
 */

static int cvSlsSparseJacBSWrapper(realtype t,
				   N_Vector yB, N_Vector fyB, 
				   SlsMat JB, void *cvode_mem,
				   N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  CVSlsMemB cvslsB_mem;
  int retval, flag;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  cvslsB_mem = (CVSlsMemB) (cvB_mem->cv_lmem);

  /* Forward solution from interpolation */
  if (ca_mem->ca_IMinterpSensi)
    flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, ca_mem->ca_yStmp);
  else 
    flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, NULL);
  if (flag != CV_SUCCESS) {
    cvProcessError(cv_mem, -1, "CVSSLS", "cvSlsSparseJacBSWrapper", MSGSP_BAD_TINTERP);
    return(-1);
  }

  /* Call user's adjoint dense djacBS routine (of type CVSlsSparseJacFnBS) */
  retval = cvslsB_mem->s_djacBS(t, ca_mem->ca_ytmp, ca_mem->ca_yStmp, 
				yB, fyB, JB, cvB_mem->cv_user_data, 
				tmp1B, tmp2B, tmp3B);

  return(retval);
}
