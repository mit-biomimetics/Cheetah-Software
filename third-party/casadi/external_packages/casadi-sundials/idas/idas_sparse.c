/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
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
 * This is the implementation file for an IDASSLS linear solver.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "idas_impl.h"
#include "idas_sparse_impl.h"
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
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

static int idaSlsSparseJacBWrapper(realtype tt, realtype c_jB,
		            N_Vector yyB, N_Vector ypB, N_Vector rBr, 
		            SlsMat JacMat, void *ida_mem,
			    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

static int idaSlsSparseJacBSWrapper(realtype tt, realtype c_jB,
			      N_Vector yyB, N_Vector ypB, N_Vector rBr, 
			      SlsMat JacMat, void *ida_mem,
			      N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * IDASlsSetSparseJacFn specifies the sparse Jacobian function.
 */
int IDASlsSetSparseJacFn(void *ida_mem, IDASlsSparseJacFn jac)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDASlsSetSparseJacFn", 
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_LMEM_NULL, "IDASSLS", 
		    "IDASlsSetSparseJacFn", MSGSP_LMEM_NULL);
    return(IDASLS_LMEM_NULL);
  }
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;

  idasls_mem->s_jaceval = jac;

  return(IDASLS_SUCCESS);
}

/*
 * IDASlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int IDASlsGetNumJacEvals(void *ida_mem, long int *njevals)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDASlsGetNumJacEvals", MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_LMEM_NULL, "IDASSLS", 
		    "IDASlsGetNumJacEvals", MSGSP_LMEM_NULL);
    return(IDASLS_LMEM_NULL);
  }
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;

  *njevals = idasls_mem->s_nje;

  return(IDASLS_SUCCESS);
}

/*
 * IDASlsGetReturnFlagName returns the name associated with a IDASLS
 * return value.
 */
char *IDASlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case IDASLS_SUCCESS:
    sprintf(name,"IDASLS_SUCCESS");
    break;   
  case IDASLS_MEM_NULL:
    sprintf(name,"IDASLS_MEM_NULL");
    break;
  case IDASLS_LMEM_NULL:
    sprintf(name,"IDASLS_LMEM_NULL");
    break;
  case IDASLS_ILL_INPUT:
    sprintf(name,"IDASLS_ILL_INPUT");
    break;
  case IDASLS_MEM_FAIL:
    sprintf(name,"IDASLS_MEM_FAIL");
    break;
  case IDASLS_JAC_NOSET:
    sprintf(name,"IDASLS_JAC_NOSET");
    break;
  case IDASLS_JACFUNC_UNRECVR:
    sprintf(name,"IDASLS_JACFUNC_UNRECVR");
    break;
  case IDASLS_JACFUNC_RECVR:
    sprintf(name,"IDASLS_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * IDASlsGetLastFlag returns the last flag set in a IDASLS function.
 */
int IDASlsGetLastFlag(void *ida_mem, long int *flag)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDASlsGetLastFlag", 
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_LMEM_NULL, "IDASSLS", 
		    "IDASlsGetLastFlag", MSGSP_LMEM_NULL);
    return(IDASLS_LMEM_NULL);
  }
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;

  *flag = idasls_mem->s_last_flag;

  return(IDASLS_SUCCESS);
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

int IDASlsSetSparseJacFnB(void *ida_mem, int which, IDASlsSparseJacFnB jacB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  IDASlsMemB idaslsB_mem;
  void *ida_memB;
  int flag;
  
  /* Is ida_mem allright? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDASlsSetSparseJacFnB", 
		    MSGSP_CAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASLS_NO_ADJ, "IDASSLS", "IDASlsSetSparseJacFnB", 
		    MSGSP_NO_ADJ);
    return(IDASLS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", 
		    "IDASlsSetSparseJacFnB", MSGSP_BAD_WHICH);
    return(IDASLS_ILL_INPUT);
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

  if (IDAB_mem->ida_lmem == NULL) {
    IDAProcessError(IDAB_mem->IDA_mem, IDASLS_LMEMB_NULL, 
                    "IDASSLS", "IDASlsSetSparseJacFnB", MSGSP_LMEMB_NULL);
    return(IDASLS_LMEMB_NULL);
  }
  idaslsB_mem = (IDASlsMemB) IDAB_mem->ida_lmem;

  idaslsB_mem->s_djacB = jacB;

  if (jacB != NULL) {
    flag = IDASlsSetSparseJacFn(ida_memB, idaSlsSparseJacBWrapper);
  } else {
    flag = IDASlsSetSparseJacFn(ida_memB, NULL);
  }

  return(flag);
}

int IDASlsSetSparseJacFnBS(void *ida_mem, int which, IDASlsSparseJacFnBS jacBS)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  IDASlsMemB idaslsB_mem;
  void *ida_memB;
  int flag;
  
  /* Is ida_mem allright? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASDLS", "IDASlsSetSparseJacFnBS", MSGSP_CAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASLS_NO_ADJ, "IDASSLS", "IDASlsSetSparseJacFnBS",  MSGSP_NO_ADJ);
    return(IDASLS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDASlsSetSparseJacFnBS", MSGSP_BAD_WHICH);
    return(IDASLS_ILL_INPUT);
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

  if (IDAB_mem->ida_lmem == NULL) {
    IDAProcessError(IDAB_mem->IDA_mem, IDASLS_LMEMB_NULL, 
                    "IDASSLS", "IDASlsSetSparseJacFnBS", MSGSP_LMEMB_NULL);
    return(IDASLS_LMEMB_NULL);
  }
  idaslsB_mem = (IDASlsMemB) IDAB_mem->ida_lmem;

  idaslsB_mem->s_djacBS = jacBS;

  if (jacBS != NULL) {
    flag = IDASlsSetSparseJacFn(ida_memB, idaSlsSparseJacBSWrapper);
  } else {
    flag = IDASlsSetSparseJacFn(ida_memB, NULL);
  }

  return(flag);
}

/*
 * -----------------------------------------------------------------
 * PRIVATE INTERFACE FUNCTIONS
 * -----------------------------------------------------------------
 */

/*
 * idaSlsSparseJacBWrapper
 *
 * This routine interfaces to the IDASlsSparseJacFnB routine provided 
 * by the user. idaSlsSparseJacBWrapper is of type IDASlsSparseJacFn.
 * NOTE: data actually contains ida_mem
 */

static int idaSlsSparseJacBWrapper(realtype tt, realtype c_jB,
                       N_Vector yyB, N_Vector ypB, N_Vector rrB,
  		       SlsMat JacMat, void *ida_mem, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  IDASlsMemB idaslsB_mem;
  int flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Get current backward problem. */
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
  
  /* Get linear solver's data for this backward problem. */
  idaslsB_mem = (IDASlsMemB) IDAB_mem->ida_lmem;

  /* Forward solution from interpolation */
  if (IDAADJ_mem->ia_noInterp == FALSE) {
    flag = IDAADJ_mem->ia_getY(IDA_mem, tt, IDAADJ_mem->ia_yyTmp, 
			       IDAADJ_mem->ia_ypTmp, NULL, NULL);
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDAB_mem->IDA_mem, -1, "IDASSLS",
		      "idaSlsSparseJacBWrapper", MSGSP_BAD_T);
      return(-1);
    }
  }

  /* Call user's adjoint sparse djacB routine */
  flag = idaslsB_mem->s_djacB(tt, c_jB, 
                              IDAADJ_mem->ia_yyTmp, IDAADJ_mem->ia_ypTmp, 
                              yyB, ypB, rrB, 
                              JacMat, IDAB_mem->ida_user_data,
                              tmp1B, tmp2B, tmp3B);
  return(flag);
}

/*
 * idaSlsSparseJacBSWrapper
 *
 * This routine interfaces to the IDASlsSparseJacFnBS routine provided 
 * by the user. idaSlsSparseJacBSWrapper is of type IDASlsSparseJacFn.
 * NOTE: data actually contains ida_mem
 */

static int idaSlsSparseJacBSWrapper(realtype tt, realtype c_jB,
                       N_Vector yyB, N_Vector ypB, N_Vector rrB,
  		       SlsMat JacMat, void *ida_mem, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDA_mem;
  IDABMem IDAB_mem;
  IDASlsMemB idaslsB_mem;
  int flag;

  IDA_mem = (IDAMem) ida_mem;
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Get current backward problem. */
  IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
  
  /* Get linear solver's data for this backward problem. */
  idaslsB_mem = (IDASlsMemB) IDAB_mem->ida_lmem;

  /* Forward solution from interpolation */
  if (IDAADJ_mem->ia_noInterp == FALSE) {
    if (IDAADJ_mem->ia_interpSensi)
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, IDAADJ_mem->ia_yyTmp, 
				 IDAADJ_mem->ia_ypTmp,
				 IDAADJ_mem->ia_yySTmp,IDAADJ_mem->ia_yySTmp);
    else
      flag = IDAADJ_mem->ia_getY(IDA_mem, tt, IDAADJ_mem->ia_yyTmp, 
			       IDAADJ_mem->ia_ypTmp, NULL, NULL);
    if (flag != IDA_SUCCESS) {
      IDAProcessError(IDAB_mem->IDA_mem, -1, "IDASSLS",
		      "idaSlsSparseJacBSWrapper", MSGSP_BAD_T);
      return(-1);
    }
  }

  /* Call user's adjoint sparse djacB routine */
  flag = idaslsB_mem->s_djacBS(tt, c_jB, 
                              IDAADJ_mem->ia_yyTmp, IDAADJ_mem->ia_ypTmp, 
			      IDAADJ_mem->ia_yySTmp, IDAADJ_mem->ia_ypSTmp,
			      yyB, ypB, rrB, 
                              JacMat, IDAB_mem->ida_user_data,
                              tmp1B, tmp2B, tmp3B);
  return(flag);
}
