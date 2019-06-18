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
 * This is the implementation file for an KINSLS linear solver.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include "kinsol_impl.h"
#include "kinsol_sparse_impl.h"
#include "sundials/sundials_math.h"

#include <stdio.h>
#include <stdlib.h>

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
 * EXPORTED FUNCTIONS 
 * =================================================================
 */
              
/*
 * KINSlsSetSparseJacFn specifies the sparse Jacobian function.
 */
int KINSlsSetSparseJacFn(void *kin_mem_v, KINSlsSparseJacFn jac)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;

  kin_mem = (KINMem) kin_mem_v;

  /* Return immediately if kin_mem is NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "KINSlsSetSparseJacFn", 
		    MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }

  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSLS_LMEM_NULL, "KINSLS", 
		    "KINSlsSetSparseJacFn", MSGSP_LMEM_NULL);
    return(KINSLS_LMEM_NULL);
  }
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;

  kinsls_mem->s_jaceval = jac;

  return(KINSLS_SUCCESS);
}

/*
 * KINSlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int KINSlsGetNumJacEvals(void *kin_mem_v, long int *njevals)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;

  kin_mem = (KINMem) kin_mem_v;

  /* Return immediately if kin_mem is NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "KINSlsGetNumJacEvals", MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }

  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSLS_LMEM_NULL, "KINSLS", 
		    "KINSlsGetNumJacEvals", MSGSP_LMEM_NULL);
    return(KINSLS_LMEM_NULL);
  }
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;

  *njevals = kinsls_mem->s_nje;

  return(KINSLS_SUCCESS);
}

/*
 * KINSlsGetReturnFlagName returns the name associated with a KINSLS
 * return value.
 */
char *KINSlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case KINSLS_SUCCESS:
    sprintf(name,"KINSLS_SUCCESS");
    break;   
  case KINSLS_MEM_NULL:
    sprintf(name,"KINSLS_MEM_NULL");
    break;
  case KINSLS_LMEM_NULL:
    sprintf(name,"KINSLS_LMEM_NULL");
    break;
  case KINSLS_ILL_INPUT:
    sprintf(name,"KINSLS_ILL_INPUT");
    break;
  case KINSLS_MEM_FAIL:
    sprintf(name,"KINSLS_MEM_FAIL");
    break;
  case KINSLS_JAC_NOSET:
    sprintf(name,"KINSLS_JAC_NOSET");
    break;
  case KINSLS_JACFUNC_UNRECVR:
    sprintf(name,"KINSLS_JACFUNC_UNRECVR");
    break;
  case KINSLS_JACFUNC_RECVR:
    sprintf(name,"KINSLS_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * KINSlsGetLastFlag returns the last flag set in a KINSLS function.
 */
int KINSlsGetLastFlag(void *kin_mem_v, long int *flag)
{
  KINMem kin_mem;
  KINSlsMem kinsls_mem;

  kin_mem= (KINMem) kin_mem_v;

  /* Return immediately if kin_mem is NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINSLS_MEM_NULL, "KINSLS", "KINSlsGetLastFlag", 
		    MSGSP_KINMEM_NULL);
    return(KINSLS_MEM_NULL);
  }

  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSLS_LMEM_NULL, "KINSLS", 
		    "KINSlsGetLastFlag", MSGSP_LMEM_NULL);
    return(KINSLS_LMEM_NULL);
  }
  kinsls_mem = (KINSlsMem) kin_mem->kin_lmem;

  *flag = kinsls_mem->s_last_flag;

  return(KINSLS_SUCCESS);
}

