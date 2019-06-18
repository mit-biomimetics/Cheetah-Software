/*
 * -----------------------------------------------------------------
 * $Revision: 4268 $
 * $Date: 2014-11-14 16:37:06 -0800 (Fri, 14 Nov 2014) $
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
 * This is the header file for the Sparse linear solver module in KINSOL.
 * -----------------------------------------------------------------
 */

#ifndef _KINSOLSPARSE_H
#define _KINSOLSPARSE_H

#include "sundials/sundials_sparse.h"
#include "sundials/sundials_nvector.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 * K I N S O L S P A R S E    C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * KINSLS return values 
 * -----------------------------------------------------------------
 */

#define KINSLS_SUCCESS           0
#define KINSLS_MEM_NULL         -1
#define KINSLS_LMEM_NULL        -2
#define KINSLS_ILL_INPUT        -3
#define KINSLS_MEM_FAIL         -4
#define KINSLS_JAC_NOSET        -5
#define KINSLS_PACKAGE_FAIL     -6

/* Additional last_flag values */

#define KINSLS_JACFUNC_UNRECVR  -7
#define KINSLS_JACFUNC_RECVR    -8

/*
 * -----------------------------------------------------------------
 * FUNCTION TYPES
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Types : KINSlsSparseJacFn
 * -----------------------------------------------------------------
 *
 * A sparse Jacobian approximation function jaceval must be of type 
 * KINSlsSparseJacFn.
 * Its parameters are:                     
 *                                                                
 * u   is the current value of the dependent variable vector     
 *                                                                
 * fval is the residual vector F(tt,yy,yp).                     
 *                                                                
 * JacMat is the compressed sparse column matrix (of type SlsMat)
 *     to be loaded by an KINSlsSparseJacFn routine with an approximation
 *     to the system Jacobian matrix
 *            J = dF/dy
 *     Note that JacMat is NOT preset to zero!
 *     Matrix data is for the nonzero entries of the Jacobian are stored in
 *     compressed column format.  Row indices of entries in 
 *     column j are stored in J->data[colptrs[j]] 
 *     through J->data[colptrs[j+i]-1]
 *     and corresponding numerical values of the Jacobian are stored 
 *     in the same entries of a.
 * 
 * user_data is a pointer to user Jacobian data 
 *                                                                
 * vtemp1 and vtemp2 are pointers to memory allocated for          
 *     N_Vectors which can be used by an KINSparseJacFn routine 
 *     as temporary storage or work space.                     
 *                                                                
 * A KINSlsSparseJacFn should return                                
 *     0 if successful,                      
 *     1 if successful but the user code indicates a new factorization
 *       of the Jacobian is required
 *     a positive int if a recoverable error occurred, or         
 *     a negative int if a nonrecoverable error occurred.         
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*KINSlsSparseJacFn)(N_Vector u, N_Vector fval, 
		     SlsMat JacMat, void *user_data,
		     N_Vector vtemp1, N_Vector vtemp2);

/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the KINSPARSE linear solver
 * -----------------------------------------------------------------
 * KINSlsSetSparseJacFn specifies the Jacobian approximation
 * routine to be used for a sparse direct linear solver.
 *
 * The return value is one of:
 *    KINSLS_SUCCESS   if successful
 *    KINSLS_MEM_NULL  if the IDA memory was NULL
 *    KINSLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSlsSetSparseJacFn(void *kin_mem_v, KINSlsSparseJacFn jac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the KINSLS linear solver
 * -----------------------------------------------------------------
 *
 * KINSlsGetNumJacEvals returns the number of calls made to the
 *                      Jacobian evaluation routine jac.
 * KINSlsGetLastFlag    returns the last error flag set by any of
 *                      the KINSLS interface functions.
 *
 * The return value of KINSlsGet* is one of:
 *    KINSLS_SUCCESS   if successful
 *    KINSLS_MEM_NULL  if the KINSOL memory was NULL
 *    KINSLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSlsGetNumJacEvals(void *kin_mem, long int *njevals);
SUNDIALS_EXPORT int KINSlsGetLastFlag(void *kin_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a IDASLS return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *KINSlsGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
