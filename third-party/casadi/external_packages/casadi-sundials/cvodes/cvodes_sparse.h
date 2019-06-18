/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer: Carol S. Woodward @ LLNL
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
 * This is the header file for the Sparse linear solver module in CVODES.
 * -----------------------------------------------------------------
 */

#ifndef _CVSSPARSE_H
#define _CVSSPARSE_H

#include <sundials/sundials_sparse.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 * C V S S P A R S E     C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * CVSSPARSE return values 
 * -----------------------------------------------------------------
 */

#define CVSLS_SUCCESS           0
#define CVSLS_MEM_NULL         -1
#define CVSLS_LMEM_NULL        -2
#define CVSLS_ILL_INPUT        -3
#define CVSLS_MEM_FAIL         -4
#define CVSLS_JAC_NOSET        -5
#define CVSLS_PACKAGE_FAIL     -6

/* Additional last_flag values */

#define CVSLS_JACFUNC_UNRECVR  -7
#define CVSLS_JACFUNC_RECVR    -8

/* Return values for the adjoint module */

#define CVSLS_NO_ADJ           -101
#define CVSLS_LMEMB_NULL       -102

/*
 * =================================================================
 * PART I:  F O R W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * FUNCTION TYPES
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Types : CVSlsSparseJacFn
 * -----------------------------------------------------------------
 *
 * A sparse Jacobian approximation function jac must be of type 
 * CVSlsSparseJacFn.
 * Its parameters are:                     
 *                                                                
 * t   is the current value of the independent variable t.        
 *                                                                
 * y   is the current value of the dependent variable vector,     
 *     namely the predicted value of y(t).                     
 *                                                                
 * fy  is the vector f(t,y).
 *     namely the predicted value of y'(t).                    
 *                                                                
 * JacMat is the compressed sparse column matrix (of type SlsMat)
 *     to be loaded by an CVSlsSparseJacFn routine with an approximation
 *     to the system Jacobian matrix
 *            J = J = (df_i/dy_j) at the point (t,y). 
 *     Note that JacMat is NOT preset to zero!
 *     Matrix data is for the nonzero entries of the Jacobian stored in
 *     compressed column format.  Row indices of entries in 
 *     column j are stored in J->rowvals[colptrs[j]] 
 *     through J->rowvals[colptrs[j+i]-1]
 *     and corresponding numerical values of the Jacobian are stored 
 *     in the same entries of J->data.
 * 
 * J_data is a pointer to user Jacobian data - the same as the    
 *     user_data parameter passed to CVodeSetFdata.                     
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          
 *     N_Vectors which can be used by an CVSparseJacFn routine 
 *     as temporary storage or work space.                     
 *                                                                
 * A CVSlsSparseJacFn should return                                
 *     0 if successful,                                           
 *     a positive int if a recoverable error occurred, or         
 *     a negative int if a nonrecoverable error occurred.         
 *
 * -----------------------------------------------------------------
 *
  * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively 
 *     (see cvode.h). The unit roundoff is available as 
 *     UNIT_ROUNDOFF defined in sundials_types.h.
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*CVSlsSparseJacFn)(realtype t,
		     N_Vector y, N_Vector fy, 
		     SlsMat JacMat, void *user_data,
		     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVSPARSE linear solver
 * -----------------------------------------------------------------
 * CVSlsSetSparseJacFn specifies the Jacobian approximation
 * routine to be used for a sparse direct linear solver.
 *
 * The return value is one of:
 *    CVSLS_SUCCESS   if successful
 *    CVSLS_MEM_NULL  if the CVODE memory was NULL
 *    CVSLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVSlsSetSparseJacFn(void *cvode_mem, CVSlsSparseJacFn jac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVSLS linear solver
 * -----------------------------------------------------------------
 *
 * CVSlsGetNumJacEvals returns the number of calls made to the
 *                      Jacobian evaluation routine jac.
 * CVSlsGetLastFlag    returns the last error flag set by any of
 *                      the CVSLS interface functions.
 *
 * The return value of CVSlsGet* is one of:
 *    CVSLS_SUCCESS   if successful
 *    CVSLS_MEM_NULL  if the IDA memory was NULL
 *    CVSLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVSlsGetNumJacEvals(void *cvode_mem, long int *njevals);
SUNDIALS_EXPORT int CVSlsGetLastFlag(void *cvode_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a CVSLS return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *CVSlsGetReturnFlagName(long int flag);

/*
 * =================================================================
 * PART II:  B A C K W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * FUNCTION TYPES
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Type: CVSlsSparseJacFnB
 * -----------------------------------------------------------------
 * A sparse Jacobian approximation function jacB for the adjoint
 * (backward) problem must have the prototype given below. 
 * -----------------------------------------------------------------
 */

typedef int (*CVSlsSparseJacFnB)(realtype t, N_Vector y, 
				 N_Vector yB, N_Vector fyB,
				 SlsMat JB, void *user_dataB, 
				 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * Type: CVSlsSparseJacFnBS
 * -----------------------------------------------------------------
 * A sparse Jacobian approximation function jacBS for the adjoint
 * (backward) problem, sensitivity-dependent case,  must have the
 *  prototype given below. 
 * -----------------------------------------------------------------
 */

typedef int (*CVSlsSparseJacFnBS)(realtype t,
				  N_Vector y, N_Vector *yS,
				  N_Vector yB, N_Vector fyB,
				  SlsMat JB, void *user_dataB, 
				  N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * EXPORTED FUNCTIONS 
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Functions: CVSlsSetSparseJacFnB and CVSlsSetSparseJacFnBS
 * -----------------------------------------------------------------
 * CVSlsSetSparseJacFnB specifies the sparse Jacobian functions to 
 * be used by a CVSPARSE linear solver for the backward integration phase
 * when the backward problem does not depend on forward sensitivities.
 * CVSlsSetSparseJacFnBS specifies the Jacobian
 * functions when the backward problem does depend on sensitivities.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVSlsSetSparseJacFnB(void *cv_mem, int which, 
					 CVSlsSparseJacFnB jacB);
SUNDIALS_EXPORT int CVSlsSetSparseJacFnBS(void *cv_mem, int which, 
					  CVSlsSparseJacFnBS jacBS);

#ifdef __cplusplus
}
#endif

#endif
