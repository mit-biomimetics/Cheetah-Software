/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
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
 * This is the header file for the Sparse linear solver module in IDAS.
 * -----------------------------------------------------------------
 */

#ifndef _IDASSPARSE_H
#define _IDASSPARSE_H

#include <sundials/sundials_sparse.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 * I D A S S P A R S E    C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * IDASSLS return values 
 * -----------------------------------------------------------------
 */

#define IDASLS_SUCCESS           0
#define IDASLS_MEM_NULL         -1
#define IDASLS_LMEM_NULL        -2
#define IDASLS_ILL_INPUT        -3
#define IDASLS_MEM_FAIL         -4
#define IDASLS_JAC_NOSET        -5
#define IDASLS_PACKAGE_FAIL     -6

/* Additional last_flag values */

#define IDASLS_JACFUNC_UNRECVR  -7
#define IDASLS_JACFUNC_RECVR    -8

/* Return values for the adjoint module */
#define IDASLS_NO_ADJ           -101
#define IDASLS_LMEMB_NULL       -102

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
 * Types : IDASlsSparseJacFn
 * -----------------------------------------------------------------
 *
 * A sparse Jacobian approximation function jaceval must be of type 
 * IDASlsSparseJacFn.
 * Its parameters are:                     
 *                                                                
 * t   is the current value of the independent variable t.        
 *                                                                
 * c_j is the scalar in the system Jacobian, proportional to 
 *     the inverse of the step size h.
 *                                                                
 * y   is the current value of the dependent variable vector,     
 *     namely the predicted value of y(t).                     
 *                                                                
 * yp  is the current value of the derivative vector y',          
 *     namely the predicted value of y'(t).                    
 *                                                                
 * f   is the residual vector F(tt,yy,yp).                     
 *                                                                
 * JacMat is the compressed sparse column matrix (of type SlsMat)
 *     to be loaded by an IDASlsSparseJacFn routine with an approximation
 *     to the system Jacobian matrix
 *            J = dF/dy' + gamma*dF/dy                            
 *     at the given point (t,y,y'), where the ODE system is    
 *     given by F(t,y,y') = 0.
 *     Note that JacMat is NOT preset to zero!
 *     Matrix data is for the nonzero entries of the Jacobian are stored in
 *     compressed column format.  Row indices of entries in 
 *     column j are stored in J->data[colptrs[j]] 
 *     through J->data[colptrs[j+i]-1]
 *     and corresponding numerical values of the Jacobian are stored 
 *     in the same entries of a.
 * 
 * user_data is a pointer to user Jacobian data - the same as the    
 *     user_data parameter passed to IDASetRdata.                     
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          
 *     N_Vectors which can be used by an IDASparseJacFn routine 
 *     as temporary storage or work space.                     
 *                                                                
 * A IDASlsSparseJacFn should return                                
 *     0 if successful,                                           
 *     a positive int if a recoverable error occurred, or         
 *     a negative int if a nonrecoverable error occurred.         
 * In the case of a recoverable error return, the integrator will 
 * attempt to recover by reducing the stepsize (which changes cj).
 *
 * -----------------------------------------------------------------
 *
  * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     IDAGetCurrentStep and IDAGetErrWeights, respectively 
 *     (see ida.h). The unit roundoff is available as 
 *     UNIT_ROUNDOFF defined in sundials_types.h.
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*IDASlsSparseJacFn)(realtype t, realtype c_j,
		     N_Vector y, N_Vector yp, N_Vector r, 
		     SlsMat JacMat, void *user_data,
		     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the IDASSPARSE linear solver
 * -----------------------------------------------------------------
 * IDASlsSetSparseJacFn specifies the Jacobian approximation
 * routine to be used for a sparse direct linear solver.
 *
 * The return value is one of:
 *    IDASLS_SUCCESS   if successful
 *    IDASLS_MEM_NULL  if the IDA memory was NULL
 *    IDASLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASlsSetSparseJacFn(void *ida_mem, IDASlsSparseJacFn jac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the IDASLS linear solver
 * -----------------------------------------------------------------
 *
 * IDASlsGetWorkSpace   returns the real and integer workspace used
 *                      by the direct linear solver.
 * IDASlsGetNumJacEvals returns the number of calls made to the
 *                      Jacobian evaluation routine jac.
 * IDASlsGetLastFlag    returns the last error flag set by any of
 *                      the IDADLS interface functions.
 *
 * The return value of IDADlsGet* is one of:
 *    IDASLS_SUCCESS   if successful
 *    IDASLS_MEM_NULL  if the IDA memory was NULL
 *    IDASLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASlsGetNumJacEvals(void *ida_mem, long int *njevals);
SUNDIALS_EXPORT int IDASlsGetLastFlag(void *ida_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a IDASLS return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *IDASlsGetReturnFlagName(long int flag);

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
 * Type: IDASlsSparseJacFnB
 * -----------------------------------------------------------------
 * A sparse Jacobian approximation function JacB for the adjoint
 * (backward) problem must have the prototype given below. 
 * -----------------------------------------------------------------
 */

typedef int (*IDASlsSparseJacFnB)(realtype tt, realtype c_jB, 
	      	     N_Vector yy, N_Vector yp,
		     N_Vector yyB, N_Vector ypB, N_Vector rrB, 
		     SlsMat JacMatB, void *user_dataB, 
		     N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);


/*
 * -----------------------------------------------------------------
 * Type: IDASlsSparseJacFnBS
 * -----------------------------------------------------------------
 * A dense Jacobian approximation function JacBS for the adjoint
 * (backward) problem, sensitivity-dependent case, must have the
 *  prototype given below. 
 * -----------------------------------------------------------------
 */

typedef int (*IDASlsSparseJacFnBS)(realtype tt, realtype c_jB, 
	      	     N_Vector yy, N_Vector yp,
		     N_Vector *yS, N_Vector *ypS,
		     N_Vector yyB, N_Vector ypB, N_Vector rrB, 
		     SlsMat JacMatB, void *user_dataB, 
		     N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);


/*
 * -----------------------------------------------------------------
 * EXPORTED FUNCTIONS 
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Functions: IDASlsSetSparseJacFnB and IDASlsSetSparseJacFnBS
 * -----------------------------------------------------------------
 * IDASlsSetSparseJacFnB specifies the sparse Jacobian functions to 
 * be used by a IDASSPARSE linear solver for the backward integration phase
 * when the backward problem does not depend on forward sensitivities.
 * IDASlsSetSparseJacFnBS specifies the Jacobian
 * functions when the backward problem does depend on sensitivities.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASlsSetSparseJacFnB(void *ida_mem, int which, 
					  IDASlsSparseJacFnB jacB);
SUNDIALS_EXPORT int IDASlsSetSparseJacFnBS(void *ida_mem, int which, 
					  IDASlsSparseJacFnBS jacBS);


#ifdef __cplusplus
}
#endif

#endif
